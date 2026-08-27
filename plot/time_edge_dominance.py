"""
time_edge_dominance.py
----------------------
Plot how edge / interior fractions of chained vs unchained bacteria evolve
over time, for each growth factor.

Directory structure expected:
    <base>/growthRateMultiplier1_0p7/Ch50_unch50/repeat0/biofilm_00010.dat
                                                          biofilm_00020.dat
                                                          ...
    <base>/growthRateMultiplier1_0p7/Ch50_unch50/repeat1/...

Time: biofilm_NNNNN.dat  →  t = NNNNN × dt  (default dt = 0.1 hrs)

Usage
-----
    python time_edge_dominance.py --data ../output/data_production/
    python time_edge_dominance.py --data ../output/data_production/ --dt 0.1
    python time_edge_dominance.py --data ../output/data_production/ --condition Ch50_unch50
    python time_edge_dominance.py --data ../output/data_production/ --save

Options
-------
    --data              Base directory containing growthRateMultiplier* subdirs
    --condition         Condition subdir name (default: Ch50_unch50)
    --dt                Time per file step in hours (default: 0.1)
    --edge-threshold    Edge band width in multiples of mean radius (default: 2.0)
    --save              Save plots as PNG instead of opening a window

Requirements
------------
    pip install matplotlib numpy pandas scipy tqdm
"""

import re
import sys
import argparse
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.spatial import ConvexHull

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Colour palette ─────────────────────────────────────────────────────────────
C_CHAINED   = "#0072B2"
C_UNCHAINED = "#E69F00"
DARK_BG     = "#0a0d0f"
PANEL_BG    = "#0f1519"
GRID_COL    = "#ffffff12"
TEXT_COL    = "#d0dde8"
SPINE       = "#ffffff25"


# ── Helpers ────────────────────────────────────────────────────────────────────
def parse_growth_factor(path):
    for part in Path(path).parts:
        m = re.search(r"growthRateMultiplier\d+_(\d+)p(\d+)", part)
        if m:
            return float(f"{m.group(1)}.{m.group(2)}")
    raise ValueError(f"No growthRateMultiplier*_NpM pattern in: {path}")


def parse_timestep(filename):
    """biofilm_00054.dat  →  54"""
    m = re.search(r"biofilm_(\d+)\.dat", Path(filename).name)
    if m:
        return int(m.group(1))
    raise ValueError(f"Cannot parse timestep from: {filename}")


def load_dat(path):
    df = pd.read_csv(
        path, sep=r"\s+",
        na_values=["None", "-", ""],
        dtype={"lower_link": "Int64", "upper_link": "Int64"},
    )
    required = {"cell_id", "length", "radius", "pos_x", "pos_y",
                "lower_link", "upper_link"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    df["_chained"] = (df["lower_link"].notna() | df["upper_link"].notna()).astype(int)
    return df


# ── Edge detection ─────────────────────────────────────────────────────────────
def dist_to_hull(points, hull):
    dists = np.full(len(points), np.inf)
    verts = hull.points[hull.vertices]
    n = len(verts)
    for i in range(n):
        A, B = verts[i], verts[(i + 1) % n]
        AB = B - A
        AB_len2 = np.dot(AB, AB)
        if AB_len2 == 0:
            d = np.linalg.norm(points - A, axis=1)
        else:
            t    = np.clip(((points - A) @ AB) / AB_len2, 0, 1)
            proj = A + np.outer(t, AB)
            d    = np.linalg.norm(points - proj, axis=1)
        dists = np.minimum(dists, d)
    return dists


def compute_metrics(df, edge_thresh_mult=2.0):
    xy      = df[["pos_x", "pos_y"]].values
    mean_r  = df["radius"].mean()
    thresh  = mean_r * edge_thresh_mult
    n_total = len(df)

    if n_total < 3:
        return None
    try:
        hull = ConvexHull(xy)
    except Exception:
        return None

    d         = dist_to_hull(xy, hull)
    edge_mask = d <= thresh
    int_mask  = ~edge_mask

    def fracs(mask):
        n = mask.sum()
        if n == 0:
            return 0.0, 0.0, 0
        nc = df.loc[mask, "_chained"].sum()
        return nc / n, 1 - nc / n, int(n)

    fc_edge, fu_edge, n_edge = fracs(edge_mask)
    fc_int,  fu_int,  n_int  = fracs(int_mask)
    fc_all  = df["_chained"].sum() / n_total
    fu_all  = 1 - fc_all

    return {
        "n_total":           n_total,
        "n_edge":            n_edge,
        "n_interior":        n_int,
        "frac_edge_chained":      fc_edge,
        "frac_edge_unchained":    fu_edge,
        "frac_int_chained":       fc_int,
        "frac_int_unchained":     fu_int,
        "frac_all_chained":       fc_all,
        "frac_all_unchained":     fu_all,
    }


# ── Discovery ──────────────────────────────────────────────────────────────────
def discover(base_dir, condition):
    """
    Returns:
        data[gf][timestep] = list of metric dicts (one per repeat)
    """
    base  = Path(base_dir)
    gf_dirs = sorted(base.glob("growthRateMultiplier*"))
    if not gf_dirs:
        raise FileNotFoundError(f"No growthRateMultiplier* dirs under {base}")

    data = defaultdict(lambda: defaultdict(list))

    for gf_dir in gf_dirs:
        try:
            gf = parse_growth_factor(str(gf_dir))
        except ValueError:
            continue

        cond_dir = gf_dir / condition
        if not cond_dir.is_dir():
            print(f"  Warning: {cond_dir} not found, skipping.")
            continue

        repeat_dirs = sorted(cond_dir.glob("repeat*"))
        for rep_dir in repeat_dirs:
            files = sorted(rep_dir.glob("biofilm_*.dat"))
            if not files:
                continue
            print(f"  gf={gf:.2f}  {rep_dir.name}  →  {len(files)} files "
                  f"({files[0].name} … {files[-1].name})")
            for fp in files:
                try:
                    ts = parse_timestep(fp)
                    df = load_dat(fp)
                    m  = compute_metrics(df)
                    if m:
                        data[gf][ts].append(m)
                except Exception as e:
                    print(f"    Skipping {fp.name}: {e}")

    return data


# ── Averaging over repeats ─────────────────────────────────────────────────────
def average_over_repeats(metric_list):
    keys = metric_list[0].keys()
    return {k: np.mean([m[k] for m in metric_list]) for k in keys}


def std_over_repeats(metric_list):
    keys = metric_list[0].keys()
    return {k: np.std([m[k] for m in metric_list]) for k in keys}


def build_timeseries(data, dt):
    """
    Returns dict:
        ts_data[gf] = {
            't':   array of times (hrs),
            'mean': {metric: array},
            'std':  {metric: array},
        }
    """
    ts_data = {}
    for gf in sorted(data.keys()):
        ts_steps = sorted(data[gf].keys())
        t_arr  = np.array(ts_steps) * dt

        mean_rows, std_rows = [], []
        for ts in ts_steps:
            reps = data[gf][ts]
            mean_rows.append(average_over_repeats(reps))
            std_rows.append(std_over_repeats(reps) if len(reps) > 1
                            else {k: 0.0 for k in reps[0]})

        keys = mean_rows[0].keys()
        ts_data[gf] = {
            "t":    t_arr,
            "mean": {k: np.array([r[k] for r in mean_rows]) for k in keys},
            "std":  {k: np.array([r[k] for r in std_rows])  for k in keys},
        }
    return ts_data


# ── Styling ────────────────────────────────────────────────────────────────────
def style_ax(ax, xlabel="", ylabel="", title=""):
    ax.set_facecolor(PANEL_BG)
    ax.grid(color=GRID_COL, linewidth=0.6, zorder=0)
    for sp in ax.spines.values():
        sp.set_edgecolor(SPINE)
    ax.tick_params(colors=TEXT_COL, labelsize=8)
    ax.xaxis.label.set_color(TEXT_COL)
    ax.yaxis.label.set_color(TEXT_COL)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, color=TEXT_COL, fontsize=9, pad=6)


# ── Plotting ───────────────────────────────────────────────────────────────────
def make_plots(ts_data, save_path=None):
    gf_list  = sorted(ts_data.keys())
    n_gf     = len(gf_list)

    # Colormap: one colour per growth factor
    cmap   = cm.get_cmap("plasma", n_gf)
    colours = {gf: cmap(i) for i, gf in enumerate(gf_list)}

    # ── Figure 1: Edge & interior fraction over time ─────────────────────────
    fig1, axes = plt.subplots(2, 2, figsize=(13, 8), facecolor=DARK_BG)
    fig1.suptitle(
        "Edge & Interior Fractions over Time — by Growth Factor",
        color=TEXT_COL, fontsize=12, y=0.99,
    )

    panels = [
        (axes[0, 0], "frac_edge_chained",   "Fraction of edge cells",     "Edge — Chained fraction"),
        (axes[0, 1], "frac_edge_unchained",  "Fraction of edge cells",     "Edge — Unchained fraction"),
        (axes[1, 0], "frac_int_chained",     "Fraction of interior cells", "Interior — Chained fraction"),
        (axes[1, 1], "frac_int_unchained",   "Fraction of interior cells", "Interior — Unchained fraction"),
    ]

    for ax, metric, ylabel, title in panels:
        for gf in gf_list:
            d    = ts_data[gf]
            t    = d["t"]
            mean = d["mean"][metric]
            std  = d["std"][metric]
            col  = colours[gf]
            ax.plot(t, mean, lw=1.8, color=col, label=f"GF={gf:.2f}")
            ax.fill_between(t, mean - std, mean + std, color=col, alpha=0.15)
        ax.axhline(0.5, color="white", lw=0.6, ls=":", alpha=0.4)
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
        style_ax(ax, "Time (hrs)", ylabel, title)

    # Shared legend via colorbar
    sm = cm.ScalarMappable(cmap=cmap,
                           norm=mcolors.Normalize(vmin=min(gf_list),
                                                  vmax=max(gf_list)))
    sm.set_array([])
    cbar = fig1.colorbar(sm, ax=axes, orientation="vertical",
                         fraction=0.02, pad=0.02)
    cbar.set_label("Growth factor", color=TEXT_COL, fontsize=9)
    cbar.ax.yaxis.set_tick_params(color=TEXT_COL, labelsize=7)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=TEXT_COL)

    plt.tight_layout(rect=[0, 0, 0.97, 0.97])

    if save_path:
        p = Path(save_path)
        out1 = p.parent / (p.stem + "_timeseries.png")
        fig1.savefig(out1, dpi=150, facecolor=fig1.get_facecolor(),
                     bbox_inches="tight")
        print(f"Saved → {out1}")
    else:
        plt.show()
    plt.close(fig1)

    # ── Figure 2: Heatmaps (GF × time) ──────────────────────────────────────
    # Build common time grid by interpolating each GF series
    all_times = sorted({t for gf in gf_list for t in ts_data[gf]["t"]})
    t_grid    = np.array(all_times)
    gf_arr    = np.array(gf_list)

    heatmap_metrics = [
        ("frac_edge_chained",   "Edge — Chained fraction",   "Blues"),
        ("frac_edge_unchained", "Edge — Unchained fraction", "Oranges"),
        ("frac_int_chained",    "Interior — Chained fraction",   "Blues"),
        ("frac_int_unchained",  "Interior — Unchained fraction", "Oranges"),
    ]

    fig2, axes2 = plt.subplots(2, 2, figsize=(13, 7), facecolor=DARK_BG)
    fig2.suptitle(
        "Heatmaps: Fraction over Time × Growth Factor",
        color=TEXT_COL, fontsize=12, y=0.99,
    )

    for ax, (metric, title, cmap_name) in zip(axes2.flat, heatmap_metrics):
        Z = np.full((len(gf_arr), len(t_grid)), np.nan)
        for i, gf in enumerate(gf_list):
            d = ts_data[gf]
            Z[i] = np.interp(t_grid, d["t"], d["mean"][metric],
                             left=np.nan, right=np.nan)

        im = ax.imshow(
            Z, aspect="auto", origin="lower",
            extent=[t_grid[0], t_grid[-1], 0, len(gf_arr)],
            vmin=0, vmax=1,
            cmap=cmap_name,
        )
        ax.set_yticks(np.arange(len(gf_arr)) + 0.5)
        ax.set_yticklabels([f"{g:.2f}" for g in gf_arr],
                           color=TEXT_COL, fontsize=7)
        ax.tick_params(colors=TEXT_COL, labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(SPINE)
        ax.set_facecolor(PANEL_BG)
        ax.set_xlabel("Time (hrs)", color=TEXT_COL, fontsize=9)
        ax.set_ylabel("Growth factor", color=TEXT_COL, fontsize=9)
        ax.set_title(title, color=TEXT_COL, fontsize=9, pad=6)

        cb = fig2.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
        cb.ax.yaxis.set_tick_params(color=TEXT_COL, labelsize=7)
        cb.set_label("Fraction", color=TEXT_COL, fontsize=8)
        plt.setp(cb.ax.yaxis.get_ticklabels(), color=TEXT_COL)
        cb.ax.yaxis.set_major_formatter(
            ticker.PercentFormatter(xmax=1, decimals=0)
        )

    plt.tight_layout(rect=[0, 0, 1, 0.97])

    if save_path:
        p = Path(save_path)
        out2 = p.parent / (p.stem + "_heatmap.png")
        fig2.savefig(out2, dpi=150, facecolor=fig2.get_facecolor(),
                     bbox_inches="tight")
        print(f"Saved → {out2}")
    else:
        plt.show()
    plt.close(fig2)


# ── Summary table ──────────────────────────────────────────────────────────────
def print_summary(ts_data):
    print(f"\n{'GF':>5}  {'Time pts':>8}  {'T_start':>8}  {'T_end':>7}  "
          f"{'Edge_C @end':>12}  {'Edge_U @end':>12}  {'Edge winner':>12}")
    print("─" * 75)
    for gf in sorted(ts_data.keys()):
        d  = ts_data[gf]
        t  = d["t"]
        ec = d["mean"]["frac_edge_chained"]
        eu = d["mean"]["frac_edge_unchained"]
        winner = "CHAINED" if ec[-1] > eu[-1] else \
                 "UNCHAINED" if eu[-1] > ec[-1] else "TIE"
        print(f"{gf:>5.2f}  {len(t):>8}  {t[0]:>8.2f}  {t[-1]:>7.2f}  "
              f"{ec[-1]:>12.1%}  {eu[-1]:>12.1%}  {winner:>12}")
    print()


# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Plot edge/interior fractions vs time for each growth factor.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--data", required=True, metavar="DIR",
                        help="Base directory with growthRateMultiplier* subdirs")
    parser.add_argument("--condition", default="Ch50_unch50", metavar="NAME",
                        help="Condition subdir name (default: Ch50_unch50)")
    parser.add_argument("--dt", type=float, default=0.1, metavar="HRS",
                        help="Time interval per file step in hours (default: 0.1)")
    parser.add_argument("--edge-threshold", type=float, default=2.0, metavar="N",
                        help="Edge band = N × mean_radius from hull (default: 2.0)")
    parser.add_argument("--save", action="store_true",
                        help="Save plots as PNG instead of opening windows")
    args = parser.parse_args()

    print(f"\nScanning : {Path(args.data).resolve()}")
    print(f"Condition: {args.condition}")
    print(f"dt       : {args.dt} hrs per file step")
    print(f"Edge band: {args.edge_threshold} × mean_radius\n")

    try:
        raw = discover(Path(args.data), args.condition)
    except FileNotFoundError as e:
        print(f"Error: {e}"); sys.exit(1)

    if not raw:
        print("No usable data found."); sys.exit(1)

    print(f"\nBuilding time series …")
    ts_data = build_timeseries(raw, args.dt)

    print_summary(ts_data)

    save_path = Path(args.data) / "edge_dominance" if args.save else None
    if save_path:
        save_path.mkdir(exist_ok=True)
        save_path = save_path / "plot"

    make_plots(ts_data, save_path=save_path)


if __name__ == "__main__":
    main()
