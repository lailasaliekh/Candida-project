"""
edge_dominance.py
-----------------
Quantify which cell type (chained vs unchained) dominates the EDGE of the
colony across growth factor values encoded in the directory path.

Growth factor is parsed from path components like:
    growthRateMultiplier1_0p7  →  0.7
    growthRateMultiplier1_0p0  →  0.0

"Edge cells" are defined as cells whose centre lies within a distance
threshold of the colony's convex hull boundary (default: 2 × mean cell radius).

Usage
-----
    # One file per growth factor:
    python edge_dominance.py \\
        /data/growthRateMultiplier1_0p0/Ch50_unch50/repeat0/biofilm_00054.dat \\
        /data/growthRateMultiplier1_0p2/Ch50_unch50/repeat0/biofilm_00054.dat \\
        ...

    # With multiple repeats — pass all files, the script groups & averages:
    python edge_dominance.py /data/growthRateMultiplier1_*/Ch50_unch50/repeat*/biofilm_00054.dat

    # Save output:
    python edge_dominance.py ... --save

    # Custom edge threshold (multiples of mean radius, default=2):
    python edge_dominance.py ... --edge-threshold 3

Requirements
------------
    pip install matplotlib numpy pandas scipy
"""

import re
import sys
import glob
import argparse
import warnings
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Colour palette (Wong – colorblind safe) ────────────────────────────────────
C_CHAINED   = "#0072B2"   # blue
C_UNCHAINED = "#E69F00"   # orange
C_EDGE_C    = "#56B4E9"
C_EDGE_U    = "#F5C842"
C_TOTAL     = "#aaaaaa"

DARK_BG  = "#0a0d0f"
PANEL_BG = "#0f1519"
GRID_COL = "#ffffff12"
TEXT_COL = "#d0dde8"
SPINE    = "#ffffff25"


# ── Path parser ────────────────────────────────────────────────────────────────
def parse_growth_factor(path):
    """
    Extract growth factor from a path containing e.g. growthRateMultiplier1_0p7
    Returns float or raises ValueError.
    """
    parts = Path(path).parts
    for part in parts:
        m = re.search(r"growthRateMultiplier\d+_(\d+)p(\d+)", part)
        if m:
            return float(f"{m.group(1)}.{m.group(2)}")
    raise ValueError(
        f"Could not find growthRateMultiplier*_NpM pattern in path:\n  {path}"
    )


# ── Data loading ───────────────────────────────────────────────────────────────
def load_dat(path):
    df = pd.read_csv(
        path,
        sep=r"\s+",
        na_values=["None", "-", ""],
        dtype={"lower_link": "Int64", "upper_link": "Int64"},
    )
    required = {"cell_id", "length", "radius", "pos_x", "pos_y",
                "lower_link", "upper_link"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    df["_chained"] = (
        df["lower_link"].notna() | df["upper_link"].notna()
    ).astype(int)
    return df


# ── Edge detection ─────────────────────────────────────────────────────────────
def dist_to_hull(points, hull):
    """
    Minimum distance from each point to any edge of a 2-D convex hull.
    Returns array of shape (n_points,).
    """
    dists = np.full(len(points), np.inf)
    verts = hull.points[hull.vertices]
    n = len(verts)
    for i in range(n):
        A = verts[i]
        B = verts[(i + 1) % n]
        AB = B - A
        AB_len2 = np.dot(AB, AB)
        if AB_len2 == 0:
            d = np.linalg.norm(points - A, axis=1)
        else:
            t = np.clip(((points - A) @ AB) / AB_len2, 0, 1)
            proj = A + np.outer(t, AB)
            d = np.linalg.norm(points - proj, axis=1)
        dists = np.minimum(dists, d)
    return dists


def compute_edge_metrics(df, edge_threshold_multiplier=2.0):
    """
    Returns dict of metrics.  Edge cells = cells within
    (mean_radius * edge_threshold_multiplier) of the convex hull boundary.
    """
    xy = df[["pos_x", "pos_y"]].values
    mean_r = df["radius"].mean()
    threshold = mean_r * edge_threshold_multiplier

    n_total    = len(df)
    n_chained  = df["_chained"].sum()
    n_unchained = n_total - n_chained

    # Colony-wide fractions
    frac_chained_all   = n_chained   / n_total if n_total else 0
    frac_unchained_all = n_unchained / n_total if n_total else 0

    # Convex hull & edge mask
    if len(xy) < 3:
        return None
    try:
        hull = ConvexHull(xy)
    except Exception:
        return None

    d = dist_to_hull(xy, hull)
    edge_mask = d <= threshold

    n_edge          = edge_mask.sum()
    n_edge_chained  = df.loc[edge_mask, "_chained"].sum()
    n_edge_unchained = n_edge - n_edge_chained

    frac_edge_chained   = n_edge_chained   / n_edge if n_edge else 0
    frac_edge_unchained = n_edge_unchained / n_edge if n_edge else 0

    # Interior cells
    n_interior           = n_total - n_edge
    n_interior_chained   = n_chained - n_edge_chained
    n_interior_unchained = n_unchained - n_edge_unchained
    frac_int_chained     = n_interior_chained   / n_interior if n_interior else 0
    frac_int_unchained   = n_interior_unchained / n_interior if n_interior else 0

    return {
        "n_total":              n_total,
        "n_chained":            int(n_chained),
        "n_unchained":          int(n_unchained),
        "frac_chained_all":     frac_chained_all,
        "frac_unchained_all":   frac_unchained_all,
        "n_edge":               int(n_edge),
        "n_edge_chained":       int(n_edge_chained),
        "n_edge_unchained":     int(n_edge_unchained),
        "frac_edge_chained":    frac_edge_chained,
        "frac_edge_unchained":  frac_edge_unchained,
        "n_interior":           int(n_interior),
        "frac_int_chained":     frac_int_chained,
        "frac_int_unchained":   frac_int_unchained,
        "hull_area":            hull.volume,   # 'volume' = area in 2-D
        "threshold_um":         threshold,
    }


def average_metrics(metric_list):
    """Average a list of metric dicts (from multiple repeats)."""
    keys = metric_list[0].keys()
    return {k: np.mean([m[k] for m in metric_list]) for k in keys}


# ── Print table ────────────────────────────────────────────────────────────────
def print_table(gf_list, avg_metrics, std_metrics=None):
    header = (f"{'GF':>5}  {'Total':>6}  {'Edge N':>6}  "
              f"{'Edge% C':>8}  {'Edge% U':>8}  {'All% C':>7}  {'All% U':>7}  "
              f"{'Edge winner':>12}")
    print("\n" + header)
    print("─" * len(header))
    for gf, m in zip(gf_list, avg_metrics):
        ec  = m["frac_edge_chained"]
        eu  = m["frac_edge_unchained"]
        winner = "CHAINED" if ec > eu else "UNCHAINED" if eu > ec else "TIE"
        print(
            f"{gf:>5.2f}  "
            f"{m['n_total']:>6.0f}  "
            f"{m['n_edge']:>6.0f}  "
            f"{ec:>8.1%}  "
            f"{eu:>8.1%}  "
            f"{m['frac_chained_all']:>7.1%}  "
            f"{m['frac_unchained_all']:>7.1%}  "
            f"{winner:>12}"
        )
    print()


# ── Plotting ───────────────────────────────────────────────────────────────────
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


def make_plot(gf_list, avg_metrics, std_metrics=None, save_path=None):
    gf  = np.array(gf_list)
    ms  = avg_metrics
    std = std_metrics  # may be None

    ec  = np.array([m["frac_edge_chained"]   for m in ms])
    eu  = np.array([m["frac_edge_unchained"]  for m in ms])
    ac  = np.array([m["frac_chained_all"]     for m in ms])
    au  = np.array([m["frac_unchained_all"]   for m in ms])
    ic  = np.array([m["frac_int_chained"]     for m in ms])
    iu  = np.array([m["frac_int_unchained"]   for m in ms])
    ne  = np.array([m["n_edge"]               for m in ms])
    nt  = np.array([m["n_total"]              for m in ms])

    lw, mks = 2.2, 6
    bar_w   = (gf[1] - gf[0]) * 0.38 if len(gf) > 1 else 0.05

    fig = plt.figure(figsize=(13, 8.5), facecolor=DARK_BG)
    fig.suptitle(
        "Colony Edge Dominance — Chained vs Unchained across Growth Factor",
        color=TEXT_COL, fontsize=12, fontweight="normal", y=0.98,
    )
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.38)

    xticks_kw = dict(rotation=35, ha="right", fontsize=7.5)

    # ── 1. Stacked area: EDGE fraction ─────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.stackplot(gf, eu, ec,
                  labels=["Unchained (edge)", "Chained (edge)"],
                  colors=[C_UNCHAINED, C_CHAINED], alpha=0.78)
    ax1.axhline(0.5, color="white", lw=0.7, ls=":", alpha=0.45)
    ax1.set_ylim(0, 1)
    ax1.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    style_ax(ax1, "Growth factor", "Fraction of edge cells",
             "Edge cell fraction (stacked)")
    ax1.legend(fontsize=7.5, facecolor=PANEL_BG, edgecolor=SPINE,
               labelcolor=TEXT_COL, loc="upper left")
    ax1.set_xticks(gf); ax1.set_xticklabels(gf, **xticks_kw)

    # ── 2. Line: edge vs interior vs whole colony ──────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    kw = dict(lw=lw, ms=mks)
    ax2.plot(gf, ec, "o-",  color=C_CHAINED,   label="Chained — edge",     **kw)
    ax2.plot(gf, ic, "o--", color=C_EDGE_C,    label="Chained — interior", **kw, alpha=0.6)
    ax2.plot(gf, eu, "s-",  color=C_UNCHAINED, label="Unchained — edge",   **kw)
    ax2.plot(gf, iu, "s--", color=C_EDGE_U,    label="Unchained — interior",**kw, alpha=0.6)
    ax2.axhline(0.5, color="white", lw=0.7, ls=":", alpha=0.45)
    ax2.set_ylim(0, 1)
    ax2.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    style_ax(ax2, "Growth factor", "Fraction", "Edge vs Interior fraction")
    ax2.legend(fontsize=6.5, facecolor=PANEL_BG, edgecolor=SPINE,
               labelcolor=TEXT_COL, ncol=2)
    ax2.set_xticks(gf); ax2.set_xticklabels(gf, **xticks_kw)

    # ── 3. Bar: edge enrichment (edge% – colony%) ─────────────────────────
    #  Positive = overrepresented on edge relative to whole colony
    ax3 = fig.add_subplot(gs[0, 2])
    enrich_c = ec - ac
    enrich_u = eu - au
    ax3.bar(gf - bar_w/2, enrich_c, bar_w, color=C_CHAINED,
            alpha=0.85, label="Chained enrichment")
    ax3.bar(gf + bar_w/2, enrich_u, bar_w, color=C_UNCHAINED,
            alpha=0.85, label="Unchained enrichment")
    ax3.axhline(0, color="white", lw=0.8, alpha=0.5)
    ax3.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, _: f"{x:+.0%}")
    )
    style_ax(ax3, "Growth factor", "Edge fraction − colony fraction",
             "Edge enrichment\n(+ = over-represented on periphery)")
    ax3.legend(fontsize=7.5, facecolor=PANEL_BG, edgecolor=SPINE, labelcolor=TEXT_COL)
    ax3.set_xticks(gf); ax3.set_xticklabels(gf, **xticks_kw)

    # ── 4. Cell counts: total / edge / interior ────────────────────────────
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(gf, nt, "^-", color=C_TOTAL,     lw=lw, ms=mks, label="Total")
    ax4.plot(gf, ne, "o-", color="#cc77ff",   lw=lw, ms=mks, label="Edge")
    ax4.plot(gf, nt - ne, "s-", color="#888", lw=lw, ms=mks, label="Interior")
    style_ax(ax4, "Growth factor", "Number of cells", "Absolute cell counts")
    ax4.legend(fontsize=7.5, facecolor=PANEL_BG, edgecolor=SPINE, labelcolor=TEXT_COL)
    ax4.set_xticks(gf); ax4.set_xticklabels(gf, **xticks_kw)

    # ── 5. Stacked bar: whole-colony composition ───────────────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.stackplot(gf, au, ac,
                  labels=["Unchained (all)", "Chained (all)"],
                  colors=[C_UNCHAINED, C_CHAINED], alpha=0.78)
    ax5.axhline(0.5, color="white", lw=0.7, ls=":", alpha=0.45)
    ax5.set_ylim(0, 1)
    ax5.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=0))
    style_ax(ax5, "Growth factor", "Fraction of all cells",
             "Whole-colony composition")
    ax5.legend(fontsize=7.5, facecolor=PANEL_BG, edgecolor=SPINE,
               labelcolor=TEXT_COL, loc="upper left")
    ax5.set_xticks(gf); ax5.set_xticklabels(gf, **xticks_kw)

    # ── 6. Edge cell count by type ─────────────────────────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    n_ec = np.array([m["n_edge_chained"]   for m in ms])
    n_eu = np.array([m["n_edge_unchained"] for m in ms])
    ax6.bar(gf - bar_w/2, n_ec, bar_w, color=C_CHAINED,   alpha=0.85, label="Chained")
    ax6.bar(gf + bar_w/2, n_eu, bar_w, color=C_UNCHAINED, alpha=0.85, label="Unchained")
    style_ax(ax6, "Growth factor", "Number of edge cells",
             "Edge cell counts by type")
    ax6.legend(fontsize=7.5, facecolor=PANEL_BG, edgecolor=SPINE, labelcolor=TEXT_COL)
    ax6.set_xticks(gf); ax6.set_xticklabels(gf, **xticks_kw)

    if save_path:
        fig.savefig(save_path, dpi=150, facecolor=fig.get_facecolor(),
                    bbox_inches="tight")
        print(f"Saved → {save_path}")
    else:
        plt.show()
    plt.close(fig)


# ── Entry point ────────────────────────────────────────────────────────────────
def find_last_biofilm(repeat_dir):
    """
    Return the last biofilm_*.dat file (by name sort) inside repeat_dir,
    or None if none exist.
    """
    files = sorted(Path(repeat_dir).glob("biofilm_*.dat"))
    return files[-1] if files else None


def discover_files(base_dir, condition="Ch50_unch50"):
    """
    Walk the directory tree:
        <base_dir>/growthRateMultiplier*/<condition>/repeat*/

    For each repeat dir, pick the last biofilm_*.dat file.
    Returns dict: gf (float) → list of Path objects (one per repeat).
    """
    base = Path(base_dir)
    groups = defaultdict(list)

    gf_dirs = sorted(base.glob("growthRateMultiplier*"))
    if not gf_dirs:
        raise FileNotFoundError(
            f"No 'growthRateMultiplier*' directories found under:\n  {base}"
        )

    for gf_dir in gf_dirs:
        try:
            gf = parse_growth_factor(str(gf_dir))
        except ValueError:
            print(f"  Skipping (no parseable GF): {gf_dir.name}")
            continue

        cond_dir = gf_dir / condition
        if not cond_dir.is_dir():
            print(f"  Warning: condition dir not found: {cond_dir}")
            continue

        repeat_dirs = sorted(cond_dir.glob("repeat*"))
        if not repeat_dirs:
            print(f"  Warning: no repeat* dirs in {cond_dir}")
            continue

        for rep_dir in repeat_dirs:
            last = find_last_biofilm(rep_dir)
            if last:
                groups[gf].append(last)
                print(f"  gf={gf:.2f}  {rep_dir.name}  →  {last.name}")
            else:
                print(f"  Warning: no biofilm_*.dat in {rep_dir}")

    return groups


def main():
    parser = argparse.ArgumentParser(
        description="Quantify edge dominance of chained vs unchained bacteria.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  # Point at the base data directory — script finds everything automatically:
  python edge_dominance.py --data ../output/data_production/

  # Different condition subdirectory name:
  python edge_dominance.py --data ../output/data_production/ --condition Ch50_unch50

  # Save plot:
  python edge_dominance.py --data ../output/data_production/ --save

  # Looser/tighter edge band (default 2 × mean radius):
  python edge_dominance.py --data ../output/data_production/ --edge-threshold 3
        """,
    )
    parser.add_argument(
        "--data", required=True, metavar="DIR",
        help="Base data directory containing growthRateMultiplier* subdirs",
    )
    parser.add_argument(
        "--condition", default="Ch50_unch50", metavar="NAME",
        help="Condition subdirectory name (default: Ch50_unch50)",
    )
    parser.add_argument(
        "--edge-threshold", type=float, default=2.0, metavar="N",
        help="Edge = cells within N × mean_radius of hull boundary (default: 2.0)",
    )
    parser.add_argument(
        "--save", action="store_true",
        help="Save plot as PNG instead of opening a window",
    )
    args = parser.parse_args()

    print(f"\nScanning: {Path(args.data).resolve()}")
    print(f"Condition: {args.condition}")
    print(f"Edge threshold: {args.edge_threshold} × mean_radius\n")

    try:
        groups = discover_files(args.data, args.condition)
    except FileNotFoundError as e:
        print(f"Error: {e}"); sys.exit(1)

    if not groups:
        print("No usable data found."); sys.exit(1)

    gf_sorted = sorted(groups.keys())
    print(f"\nGrowth factors found: {gf_sorted}\n")

    avg_metrics_list = []
    std_metrics_list = []
    valid_gf         = []

    for gf in gf_sorted:
        file_list = groups[gf]
        repeat_metrics = []
        for fp in file_list:
            try:
                df = load_dat(fp)
                m  = compute_edge_metrics(df, args.edge_threshold)
                if m:
                    repeat_metrics.append(m)
            except Exception as e:
                print(f"  Error loading {fp.name}: {e}")

        if not repeat_metrics:
            print(f"  gf={gf:.2f}: no valid data, skipping.")
            continue

        n_rep = len(repeat_metrics)
        avg   = average_metrics(repeat_metrics)
        std   = (
            {k: np.std([m[k] for m in repeat_metrics]) for k in avg}
            if n_rep > 1
            else {k: 0.0 for k in avg}
        )

        avg_metrics_list.append(avg)
        std_metrics_list.append(std)
        valid_gf.append(gf)

        print(f"  gf={gf:.2f}  |  {n_rep} repeat(s)  |  "
              f"avg total={avg['n_total']:.0f}  |  "
              f"edge chained={avg['frac_edge_chained']:.1%}  "
              f"edge unchained={avg['frac_edge_unchained']:.1%}")

    if not avg_metrics_list:
        print("No valid data to plot."); sys.exit(1)

    print_table(valid_gf, avg_metrics_list)

    save_path = (Path(args.data) / "edge_dominance.png") if args.save else None
    make_plot(valid_gf, avg_metrics_list, std_metrics_list, save_path=save_path)


if __name__ == "__main__":
    main()
