"""
snapshots.py
------------
Render every biofilm_*.dat file in a directory as a high-quality PNG snapshot.
Uses the same capsule drawing style as plot_biofilm.py.

Once you've seen them all you can pick the ones that best show buckling
and the final loopy colony.

Usage
-----
    python snapshots.py --repeat-dir path/to/repeat0/

    # Custom output folder:
    python snapshots.py --repeat-dir path/to/repeat0/ --output-dir ./snapshots/

    # Only every Nth file (if you have hundreds of frames):
    python snapshots.py --repeat-dir path/to/repeat0/ --step 5

    # Higher resolution (default 200 dpi):
    python snapshots.py --repeat-dir path/to/repeat0/ --dpi 300

    # Time per step in hrs (default 0.1):
    python snapshots.py --repeat-dir path/to/repeat0/ --dt 0.1

Requirements
------------
    pip install matplotlib numpy pandas tqdm
"""

import re
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms

# ── Colour palette (Wong – colorblind safe) ────────────────────────────────────
COLOUR = {"unchained": "#E69F00", "chained": "#0072B2"}
EDGE   = {"unchained": "#F5C842", "chained": "#56B4E9"}
DARK_BG  = "#050d0a"
PANEL_BG = "#061a10"
TEXT_COL = "#c0ffe8"
GRID_COL = "#0af5a015"


# ── Data loading ───────────────────────────────────────────────────────────────
def load_dat(path):
    df = pd.read_csv(
        path, sep=r"\s+",
        na_values=["None", "-", ""],
        dtype={"lower_link": "Int64", "upper_link": "Int64"},
    )
    required = {"cell_id", "length", "radius", "pos_x", "pos_y",
                "ori_x", "ori_y", "lower_link", "upper_link"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    df["_type"] = df.apply(
        lambda r: "chained"
        if pd.notna(r["lower_link"]) or pd.notna(r["upper_link"])
        else "unchained",
        axis=1,
    )
    return df


def parse_timestep(filename):
    m = re.search(r"biofilm_(\d+)\.dat", Path(filename).name)
    return int(m.group(1)) if m else 0


# ── Capsule drawing ────────────────────────────────────────────────────────────
def draw_capsule(ax, cx, cy, length, radius, angle_rad, facecolor, edgecolor):
    rect = mpatches.FancyBboxPatch(
        (-length / 2, -radius),
        length, 2 * radius,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=0.5,
        edgecolor=edgecolor,
        facecolor=facecolor,
        alpha=0.93,
        zorder=2,
    )
    t = (
        mtransforms.Affine2D()
        .rotate(angle_rad)
        .translate(cx, cy)
        + ax.transData
    )
    rect.set_transform(t)
    ax.add_patch(rect)


def draw_chain_links(ax, df):
    id_map = df.set_index("cell_id")
    seen   = set()
    for _, row in df.iterrows():
        for link in [row["lower_link"], row["upper_link"]]:
            if pd.isna(link):
                continue
            link = int(link)
            key  = tuple(sorted([int(row["cell_id"]), link]))
            if key in seen or link not in id_map.index:
                continue
            seen.add(key)
            other = id_map.loc[link]
            ax.plot(
                [row["pos_x"], other["pos_x"]],
                [row["pos_y"], other["pos_y"]],
                color="#5599ff", lw=0.5, ls="--", alpha=0.35, zorder=1,
            )


# ── Compute global bounds across all frames ────────────────────────────────────
def global_bounds(files, pad=2.0):
    xs_min, xs_max, ys_min, ys_max = [], [], [], []
    print("  Pre-scanning bounds across all frames …")
    for fp in files:
        try:
            df = load_dat(fp)
            half = df["length"].max() / 2 + df["radius"].max()
            xs_min.append(df["pos_x"].min() - half)
            xs_max.append(df["pos_x"].max() + half)
            ys_min.append(df["pos_y"].min() - half)
            ys_max.append(df["pos_y"].max() + half)
        except Exception:
            pass
    return (
        min(xs_min) - pad, max(xs_max) + pad,
        min(ys_min) - pad, max(ys_max) + pad,
    )


# ── Render one snapshot ────────────────────────────────────────────────────────
def render_snapshot(df, ts, dt, xlim, ylim, out_path, dpi):
    time_hrs   = ts * dt
    n_total    = len(df)
    n_chained  = (df["_type"] == "chained").sum()
    n_unchained = n_total - n_chained

    # Figure size proportional to data aspect ratio
    data_w = xlim[1] - xlim[0]
    data_h = ylim[1] - ylim[0]
    aspect = data_h / data_w
    fig_w  = 8.0
    fig_h  = max(5.0, fig_w * aspect)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=DARK_BG)
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.09, right=0.97)

    ax.set_facecolor(PANEL_BG)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.grid(color=GRID_COL, linewidth=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_edgecolor("#ffffff20")
    ax.tick_params(colors="#0af5a070", labelsize=8)
    ax.set_xlabel("x (µm)", color="#0af5a090", fontsize=9)
    ax.set_ylabel("y (µm)", color="#0af5a090", fontsize=9)

    # Draw chain links then capsules
    draw_chain_links(ax, df)
    for _, row in df.iterrows():
        ct = row["_type"]
        draw_capsule(
            ax,
            cx=row["pos_x"], cy=row["pos_y"],
            length=row["length"], radius=row["radius"],
            angle_rad=np.arctan2(row["ori_y"], row["ori_x"]),
            facecolor=COLOUR[ct],
            edgecolor=EDGE[ct],
        )

    # Legend
    legend_handles = [
        mpatches.Patch(facecolor=COLOUR["chained"],   edgecolor=EDGE["chained"],
                       label=f"Chained  ({n_chained})"),
        mpatches.Patch(facecolor=COLOUR["unchained"], edgecolor=EDGE["unchained"],
                       label=f"Unchained  ({n_unchained})"),
        plt.Line2D([0], [0], color="#5599ff", lw=0.8, ls="--", label="Chain link"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right", fontsize=8,
        facecolor="#071a10", edgecolor="#0af5a040",
        labelcolor=TEXT_COL,
    )

    # Title
    fig.suptitle(
        f"t = {time_hrs:.1f} hrs   ·   {n_total} cells",
        color=TEXT_COL, fontsize=11, fontweight="normal",
    )

    fig.savefig(out_path, dpi=dpi, facecolor=DARK_BG, bbox_inches="tight")
    plt.close(fig)


# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Render every biofilm_*.dat frame as a high-quality PNG."
    )
    parser.add_argument("--repeat-dir", required=True, metavar="DIR",
                        help="Directory containing biofilm_*.dat files")
    parser.add_argument("--output-dir", default=None, metavar="DIR",
                        help="Where to save PNGs (default: <repeat-dir>/snapshots/)")
    parser.add_argument("--step", type=int, default=1,
                        help="Use every Nth file (default: 1 = all)")
    parser.add_argument("--dt", type=float, default=0.1,
                        help="Hours per timestep (default: 0.1)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="Image resolution in DPI (default: 200)")
    parser.add_argument("--no-fixed-bounds", action="store_true",
                        help="Let axes rescale each frame (default: fixed bounds)")
    args = parser.parse_args()

    repeat_dir = Path(args.repeat_dir)
    if not repeat_dir.is_dir():
        print(f"Error: not a directory: {repeat_dir}"); sys.exit(1)

    # Output folder
    out_dir = Path(args.output_dir) if args.output_dir else repeat_dir / "snapshots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Collect files
    files = sorted(repeat_dir.glob("biofilm_*.dat"),
                   key=lambda f: parse_timestep(f))[::args.step]
    if not files:
        print(f"No biofilm_*.dat files found in {repeat_dir}"); sys.exit(1)

    print(f"\nRepeat dir : {repeat_dir}")
    print(f"Output dir : {out_dir}")
    print(f"Files found: {len(files)}  (step={args.step})")
    print(f"DPI        : {args.dpi}")
    print(f"dt         : {args.dt} hrs\n")

    # Global bounds (so all frames are on the same axes)
    if not args.no_fixed_bounds:
        x0, x1, y0, y1 = global_bounds(files)
        xlim, ylim = (x0, x1), (y0, y1)
        print(f"  Bounds: x=[{x0:.1f}, {x1:.1f}]  y=[{y0:.1f}, {y1:.1f}]\n")
    else:
        xlim = ylim = None

    # Render
    print(f"  Rendering {len(files)} snapshots …\n")
    for i, fp in enumerate(files):
        ts = parse_timestep(fp)
        out_path = out_dir / f"snapshot_{ts:06d}.png"

        try:
            df = load_dat(fp)

            # Per-frame bounds if not fixed
            if xlim is None:
                pad  = 1.5
                half = df["length"].max() / 2 + df["radius"].max()
                xlim_ = (df["pos_x"].min() - half - pad,
                          df["pos_x"].max() + half + pad)
                ylim_ = (df["pos_y"].min() - half - pad,
                          df["pos_y"].max() + half + pad)
            else:
                xlim_, ylim_ = xlim, ylim

            render_snapshot(df, ts, args.dt, xlim_, ylim_, out_path, args.dpi)

            print(f"  [{i+1:>4}/{len(files)}]  t={ts*args.dt:>6.1f} hrs  "
                  f"cells={len(df):>5}  →  {out_path.name}")

        except Exception as e:
            print(f"  [{i+1:>4}/{len(files)}]  SKIP {fp.name}: {e}")

    print(f"\nDone!  {len(files)} PNGs saved to:\n  {out_dir}\n")


if __name__ == "__main__":
    main()
