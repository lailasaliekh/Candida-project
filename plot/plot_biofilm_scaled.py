"""
plot_biofilm.py
---------------
Plot bacterial colonies from biofilm*.dat files.

Usage:
    python plot_biofilm.py biofilm001.dat
    python plot_biofilm.py biofilm*.dat          # multiple files → one plot each
    python plot_biofilm.py biofilm001.dat --save # save as PNG instead of showing

Requirements:
    pip install matplotlib numpy pandas
"""

import sys
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from matplotlib.collections import PatchCollection


# ── Colours ────────────────────────────────────────────────────────────────────
COLOUR = {
    "unchained": "#E69F00",   # orange  (Wong palette)
    "chained":   "#0072B2",   # blue    (Wong palette)
}
EDGE = {
    "unchained": "#F5C842",
    "chained":   "#56B4E9",
}


def chain_type(row):
    links = sum(1 for v in [row["lower_link"], row["upper_link"]] if pd.notna(v))
    return "chained" if links > 0 else "unchained"


# ── Parser ─────────────────────────────────────────────────────────────────────
def load_dat(path):
    df = pd.read_csv(
        path,
        sep=r"\s+",
        na_values=["None", "-", ""],
        dtype={
            "lower_link": "Int64",
            "upper_link": "Int64",
        },
    )
    required = {"cell_id", "length", "radius", "pos_x", "pos_y",
                "ori_x", "ori_y", "lower_link", "upper_link"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    df["_type"] = df.apply(chain_type, axis=1)
    return df


# ── Capsule drawing ────────────────────────────────────────────────────────────
def draw_capsule(ax, cx, cy, length, radius, angle_rad, facecolor, edgecolor, alpha=0.9, zorder=2):
    """
    Draw a capsule (stadium shape) centred at (cx, cy), oriented by angle_rad.
    Built from a rectangle + two semicircles.
    """
    import matplotlib.transforms as transforms

    dx = np.cos(angle_rad)
    dy = np.sin(angle_rad)
    half = length / 2

    # Rectangle body
    rect = mpatches.FancyBboxPatch(
        (-half, -radius),
        length, 2 * radius,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=0.8,
        edgecolor=edgecolor,
        facecolor=facecolor,
        alpha=alpha,
        zorder=zorder,
    )
    t = (
        transforms.Affine2D()
        .rotate(angle_rad)
        .translate(cx, cy)
        + ax.transData
    )
    rect.set_transform(t)
    ax.add_patch(rect)


# ── Chain links ────────────────────────────────────────────────────────────────
def draw_chains(ax, df):
    id_to_row = df.set_index("cell_id")
    seen = set()
    for _, row in df.iterrows():
        for link in [row["lower_link"], row["upper_link"]]:
            if pd.isna(link):
                continue
            link = int(link)
            key = tuple(sorted([int(row["cell_id"]), link]))
            if key in seen or link not in id_to_row.index:
                continue
            seen.add(key)
            other = id_to_row.loc[link]
            ax.plot(
                [row["pos_x"], other["pos_x"]],
                [row["pos_y"], other["pos_y"]],
                color="#5599ff",
                linewidth=0.8,
                linestyle="--",
                alpha=0.5,
                zorder=1,
            )


# ── Bounding box helper ─────────────────────────────────────────────────────────
def cell_bbox_half_extent(row):
    """
    Axis-aligned half-extent (dx, dy) of a rotated capsule, so we can pad the
    plot enough that no cell gets clipped regardless of its orientation.
    """
    angle = np.arctan2(row["ori_y"], row["ori_x"])
    half = row["length"] / 2
    r = row["radius"]
    # Half-extent of the rotated rectangle body, plus the radius in every
    # direction (covers the semicircular caps too).
    dx = half * abs(np.cos(angle)) + r
    dy = half * abs(np.sin(angle)) + r
    return dx, dy


# ── Main plot ──────────────────────────────────────────────────────────────────
def plot_colony(df, title="Bacterial Colony", save_path=None):
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Chain links first (behind cells)
    draw_chains(ax, df)

    # Cells
    for _, row in df.iterrows():
        angle = np.arctan2(row["ori_y"], row["ori_x"])
        ct = row["_type"]
        draw_capsule(
            ax,
            cx=row["pos_x"],
            cy=row["pos_y"],
            length=row["length"],
            radius=row["radius"],
            angle_rad=angle,
            facecolor=COLOUR[ct],
            edgecolor=EDGE[ct],
        )

    # No axes — just the cells
    ax.set_axis_off()
    ax.set_aspect("equal")

    # ── Auto margin (fixed) ──────────────────────────────────────────────
    # Instead of a fixed pad, compute the true rotated bounding box of every
    # capsule (length/2 + radius, adjusted for orientation) and use the
    # largest per-cell extent found anywhere in the colony. This guarantees
    # cells near the edge of the frame — especially long or diagonally
    # oriented ones during buckling — are never clipped, while still
    # cropping tightly for small/round colonies.
    half_extents = df.apply(cell_bbox_half_extent, axis=1, result_type="expand")
    half_extents.columns = ["dx", "dy"]

    x_min = (df["pos_x"] - half_extents["dx"]).min()
    x_max = (df["pos_x"] + half_extents["dx"]).max()
    y_min = (df["pos_y"] - half_extents["dy"]).min()
    y_max = (df["pos_y"] + half_extents["dy"]).max()

    # Small extra margin so edges aren't flush against the frame
    margin = 0.05 * max(x_max - x_min, y_max - y_min, 1.0)
    ax.set_xlim(x_min - margin, x_max + margin)
    ax.set_ylim(y_min - margin, y_max + margin)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, facecolor=fig.get_facecolor())
        print(f"Saved → {save_path}")
    else:
        plt.show()

    plt.close(fig)


# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Plot biofilm*.dat bacterial colonies.")
    parser.add_argument("files", nargs="+", help="One or more .dat files (globs accepted)")
    parser.add_argument("--save", action="store_true",
                        help="Save plots as PNG files instead of opening a window")
    args = parser.parse_args()

    # Expand globs (needed on Windows where the shell doesn't)
    paths = []
    for pattern in args.files:
        expanded = glob.glob(pattern)
        paths.extend(expanded if expanded else [pattern])

    if not paths:
        print("No files found.")
        sys.exit(1)

    for path in paths:
        p = Path(path)
        print(f"Loading {p.name} …")
        try:
            df = load_dat(p)
        except Exception as e:
            print(f"  Error: {e}")
            continue

        print(f"  {len(df)} cells  |  "
              f"{(df['_type']=='chained').sum()} chained  |  "
              f"{(df['_type']=='unchained').sum()} unchained")

        save_path = p.with_suffix(".pdf") if args.save else None
        plot_colony(df, title=p.name, save_path=save_path)


if __name__ == "__main__":
    main()