"""
biofilm_movie.py
----------------
Render all biofilm_*.dat files in a repeat directory into a movie (MP4 or GIF).

Usage
-----
    # MP4 (requires ffmpeg):
    python biofilm_movie.py --repeat-dir ../output/data_production/growthRateMultiplier1_0p7/Ch50_unch50/repeat0/

    # GIF (no ffmpeg needed):
    python biofilm_movie.py --repeat-dir .../repeat0/ --format gif

    # Faster / slower playback:
    python biofilm_movie.py --repeat-dir .../repeat0/ --fps 10

    # Only every Nth file (speeds up rendering for large datasets):
    python biofilm_movie.py --repeat-dir .../repeat0/ --step 5

    # Custom output path:
    python biofilm_movie.py --repeat-dir .../repeat0/ --output my_colony.mp4

Requirements
------------
    pip install matplotlib numpy pandas scipy tqdm pillow
    # For MP4: ffmpeg must be installed and on PATH
    #   macOS:   brew install ffmpeg
    #   Ubuntu:  sudo apt install ffmpeg
    #   Windows: https://ffmpeg.org/download.html
"""

import re
import sys
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for rendering
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as transforms
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

warnings.filterwarnings("ignore")

# ── Colour palette (Wong – colorblind safe) ────────────────────────────────────
COLOUR = {"unchained": "#E69F00", "chained": "#0072B2"}
EDGE   = {"unchained": "#F5C842", "chained": "#56B4E9"}
DARK_BG   = "#050d0a"
PANEL_BG  = "#061a10"
TEXT_COL  = "#c0ffe8"
GRID_COL  = "#0af5a015"


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


# ── Drawing ────────────────────────────────────────────────────────────────────
def draw_capsule(ax, cx, cy, length, radius, angle_rad, facecolor, edgecolor):
    rect = mpatches.FancyBboxPatch(
        (-length / 2, -radius),
        length, 2 * radius,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=0.6,
        edgecolor=edgecolor,
        facecolor=facecolor,
        alpha=0.92,
        zorder=2,
    )
    t = (
        transforms.Affine2D()
        .rotate(angle_rad)
        .translate(cx, cy)
        + ax.transData
    )
    rect.set_transform(t)
    ax.add_patch(rect)
    return rect


def draw_frame(ax, df, xlim, ylim):
    """Clear axis and redraw all cells."""
    ax.cla()
    ax.set_facecolor(PANEL_BG)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.grid(color=GRID_COL, linewidth=0.5, zorder=0)
    for sp in ax.spines.values():
        sp.set_edgecolor("#ffffff20")
    ax.tick_params(colors="#0af5a060", labelsize=7)
    ax.xaxis.label.set_color("#0af5a080")
    ax.yaxis.label.set_color("#0af5a080")
    ax.set_xlabel("x (µm)", fontsize=8)
    ax.set_ylabel("y (µm)", fontsize=8)

    # Chain links
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
                color="#5599ff", lw=0.6, ls="--", alpha=0.4, zorder=1,
            )

    # Cells
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


# ── Compute fixed axes bounds across all frames ────────────────────────────────
def global_bounds(files, pad=1.5):
    xmins, xmaxs, ymins, ymaxs = [], [], [], []
    print("  Computing global bounds …")
    for fp in files:
        try:
            df = load_dat(fp)
            xmins.append(df["pos_x"].min() - df["length"].max() / 2)
            xmaxs.append(df["pos_x"].max() + df["length"].max() / 2)
            ymins.append(df["pos_y"].min() - df["length"].max() / 2)
            ymaxs.append(df["pos_y"].max() + df["length"].max() / 2)
        except Exception:
            pass
    return (min(xmins) - pad, max(xmaxs) + pad,
            min(ymins) - pad, max(ymaxs) + pad)


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Render biofilm_*.dat frames into a movie."
    )
    parser.add_argument("--repeat-dir", required=True, metavar="DIR",
                        help="Directory containing biofilm_*.dat files")
    parser.add_argument("--output", default=None, metavar="FILE",
                        help="Output file path (default: <repeat-dir>/biofilm_movie.mp4)")
    parser.add_argument("--format", choices=["mp4", "gif"], default="mp4",
                        help="Output format (default: mp4)")
    parser.add_argument("--fps", type=int, default=15,
                        help="Frames per second (default: 15)")
    parser.add_argument("--step", type=int, default=1,
                        help="Use every Nth file (default: 1 = all files)")
    parser.add_argument("--dt", type=float, default=0.1,
                        help="Time per file step in hours (default: 0.1)")
    parser.add_argument("--dpi", type=int, default=120,
                        help="Resolution in DPI (default: 120)")
    args = parser.parse_args()

    repeat_dir = Path(args.repeat_dir)
    if not repeat_dir.is_dir():
        print(f"Error: directory not found: {repeat_dir}")
        sys.exit(1)

    # Collect and sort files
    files = sorted(repeat_dir.glob("biofilm_*.dat"),
                   key=lambda f: parse_timestep(f))
    files = files[::args.step]

    if not files:
        print(f"No biofilm_*.dat files found in {repeat_dir}")
        sys.exit(1)

    print(f"\nFound {len(files)} files  (step={args.step})")
    print(f"  First: {files[0].name}   Last: {files[-1].name}")
    print(f"  Duration: {parse_timestep(files[-1]) * args.dt:.1f} hrs at {args.fps} fps\n")

    # Output path
    fmt = args.format
    if args.output:
        out_path = Path(args.output)
    else:
        out_path = repeat_dir / f"biofilm_movie.{fmt}"

    # Global bounds (so axes don't jump between frames)
    x0, x1, y0, y1 = global_bounds(files)
    xlim, ylim = (x0, x1), (y0, y1)

    # Set up figure
    fig, ax = plt.subplots(figsize=(7, 6), facecolor=DARK_BG)
    fig.subplots_adjust(top=0.91, bottom=0.09, left=0.10, right=0.97)

    # Legend
    legend_handles = [
        mpatches.Patch(facecolor=COLOUR["chained"],   edgecolor=EDGE["chained"],
                       label="Chained"),
        mpatches.Patch(facecolor=COLOUR["unchained"], edgecolor=EDGE["unchained"],
                       label="Unchained"),
        plt.Line2D([0], [0], color="#5599ff", lw=0.8, ls="--", label="Chain link"),
    ]
    leg = ax.legend(
        handles=legend_handles, loc="upper right",
        fontsize=7, facecolor="#071a10", edgecolor="#0af5a040",
        labelcolor=TEXT_COL,
    )

    title_obj = fig.suptitle("", color=TEXT_COL, fontsize=10, fontweight="normal")

    # Cache DataFrames to avoid re-reading during animation
    print("  Loading all frames …")
    frames_data = []
    for i, fp in enumerate(files):
        try:
            df = load_dat(fp)
            ts = parse_timestep(fp)
            frames_data.append((ts, df))
            if (i + 1) % 20 == 0 or i == len(files) - 1:
                print(f"    {i + 1}/{len(files)} loaded", end="\r")
        except Exception as e:
            print(f"\n  Warning: skipping {fp.name}: {e}")
    print(f"\n  {len(frames_data)} frames ready.\n")

    def update(frame_idx):
        ts, df = frames_data[frame_idx]
        time_hrs = ts * args.dt
        n_c  = (df["_type"] == "chained").sum()
        n_u  = (df["_type"] == "unchained").sum()
        draw_frame(ax, df, xlim, ylim)
        # Re-add legend (draw_frame clears the axis)
        ax.legend(handles=legend_handles, loc="upper right",
                  fontsize=7, facecolor="#071a10",
                  edgecolor="#0af5a040", labelcolor=TEXT_COL)
        title_obj.set_text(
            f"t = {time_hrs:.1f} hrs   |   "
            f"cells: {len(df)}   chained: {n_c}   unchained: {n_u}"
        )
        return ax.patches

    print(f"  Rendering {len(frames_data)} frames …")
    anim = FuncAnimation(
        fig, update,
        frames=len(frames_data),
        interval=1000 / args.fps,
        blit=False,
    )

    if fmt == "mp4":
        writer = FFMpegWriter(
            fps=args.fps,
            metadata={"title": "Biofilm colony"},
            bitrate=1800,
            extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"],
        )
        print(f"  Writing MP4 → {out_path}")
        anim.save(str(out_path), writer=writer, dpi=args.dpi,
                  savefig_kwargs={"facecolor": DARK_BG})
    else:
        writer = PillowWriter(fps=args.fps)
        print(f"  Writing GIF → {out_path}")
        anim.save(str(out_path), writer=writer, dpi=args.dpi,
                  savefig_kwargs={"facecolor": DARK_BG})

    plt.close(fig)
    print(f"\nDone!  →  {out_path}\n")


if __name__ == "__main__":
    main()
