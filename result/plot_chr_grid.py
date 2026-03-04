#!/usr/bin/env python3
import os

# --- SLURM / HPC stability: avoid extra threads from BLAS ---
#os.environ["OMP_NUM_THREADS"] = "1"
#os.environ["OPENBLAS_NUM_THREADS"] = "1"
#os.environ["MKL_NUM_THREADS"] = "1"
#os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")  # non-interactive backend: save only
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COLORS = {
    "TR1": "#ff3b9d",
    "knob180": "#33cc33",
    "CentC": "#f4c430",
    "Contig": "#bdbdbd",
}

def read_bed(path: str) -> pd.DataFrame:
    """Read BED-like file safely. Uses first 3 columns: chr, start, end."""
    cols = ["chr", "start", "end"]
    if (not os.path.exists(path)) or (os.path.getsize(path) == 0):
        return pd.DataFrame(columns=cols)

    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] < 3:
        return pd.DataFrame(columns=cols)

    df = df.iloc[:, :3].copy()
    df.columns = cols
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start", "end"])
    return df

def read_agp(path: str) -> pd.DataFrame:
    """Read AGP file safely. Keeps chr, start, end, component ID."""
    cols = ["chr", "start", "end", "ID"]
    if (not os.path.exists(path)) or (os.path.getsize(path) == 0):
        return pd.DataFrame(columns=cols)

    agp = pd.read_csv(path, sep="\t", header=None, comment="#")
    if agp.shape[1] < 6:
        return pd.DataFrame(columns=cols)

    agp = agp.iloc[:, [0, 1, 2, 5]].copy()
    agp.columns = cols
    agp["start"] = pd.to_numeric(agp["start"], errors="coerce")
    agp["end"]   = pd.to_numeric(agp["end"], errors="coerce")
    agp = agp.dropna(subset=["start", "end"])
    return agp

def get_chr_length_mb_perchr(base_dir: str, chr_name: str) -> float:
    """
    Reads length from base_dir/{chr}.txt
    Accepts:
      chr7 <tab> 218916177
    or:
      218916177
    """
    length_path = os.path.join(base_dir, f"{chr_name}.txt")
    if (not os.path.exists(length_path)) or (os.path.getsize(length_path) == 0):
        raise FileNotFoundError(f"Missing length file: {length_path}")

    df = pd.read_csv(length_path, sep=r"\s+", header=None)  # tabs OR spaces

    if df.shape[1] >= 2:
        bp = int(df.iloc[0, 1])
    else:
        bp = int(df.iloc[0, 0])

    return bp / 1e6

def add_track(ax, df: pd.DataFrame, color: str, repeat_y: float, repeat_h: float,
              z: int = 10, alpha: float = 0.95):
    """Add repeat rectangles to axis."""
    if df.empty:
        return
    for _, r in df.iterrows():
        x0 = float(r["start"]) / 1e6
        w  = (float(r["end"]) - float(r["start"])) / 1e6
        if w <= 0:
            continue
        ax.add_patch(
            Rectangle(
                (x0, repeat_y),
                w,
                repeat_h,
                facecolor=color,
                edgecolor="none",
                alpha=alpha,
                zorder=z,
            )
        )

def plot_one_sample(ax, base_dir: str, sample: str, chr_name: str, tick_step: int = 10):
    """
    base_dir is the repeat_blast directory inside the sample folder.
    Example: ragtag_B73D1_v3/repeat_blast
    """
    # repeat beds
    tr1_path   = os.path.join(base_dir, f"intervals_{sample}_TR1.bed")
    knob_path  = os.path.join(base_dir, f"intervals_{sample}_knob180.bed")
    centc_path = os.path.join(base_dir, f"intervals_{sample}_CentC.bed")

    # agp
    agp_path = os.path.join(base_dir, f"{chr_name}.agp")

    tr1   = read_bed(tr1_path)
    knob  = read_bed(knob_path)
    centc = read_bed(centc_path)
    agp   = read_agp(agp_path)

    total_length_mb = get_chr_length_mb_perchr(base_dir, chr_name)

    # geometry
    bar_height = 0.35
    y_center = 0.0
    y_bottom = y_center - bar_height / 2

    contig_h = bar_height * 0.35
    repeat_h = bar_height * 0.95
    gap = 0.05

    contig_y = y_bottom - contig_h - gap
    repeat_y = y_center - repeat_h / 2

    # chromosome backbone
    ax.add_patch(
        Rectangle(
            (0, y_bottom),
            total_length_mb,
            bar_height,
            facecolor="#1673A3",
            edgecolor="none",
            zorder=1,
        )
    )

    # AGP contigs
    if not agp.empty:
        for _, row in agp.iterrows():
            x0 = float(row["start"]) / 1e6
            w  = (float(row["end"]) - float(row["start"])) / 1e6
            if w <= 0:
                continue
            ax.add_patch(
                Rectangle(
                    (x0, contig_y),
                    w,
                    contig_h,
                    facecolor=COLORS["Contig"],
                    edgecolor="black",
                    linewidth=0.35,
                    zorder=2,
                )
            )

    # repeats on backbone
    add_track(ax, tr1,   COLORS["TR1"],     repeat_y, repeat_h, z=10)
    add_track(ax, knob,  COLORS["knob180"], repeat_y, repeat_h, z=10)
    add_track(ax, centc, COLORS["CentC"],   repeat_y, repeat_h, z=10)

    # formatting
    ax.set_xlim(0, total_length_mb * 1.02)
    ax.set_ylim(contig_y - 0.12, y_bottom + bar_height + 0.12)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, total_length_mb + tick_step, tick_step))
    ax.set_title(sample, fontsize=10)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="linename.lst (must contain exactly 8 sample names)")
    ap.add_argument("--chr", required=True, help="chr1/chr2/chr7/chr8")
    ap.add_argument("--root_pattern", default="ragtag_{sample}_v3",
                    help='sample directory pattern (default "ragtag_{sample}_v3")')
    ap.add_argument("--subdir", default="repeat_blast",
                    help="subdir under sample dir where files live (default repeat_blast)")
    ap.add_argument("--outdir", default="chr_plots", help="output directory for figures")
    ap.add_argument("--tick_step", type=int, default=10)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    # load sample names
    with open(args.samples) as f:
        samples = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]

    if len(samples) != 8:
        raise ValueError(f"Expected 8 samples in {args.samples}, found {len(samples)}")

    os.makedirs(args.outdir, exist_ok=True)

    fig, axes = plt.subplots(4, 2, figsize=(14, 10))
    axes = axes.flatten()

    for i, sample in enumerate(samples):
        ax = axes[i]
        sample_root = args.root_pattern.format(sample=sample)
        base_dir = os.path.join(sample_root, args.subdir)

        plot_one_sample(ax, base_dir, sample, args.chr, tick_step=args.tick_step)

        # only bottom row shows x labels
        if i < 6:
            ax.set_xticklabels([])
            ax.set_xlabel("")
        else:
            ax.set_xlabel("Position (Mb)")

    fig.suptitle(f"{args.chr} repeats + AGP contigs (8 individuals; 2×4 grid)", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out_pdf = os.path.join(args.outdir, f"{args.chr}.8samples.grid.pdf")
    out_png = os.path.join(args.outdir, f"{args.chr}.8samples.grid.png")

    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")

    print("Saved:", out_pdf)
    print("Saved:", out_png)

    plt.close(fig)  # IMPORTANT: free memory on SLURM

if __name__ == "__main__":
    main()
