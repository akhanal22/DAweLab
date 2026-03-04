#!/usr/bin/env python3
import os

# non-interactive for SLURM (save only)
import matplotlib
matplotlib.use("Agg")

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COLORS = {
    "TR1": "#ff3b9d",
    "knob180": "#33cc33",
    "CentC": "#f4c430",
    "Contig": "#bdbdbd",
}

def read_bed(path):
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
    return df.dropna(subset=["start", "end"])

def read_agp(path):
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
    return agp.dropna(subset=["start", "end"])

def chr_len_mb(repeat_blast_dir, chr_name):
    # uses your per-chrom length files: chr7.txt, chr8.txt, etc.
    p = os.path.join(repeat_blast_dir, f"{chr_name}.txt")
    df = pd.read_csv(p, sep=r"\s+", header=None)  # tabs OR spaces
    bp = int(df.iloc[0, 1]) if df.shape[1] >= 2 else int(df.iloc[0, 0])
    return bp / 1e6

def add_track(ax, df, color, y, h, z=10, alpha=0.95):
    if df.empty:
        return
    for r in df.itertuples(index=False):
        x0 = float(r.start) / 1e6
        w  = (float(r.end) - float(r.start)) / 1e6
        if w <= 0:
            continue
        ax.add_patch(Rectangle((x0, y), w, h, facecolor=color, edgecolor="none",
                               alpha=alpha, zorder=z))

def plot_panel(ax, sample, chr_name, root_pattern, subdir, tick_step=10):
    sample_dir = root_pattern.format(sample=sample)
    base_dir = os.path.join(sample_dir, subdir)  # .../repeat_blast

    tr1   = read_bed(os.path.join(base_dir, f"intervals_{sample}_TR1.bed"))
    knob  = read_bed(os.path.join(base_dir, f"intervals_{sample}_knob180.bed"))
    centc = read_bed(os.path.join(base_dir, f"intervals_{sample}_CentC.bed"))
    agp   = read_agp(os.path.join(base_dir, f"{chr_name}.agp"))

    L = chr_len_mb(base_dir, chr_name)

    # geometry (same idea as your working script)
    bar_height = 0.4
    y_center = 0.0
    y_bottom = y_center - bar_height/2

    contig_h = bar_height * 0.35
    repeat_h = 0.38
    gap = 0.05
    contig_y = y_bottom - contig_h - gap
    repeat_y = y_center - repeat_h/2

    # backbone
    ax.add_patch(Rectangle((0, y_bottom), L, bar_height,
                           facecolor="#1673A3", edgecolor="none", zorder=1))

    # contigs (AGP)
    if not agp.empty:
        for r in agp.itertuples(index=False):
            x0 = float(r.start)/1e6
            w  = (float(r.end)-float(r.start))/1e6
            if w <= 0:
                continue
            ax.add_patch(Rectangle((x0, contig_y), w, contig_h,
                                   facecolor=COLORS["Contig"], edgecolor="black",
                                   linewidth=0.35, zorder=2))

    # repeats
    add_track(ax, tr1,   COLORS["TR1"],     repeat_y, repeat_h, z=10)
    add_track(ax, knob,  COLORS["knob180"], repeat_y, repeat_h, z=10)
    add_track(ax, centc, COLORS["CentC"],   repeat_y, repeat_h, z=10)

    ax.set_xlim(0, L * 1.02)
    ax.set_ylim(contig_y - 0.15, y_bottom + bar_height + 0.15)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0, L + tick_step, tick_step))
    ax.set_title(sample, fontsize=10)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)          # linename.lst
    ap.add_argument("--chr", required=True)              # chr1/chr2/chr7/chr8
    ap.add_argument("--root_pattern", default="ragtag_{sample}_v3")
    ap.add_argument("--subdir", default="repeat_blast")
    ap.add_argument("--outdir", default="chr_plots")
    ap.add_argument("--tick_step", type=int, default=10)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    with open(args.samples) as f:
        samples = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]
    if len(samples) != 8:
        raise ValueError(f"Expected 8 samples, found {len(samples)} in {args.samples}")

    os.makedirs(args.outdir, exist_ok=True)

    fig, axes = plt.subplots(4, 2, figsize=(14, 10))
    axes = axes.flatten()

    for i, sample in enumerate(samples):
        plot_panel(axes[i], sample, args.chr, args.root_pattern, args.subdir, tick_step=args.tick_step)

        # show x labels only on bottom row
        if i < 6:
            axes[i].set_xticklabels([])
            axes[i].set_xlabel("")
        else:
            axes[i].set_xlabel("Position (Mb)")

    fig.suptitle(f"{args.chr} (8 individuals; 2×4 grid)", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out_pdf = os.path.join(args.outdir, f"{args.chr}.8samples.grid.pdf")
    out_png = os.path.join(args.outdir, f"{args.chr}.8samples.grid.png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", out_pdf)
    print("Saved:", out_png)

if __name__ == "__main__":
    main()
