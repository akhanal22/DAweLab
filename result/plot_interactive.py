#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COLORS = {
    "Backbone": "#1673A3",
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
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
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
    agp["end"] = pd.to_numeric(agp["end"], errors="coerce")
    return agp.dropna(subset=["start", "end"])

def chr_len_mb(base_dir, chr_name):
    p = os.path.join(base_dir, f"{chr_name}.txt")
    df = pd.read_csv(p, sep=r"\s+", header=None)
    bp = int(df.iloc[0, 1]) if df.shape[1] >= 2 else int(df.iloc[0, 0])
    return bp / 1e6

def add_track(ax, df, y, h, color, alpha=0.95, z=10):
    if df.empty:
        return
    for r in df.itertuples(index=False):
        x0 = float(r.start) / 1e6
        w  = (float(r.end) - float(r.start)) / 1e6
        if w <= 0:
            continue
        ax.add_patch(Rectangle((x0, y), w, h,
                               facecolor=color,
                               edgecolor="none",
                               alpha=alpha,
                               zorder=z))

def add_contigs(ax, agp, y, h, color, z=5):
    if agp.empty:
        return
    for r in agp.itertuples(index=False):
        x0 = float(r.start) / 1e6
        w  = (float(r.end) - float(r.start)) / 1e6
        if w <= 0:
            continue
        ax.add_patch(Rectangle((x0, y), w, h,
                               facecolor=color,
                               edgecolor="black",
                               linewidth=0.3,
                               zorder=z))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--chr", required=True)
    ap.add_argument("--root_pattern", default="ragtag_{sample}_v3")
    ap.add_argument("--subdir", default="repeat_blast")
    ap.add_argument("--tick_step", type=int, default=10)
    ap.add_argument("--no_agp", action="store_true")
    args = ap.parse_args()

    with open(args.samples) as f:
        samples = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith("#")]

    # determine chromosome length from first sample
    first_dir = os.path.join(args.root_pattern.format(sample=samples[0]), args.subdir)
    L = chr_len_mb(first_dir, args.chr)

    fig_h = max(4, 0.6 * len(samples))
    fig, ax = plt.subplots(figsize=(14, fig_h))

    row_gap = 0.35
    backbone_h = 0.20
    repeat_h = 0.22
    contig_h = 0.10

    for i, sample in enumerate(samples):
        y0 = (len(samples) - 1 - i) * (backbone_h + row_gap)

        base_dir = os.path.join(args.root_pattern.format(sample=sample), args.subdir)

        tr1   = read_bed(os.path.join(base_dir, f"intervals_{sample}_TR1.bed"))
        knob  = read_bed(os.path.join(base_dir, f"intervals_{sample}_knob180.bed"))
        centc = read_bed(os.path.join(base_dir, f"intervals_{sample}_CentC.bed"))

        # backbone
        ax.add_patch(Rectangle((0, y0), L, backbone_h,
                               facecolor=COLORS["Backbone"],
                               edgecolor="none",
                               zorder=1))

        # contigs
        if not args.no_agp:
            agp = read_agp(os.path.join(base_dir, f"{args.chr}.agp"))
            add_contigs(ax, agp, y0 - contig_h - 0.03, contig_h,
                        COLORS["Contig"], z=2)

        # repeats
        add_track(ax, tr1,   y0, repeat_h, COLORS["TR1"])
        add_track(ax, knob,  y0, repeat_h, COLORS["knob180"])
        add_track(ax, centc, y0, repeat_h, COLORS["CentC"])

        ax.text(-0.02 * L, y0 + backbone_h/2, sample,
                ha="right", va="center", fontsize=10)

    ax.set_xlim(0, L * 1.02)
    ax.set_yticks([])
    ax.set_xlabel("Position (Mb)")
    ax.set_xticks(np.arange(0, L + args.tick_step, args.tick_step))
    ax.set_title(f"{args.chr} stacked features across samples", fontsize=14)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
