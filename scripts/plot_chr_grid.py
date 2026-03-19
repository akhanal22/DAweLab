#!/usr/bin/env python3
import argparse
import math
import os
from collections import defaultdict

import matplotlib.pyplot as plt


def iter_depth_bins(depth_path, chrom, bin_bp):
    """Read samtools depth (chr pos depth), bin by bin_bp, return x(Mb), y(mean depth per bin)."""
    sums = defaultdict(float)
    cnts = defaultdict(int)

    with open(depth_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            c, pos_s, dep_s = line.rstrip("\n").split("\t")
            if c != chrom:
                continue
            pos = int(pos_s)
            dep = float(dep_s)
            b = (pos - 1) // bin_bp
            sums[b] += dep
            cnts[b] += 1

    if not sums:
        return [], []

    bins_sorted = sorted(sums.keys())
    x = []
    y = []
    for b in bins_sorted:
        mid_bp = b * bin_bp + bin_bp / 2.0
        x.append(mid_bp / 1e6)         # Mb
        y.append(sums[b] / cnts[b])
    return x, y


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="File with sample names, one per line")
    ap.add_argument("--chr", required=True, help="Chromosome name, e.g. chr4_RagTag or chr4")
    ap.add_argument("--depth_pattern", default="./{sample}/depth/{sample}.depth",
                    help="Depth file pattern. Use {sample}.")
    ap.add_argument("--bin_bp", type=int, default=100000, help="Bin size in bp (default 100kb)")
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--title", default=None)
    ap.add_argument("--sharey", action="store_true", help="Use same y-scale across all panels")
    args = ap.parse_args()

    # Read sample list
    samples = [s.strip() for s in open(args.samples) if s.strip()]
    n = len(samples)
    if n == 0:
        raise SystemExit("ERROR: samples file is empty")

    # Fixed layout: 2 columns, rows computed
    ncols = 2
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 10),
                             sharex=True, sharey=args.sharey)
    axes = axes.flatten() if n > 1 else [axes]

    # Track global y max if sharey requested (optional)
    global_ymax = 0.0

    for i, sample in enumerate(samples):
        ax = axes[i]
        depth_path = args.depth_pattern.format(sample=sample)

        if not os.path.isfile(depth_path) or os.path.getsize(depth_path) == 0:
            ax.text(0.5, 0.5, f"Missing depth\n{sample}", ha="center", va="center")
            ax.set_title(sample)
            continue

        x, y = iter_depth_bins(depth_path, args.chr, args.bin_bp)
        if not x:
            ax.text(0.5, 0.5, f"No data for {args.chr}\n{sample}", ha="center", va="center")
            ax.set_title(sample)
            continue

        ax.plot(x, y, linewidth=0.8)
        ax.set_title(sample)
        global_ymax = max(global_ymax, max(y))

    # Turn off any unused panels
    for j in range(n, len(axes)):
        axes[j].axis("off")

    # Labels
    for ax in axes[:n]:
        ax.set_ylabel("Mean depth per bin")
    for ax in axes[max(0, n - ncols):n]:
        ax.set_xlabel("Position (Mb)")

    if args.title:
        fig.suptitle(args.title, y=0.995)
    else:
        fig.suptitle(f"{args.chr} depth (binned)", y=0.995)

    if args.sharey and global_ymax > 0:
        for ax in axes[:n]:
            ax.set_ylim(0, global_ymax * 1.05)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.out_png, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
