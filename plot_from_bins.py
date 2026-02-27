#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt

def read_bins(tsv_path):
    # TSV: midpoint_kbp \t mean_depth
    x = []
    y = []
    with open(tsv_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            a, b = line.rstrip("\n").split("\t")
            x.append(float(a))
            y.append(float(b))
    return x, y

def rolling_mean(x, y, win_bins, step_bins):
    """
    Rolling mean over binned depths.
    x: midpoints in kbp (one per bin, evenly spaced)
    y: mean depth per bin
    """
    if win_bins <= 0 or step_bins <= 0:
        raise ValueError("win_bins and step_bins must be positive")

    outx, outy = [], []
    n = len(y)
    if n < win_bins:
        return outx, outy

    # initial window sum
    s = sum(y[:win_bins])
    i = 0
    while i + win_bins <= n:
        mid_index = i + win_bins // 2
        outx.append(x[mid_index])
        outy.append(s / win_bins)

        # advance
        next_i = i + step_bins
        if next_i + win_bins > n:
            break

        # update sum efficiently by removing/adding step_bins elements
        for k in range(i, next_i):
            s -= y[k]
        for k in range(i + win_bins, next_i + win_bins):
            s += y[k]

        i = next_i

    return outx, outy

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary_bins", required=True, help="primary binned TSV (mid_kbp\\tdepth)")
    ap.add_argument("--secondary_bins", required=True, help="secondary binned TSV (mid_kbp\\tdepth)")
    ap.add_argument("--bin_bp", type=int, default=5000, help="bin size in bp (default 5000)")
    ap.add_argument("--window_bp", type=int, default=20000, help="sliding window size in bp (default 20000)")
    ap.add_argument("--step_bp", type=int, default=5000, help="sliding step size in bp (default 5000)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--out_tsv", default=None, help="optional output TSV of rolling means")
    ap.add_argument("--maxy", type=float, default=None)
    args = ap.parse_args()

    # Convert bp params to bins
    if args.window_bp % args.bin_bp != 0 or args.step_bp % args.bin_bp != 0:
        raise ValueError("window_bp and step_bp must be multiples of bin_bp")

    win_bins = args.window_bp // args.bin_bp
    step_bins = args.step_bp // args.bin_bp

    xp, yp = read_bins(args.primary_bins)
    xs, ys = read_bins(args.secondary_bins)

    # Rolling mean on binned data
    rxp, ryp = rolling_mean(xp, yp, win_bins, step_bins)
    rxs, rys = rolling_mean(xs, ys, win_bins, step_bins)

    title = args.title

    # Optional TSV
    if args.out_tsv:
        n = min(len(rxp), len(rxs), len(ryp), len(rys))
        with open(args.out_tsv, "w") as out:
            out.write("midpoint_kbp\tprimary_mean_depth\tsecondary_mean_depth\n")
            for i in range(n):
                out.write(f"{rxp[i]}\t{ryp[i]}\t{rys[i]}\n")

    # Plot
    plt.figure(figsize=(14, 3.5))
    plt.plot(rxp, ryp, linewidth=1.0, label="Primary", color="black")
    plt.plot(rxs, rys, linewidth=1.0, label="Secondary", color="red")
    plt.xlabel("Position (kbp)")
    plt.ylabel("Mean depth (sliding on bins)")
    plt.title(title)
    plt.legend(frameon=False)
    # Force y-axis range and ticks
    plt.ylim(-5, 200)
    plt.yticks([0, 50, 100, 150, 200])
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=200)
    plt.close()

if __name__ == "__main__":
    main()
