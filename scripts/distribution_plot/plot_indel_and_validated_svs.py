#!/usr/bin/env python3
"""Side-by-side size distribution of indels/SVs called by DeepVariant and SyRI.

Both callers are shown over the full size range:
  - DeepVariant: all private biallelic indels (any size)
  - SyRI: all INS/DEL <=50 bp (no read validation exists for these), and only
    validated SVs >50 bp (depth ratio >= 0.70, all junction checks PASS, as in
    plot_j702_variant_ideogram.py)
Variants with N in REF/ALT are excluded.

x axis: one bar per bp from 1 to 50 bp. Above 50 bp the size ranges double
(50-100, 100-200, 200-400, 400-800, ...) and each doubling range is split into
SUB equal-width bins (default 5: 50-60, 60-70, ..., 100-120, 120-140, ...).
Bins include their upper bound (e.g. 51-60). Every bin gets the same width on
the axis, and within it
DeepVariant (blue, left) and SyRI (deep orange, right) bars sit side by side with
identical widths. y axis: count, log scale. Colours: Okabe-Ito.
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_j702_variant_ideogram import CHROMS, read_vcf, validated_sv_keys

BLUE, ORANGE = "#0072B2", "#D55E00"  # Okabe-Ito blue and vermillion (deep orange)
LIN_MAX = 50


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--syri", default="syri.after_tri_1_50.vcf.gz")
    ap.add_argument("--indels", default="private_biallelic_indels.all_sizes.exactly_one_line.chr1-10.vcf.gz")
    ap.add_argument("--ins-val", default="insertion_validation_with_cigar_with500kbfilteredvcfagainstj702.tsv")
    ap.add_argument("--del-val", default="deletion_validation_withcigar_for500kbfilteredvcfagainstj702.tsv")
    ap.add_argument("--depth-ratio", type=float, default=0.70)
    ap.add_argument("--out", default="indel_sv_size_distribution")
    a = ap.parse_args()

    _, indels = read_vcf(a.indels)
    _, syri = read_vcf(a.syri)
    val = validated_sv_keys(a.ins_val, a.del_val, a.depth_ratio)
    size = lambda r, al: abs(len(al) - len(r))
    hasN = lambda r, al: "N" in (r + al).upper()

    dv = np.array([size(r, al) for c, p, vid, r, al, smp in indels
                   if c in CHROMS and size(r, al) > 0 and not hasN(r, al)])
    sy = np.array([size(r, al) for c, p, vid, r, al, smp in syri
                   if c in CHROMS and not al.startswith("<") and size(r, al) > 0 and not hasN(r, al)
                   and (size(r, al) <= LIN_MAX or any((c, p, vid, s) in val for s in smp))])
    for name, x in (("DeepVariant", dv), ("SyRI", sy)):
        print(f"{name}: n={len(x)}  <=50 bp: {(x <= LIN_MAX).sum()}  >50 bp: {(x > LIN_MAX).sum()}  max {x.max():,} bp")

    # bins: one bar per bp for 1..50; above 50 bp, doubling ranges
    # (50-100, 100-200, 200-400, ...) each split into SUB equal-width bins
    SUB = 5

    def doubling(x):  # d such that 50*2^d < x <= 50*2^(d+1)
        d = 0
        while x > LIN_MAX * 2 ** (d + 1):
            d += 1
        return d

    ndbl = doubling(max(dv.max(), sy.max())) + 1
    nbins = LIN_MAX + SUB * ndbl

    def to_bin(xs):
        out = []
        for x in xs:
            if x <= LIN_MAX:
                out.append(x - 1)
            else:
                d = doubling(x)
                lo, step = LIN_MAX * 2 ** d, LIN_MAX * 2 ** d / SUB
                out.append(LIN_MAX + SUB * d + int(np.ceil((x - lo) / step)) - 1)
        return np.array(out, dtype=int)

    cdv = np.bincount(to_bin(dv), minlength=nbins)
    csy = np.bincount(to_bin(sy), minlength=nbins)
    with open(f"{a.out}.bins.tsv", "w") as fh:
        fh.write("bin_bp\tDeepVariant\tSyRI\n")
        for i in range(nbins):
            if i < LIN_MAX:
                lab = str(i + 1)
            else:
                d, k = divmod(i - LIN_MAX, SUB)
                lo, step = LIN_MAX * 2 ** d, LIN_MAX * 2 ** d // SUB
                lab = f"{lo + k * step + 1}-{lo + (k + 1) * step}"
            fh.write(f"{lab}\t{cdv[i]}\t{csy[i]}\n")

    fig, ax = plt.subplots(figsize=(17, 5))
    centers = np.arange(nbins) + 0.5
    w = 0.42
    ax.bar(centers - w / 2, cdv, width=w, color=BLUE, ec="none",
           label=f"DeepVariant indels (n={len(dv):,}; {(dv > LIN_MAX).sum()} >50 bp)")
    ax.bar(centers + w / 2, csy, width=w, color=ORANGE, ec="none",
           label=f"SyRI indels/SVs (n={len(sy):,}; ≤50 bp all, >50 bp validated only)")
    ax.set_yscale("log")
    ax.axvline(LIN_MAX, color="#555555", lw=0.8, ls="--")
    ymax = max(cdv.max(), csy.max())
    ax.text(LIN_MAX / 2, ymax * 2.5, "1 bp per bar", ha="center", fontsize=9, color="#555555")
    ax.text(LIN_MAX + SUB * ndbl / 2, ymax * 2.5,
            f"{SUB} bars per doubling (50–100, 100–200, 200–400, ...)", ha="center",
            fontsize=9, color="#555555")

    def u_of(x):  # axis position of size x above 50 bp (interpolated inside its doubling)
        d = doubling(x)
        lo = LIN_MAX * 2 ** d
        return LIN_MAX + SUB * d + (x - lo) / (lo / SUB)

    lin_ticks = [1, 10, 20, 30, 40]
    top = LIN_MAX * 2 ** ndbl
    big = [t for t in (100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000) if t <= top]
    fmt = lambda v: f"{v // 1000:,}k" if v >= 1000 else str(v)
    ax.set_xticks([t - 0.5 for t in lin_ticks] + [LIN_MAX] + [u_of(t) for t in big])
    ax.set_xticklabels([str(t) for t in lin_ticks] + ["50"] + [fmt(t) for t in big])
    # minor ticks at every bin edge above 50 bp
    ax.set_xticks([LIN_MAX + i for i in range(SUB * ndbl + 1)], minor=True)
    for lab in ax.get_xticklabels()[len(lin_ticks) + 1:]:
        lab.set_rotation(90)
    for t in big:
        ax.axvline(u_of(t), color="#dddddd", lw=0.6, zorder=0)
    ax.set_xlim(-0.5, nbins + 0.5)
    ax.set_ylim(0.7, ymax * 6)
    ax.set_xlabel("Indel / SV size (bp)")
    ax.set_ylabel("Number of variants")
    ax.grid(axis="y", color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0.20, 0.88))
    ax.set_title("Size distribution of private indels and SVs: DeepVariant vs SyRI", loc="left", fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=300)


if __name__ == "__main__":
    main()
