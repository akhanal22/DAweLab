#!/usr/bin/env python3
"""Plot filtered SyRI indels from 50 bp through the largest retained call.

Run from the directory containing vcf_to_compare with
indel_filter_steps_by_column.py available for import. From 50 to 100 bp,
10 bp ranges have equal spacing. Beyond 100 bp, each decade is divided into
five logarithmically spaced ranges, so 50–100 bp and every subsequent
tenfold span occupy the same horizontal width. Bars have equal visual width.
All calls with 10 consecutive Ns are excluded; calls above 100,000 bp are
also excluded if REF or ALT has three consecutive Ns.
"""

from collections import Counter
from math import floor, log10
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import indel_filter_steps_by_column as s


MIN_SIZE = 50
MIN_AXIS_END = 100000
LOG_BINS_PER_DECADE = 5
N_RUN = re.compile(r"N{10,}", re.IGNORECASE)
LARGE_N_RUN = re.compile(r"N{3,}", re.IGNORECASE)
SYRI_COLOR = "#D55E00"  # same orange as the two-caller figure


def valid_alleles(ref, alt):
    return not (N_RUN.search(ref) or N_RUN.search(alt))


def filtered_syri_counts():
    """Apply the same per-sample 500 kb and repeat filters as before."""
    cg_path = s.REPEAT_ROOT / "reference.CG_homopolymers.min5.bed"
    cg_core = s.read_bed(cg_path)
    cg_runs = s.read_core_runs(cg_path)
    cg_flank = s.read_bed(
        s.REPEAT_ROOT / "reference.CG_homopolymers.min5.plus1bp.bed"
    )
    di = s.read_bed(
        s.REPEAT_ROOT / "reference.perfect_dinucleotide_repeats.min3copies.bed"
    )
    tri = s.read_bed(
        s.REPEAT_ROOT / "reference.perfect_trinucleotide_repeats.min3copies.bed"
    )

    counts = Counter()
    bad_bins = {}
    for item in s.call_records("syri", s.INPUTS["syri"]):
        if item[0] == "names":
            bad_bins = {
                sample: s.read_bed(
                    s.WINDOW_ROOT / "syri" / "bad_500kb_bins"
                    / f"{sample}.bad_500kb_bins.bed"
                )
                for sample in item[1]
            }
            continue
        _, sample, chrom, pos, ref, alt, _group = item
        size = abs(len(alt) - len(ref))
        if size < MIN_SIZE or not valid_alleles(ref, alt):
            continue
        if size > MIN_AXIS_END and (LARGE_N_RUN.search(ref) or LARGE_N_RUN.search(alt)):
            continue
        if s.touches(bad_bins[sample], chrom, pos - 1, pos):
            continue
        start, end, inserted = s.event_interval(pos, ref, alt)
        if s.in_cg_core(cg_core, cg_runs, chrom, start, end, inserted):
            continue
        if s.touches(cg_flank, chrom, start, end):
            continue
        if size <= s.MAX_REPEAT_INDEL_BP and s.touches(di, chrom, start, end):
            continue
        if size <= s.MAX_REPEAT_INDEL_BP and s.touches(tri, chrom, start, end):
            continue
        counts[size] += 1
    return counts


def x_position(size):
    """Linear through 100 bp; logarithmic beyond 100 bp."""
    if size <= 100:
        return (size - 50.0) / 10.0
    return 5.0 + LOG_BINS_PER_DECADE * log10(size / 100.0)


def bin_index(size):
    """50–59, ..., 90–99, then five logarithmic ranges per decade."""
    if size < 100:
        return (size - MIN_SIZE) // 10
    return 5 + floor(LOG_BINS_PER_DECADE * log10(size / 100.0))


def x_ticks(axis_end):
    ticks = [50, 60, 70, 80, 90, 100]
    decade = 1000
    while decade <= axis_end:
        ticks.append(decade)
        decade *= 10
    return ticks


def main():
    s.REPEAT_ROOT.mkdir(parents=True, exist_ok=True)
    counts = filtered_syri_counts()
    sizes = sorted(counts)
    if not sizes:
        raise SystemExit(f"No filtered SyRI indels of at least {MIN_SIZE} bp")

    binned = Counter()
    for size, n in counts.items():
        binned[bin_index(size)] += n

    axis_end = max(MIN_AXIS_END, sizes[-1])

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.bar(
        [index + 0.5 for index in sorted(binned)],
        [binned[index] for index in sorted(binned)],
        width=0.88, color=SYRI_COLOR, linewidth=0,
        label=f"SyRI filtered ({sum(counts.values()):,} calls)",
    )
    ticks = x_ticks(axis_end)
    ax.set_xticks([x_position(size) for size in ticks],
                  [f"{size:,}" for size in ticks])
    ax.set_xlim(-0.05, max(x_position(axis_end), max(binned) + 1) + 0.05)
    ax.set_yscale("symlog", linthresh=10, base=10)
    ax.tick_params(axis="both", labelsize=15)
    ax.set_xlabel("Indel size (bp; linear 50–100, logarithmic above 100)",
                  fontsize=19, labelpad=14)
    ax.set_ylabel("Indel calls per size range (symlog scale)",
                  fontsize=19, labelpad=12)
    ax.set_title("Private SyRI indels ≥50 bp (N10 excluded)",
                 fontsize=20, fontweight="bold", pad=16)
    ax.grid(axis="y", which="both", alpha=0.18)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=15)
    fig.tight_layout()

    output = s.REPEAT_ROOT / "filtered_syri_indels_50bp_onward.dynamic_axis.no_N10_large_no_N3.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"SyRI: {sum(counts.values()):,} calls ≥{MIN_SIZE} bp; "
          f"{sum(n for size, n in counts.items() if size > MIN_AXIS_END):,} "
          f"above {MIN_AXIS_END:,} bp; largest = {sizes[-1]:,} bp", flush=True)
    print(f"Saved: {output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
