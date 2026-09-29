#!/usr/bin/env python3
"""Plot bars and a dot alternative on a true logarithmic x-axis.

Run from the directory containing vcf_to_compare, with
indel_filter_steps_by_column.py available for import. The x-axis is log10
starting at 1 bp, so even sizes 1–9 get progressively closer together.
Each exact size gets a pair of approximately equal visual-width bars.
Large nearby sizes can overlap visually. Calls with
at least 10 consecutive Ns in REF or ALT are excluded from both series.
The dot plot shows the same exact-size counts without wide, overlapping bars.
"""

from collections import Counter
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import indel_filter_steps_by_column as s


EXTRA_DV_VCF = (
    s.ROOT / "filtered_0.95hom_40-60her_dv"
    / "private_biallelic_indels.all_sizes.exactly_one_line.chr1-10.vcf.gz"
)
N_RUN = re.compile(r"N{10,}", re.IGNORECASE)
SYRI_COLOR = "#D55E00"  # orange
DV_COLOR = "#0072B2"    # blue


def valid_alleles(ref, alt):
    return not (N_RUN.search(ref) or N_RUN.search(alt))


def filtered_syri_counts():
    """Use the existing per-sample 500 kb and repeat filters."""
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
        if s.touches(bad_bins[sample], chrom, pos - 1, pos):
            continue
        start, end, inserted = s.event_interval(pos, ref, alt)
        if s.in_cg_core(cg_core, cg_runs, chrom, start, end, inserted):
            continue
        if s.touches(cg_flank, chrom, start, end):
            continue
        size = abs(len(alt) - len(ref))
        if size <= s.MAX_REPEAT_INDEL_BP and s.touches(di, chrom, start, end):
            continue
        if size <= s.MAX_REPEAT_INDEL_BP and s.touches(tri, chrom, start, end):
            continue
        if valid_alleles(ref, alt):
            counts[size] += 1
    return counts


def supplied_dv_counts():
    """Count supplied DV calls without adding the SyRI-specific filters."""
    counts = Counter()
    for item in s.call_records("deepvariant", EXTRA_DV_VCF):
        if item[0] == "names":
            continue
        _, _sample, _chrom, _pos, ref, alt, _group = item
        if valid_alleles(ref, alt):
            counts[abs(len(alt) - len(ref))] += 1
    return counts


def x_ticks():
    """Only the requested size labels, positioned on the true log axis."""
    return [1, 10, 1000, 10000, 100000]


def save_dot_plot(syri, dv, largest):
    """Show only observed exact sizes; gaps mean no retained calls."""
    fig, ax = plt.subplots(figsize=(20, 8))
    for counts, label, color, marker, scale in (
        (syri, "SyRI filtered", SYRI_COLOR, "o", 0.98),
        (dv, "DeepVariant supplied VCF", DV_COLOR, "s", 1.02),
    ):
        sizes = sorted(size for size, n in counts.items() if size >= 1 and n)
        ax.scatter(
            [size * scale for size in sizes],
            [counts[size] for size in sizes],
            s=36, marker=marker, color=color, alpha=0.85,
            edgecolors="white", linewidths=0.35,
            label=f"{label} ({sum(counts[size] for size in sizes):,} calls)",
            zorder=3,
        )

    ticks = x_ticks()
    ax.set_xscale("log", base=10)
    ax.set_xticks(ticks,
                  [f"{size:,}" for size in ticks])
    ax.set_xlim(0.80, max(100000, largest) * 1.20)
    ax.set_yscale("log", base=10)
    ax.set_ylim(bottom=0.8)
    ax.tick_params(axis="both", labelsize=15)
    ax.set_xlabel("Exact indel size (bp; log10 scale)",
                  fontsize=19, labelpad=14)
    ax.set_ylabel("Calls at each observed exact size (log10 scale)",
                  fontsize=19, labelpad=12)
    ax.set_title("Private indel sizes: SyRI vs supplied DeepVariant "
                 "(N10 excluded)", fontsize=20, fontweight="bold", pad=16)
    ax.grid(axis="y", which="both", alpha=0.18)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=15)
    fig.tight_layout()
    output = (s.REPEAT_ROOT /
              "syri_vs_supplied_DV_all_indel_sizes_dotplot.png")
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved: {output.resolve()}", flush=True)


def main():
    s.REPEAT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"Output folder: {s.REPEAT_ROOT.resolve()}", flush=True)
    syri = filtered_syri_counts()
    dv = supplied_dv_counts()
    largest = max(
        (size for counts in (syri, dv) for size in counts if size >= 1),
        default=1,
    )
    fig, ax = plt.subplots(figsize=(20, 8))
    for counts, label, color, scale in (
        (syri, "SyRI filtered", SYRI_COLOR, 0.93),
        (dv, "DeepVariant supplied VCF", DV_COLOR, 1.07),
    ):
        sizes = sorted(size for size, n in counts.items() if size >= 1 and n)
        positions = [size * scale for size in sizes]
        values = [counts[size] for size in sizes]
        ax.bar(
            positions, values, width=[size * 0.12 for size in sizes], color=color,
            label=f"{label} ({sum(values):,} calls)", linewidth=0,
        )
        print(
            f"{label}: {sum(values):,} indel calls; "
            f"largest = {max(sizes, default=0):,} bp",
            flush=True,
        )

    ticks = x_ticks()
    ax.set_xscale("log", base=10)
    ax.set_xticks(ticks,
                  [f"{size:,}" for size in ticks])
    ax.set_xlim(0.80, max(100000, largest) * 1.20)
    ax.set_yscale("log", base=10)
    ax.set_ylim(bottom=0.8)
    ax.tick_params(axis="both", labelsize=15)
    ax.set_xlabel(
        "Exact indel size (bp; log10 scale)",
        fontsize=19, labelpad=14,
    )
    ax.set_ylabel("Indel calls per exact size (log10 scale)",
                  fontsize=19, labelpad=12)
    ax.set_title(
        "Private indels of all sizes: SyRI vs supplied DeepVariant "
        "(N10 excluded)",
        fontsize=20, fontweight="bold", pad=16,
    )
    ax.grid(axis="y", which="both", alpha=0.18)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=15)
    fig.tight_layout()
    output = (
        s.REPEAT_ROOT
        / "syri_vs_supplied_DV_all_indel_sizes.no_N10.true_log_bars.five_ticks.png"
    )
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved: {output.resolve()}", flush=True)
    save_dot_plot(syri, dv, largest)


if __name__ == "__main__":
    main()
