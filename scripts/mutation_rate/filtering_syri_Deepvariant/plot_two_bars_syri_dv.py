#!/usr/bin/env python3
"""One exact-size indel chart: filtered SyRI and supplied DeepVariant VCF.

Run beside indel_filter_steps_by_column.py from the directory containing
vcf_to_compare. Includes all observed sizes on the x-axis. Calls with >=10
consecutive Ns in REF or ALT are excluded from both series.
"""

import re
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator

import indel_filter_steps_by_column as s


EXTRA_DV_VCF = (
    s.ROOT / "filtered_0.95hom_40-60her_dv"
    / "private_biallelic_indels.all_sizes.exactly_one_line.chr1-10.vcf.gz"
)
N_RUN = re.compile(r"N{10,}", re.IGNORECASE)


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
    """Count the supplied DV VCF without applying the SyRI filters."""
    counts = Counter()
    for item in s.call_records("deepvariant", EXTRA_DV_VCF):
        if item[0] == "names":
            continue
        _, _sample, _chrom, _pos, ref, alt, _group = item
        if valid_alleles(ref, alt):
            counts[abs(len(alt) - len(ref))] += 1
    return counts


def main():
    s.REPEAT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"Output folder: {s.REPEAT_ROOT.resolve()}", flush=True)
    syri = filtered_syri_counts()
    dv = supplied_dv_counts()
    largest = max((size for c in (syri, dv) for size in c if size >= 1),
                  default=1)

    fig, ax = plt.subplots(figsize=(19, 7))
    for counts, label, color, offset in (
        (syri, "SyRI filtered", "#D95F02", -0.20),
        (dv, "DeepVariant supplied VCF", "#CC79A7", 0.20),
    ):
        sizes = sorted(size for size, n in counts.items() if size >= 1 and n)
        # Small indels have fixed-width adjacent bars. Above 10 bp, scale
        # widths with size so singleton calls are visible on the log axis.
        positions = [
            size + offset if size <= 10 else size * (1 + offset * 0.10)
            for size in sizes
        ]
        widths = [0.39 if size <= 10 else 0.035 * size for size in sizes]
        ax.bar(positions, [counts[size] for size in sizes], width=widths,
               color=color, label=f"{label} ({sum(counts.values()):,} calls)")
        print(f"{label}: {sum(counts.values()):,} indel calls; "
              f"largest = {max(sizes, default=0):,} bp", flush=True)

    ax.set_xscale("symlog", linthresh=10, base=10)
    ax.set_xlim(0.5, max(largest * 1.06, 1.5))
    ticks = [1, 2, 5, 10]
    for decade in (10 ** n for n in range(1, 12)):
        ticks.extend(value for value in (2 * decade, 5 * decade, 10 * decade)
                     if value <= largest)
    ax.xaxis.set_major_locator(FixedLocator(sorted(set(ticks))))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax.tick_params(axis="x", labelrotation=40)
    ax.set_yscale("symlog", linthresh=10, base=10)
    ax.set_xlabel("Absolute indel size (bp; linear through 10, log above)")
    ax.set_ylabel("Indel calls per exact size (linear through 10, log above)")
    ax.set_title("Private indels of all sizes: SyRI vs supplied DeepVariant "
                 "(N10 excluded)", fontweight="bold")
    ax.grid(axis="y", which="both", alpha=0.18)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)
    fig.tight_layout()
    output = s.REPEAT_ROOT / "syri_vs_supplied_DV_all_indel_sizes.no_N10.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved: {output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
