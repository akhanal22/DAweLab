#!/usr/bin/env python3
"""Plot exact 1–50 bp indel sizes from the supplied DeepVariant VCF."""

from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import indel_filter_steps_by_column as s


VCF = (
    s.ROOT
    / "filtered_0.95hom_40-60her_dv"
    / "private_biallelic_indels.all_sizes.exactly_one_line.chr1-10.vcf.gz"
)
OUTPUT = s.REPEAT_ROOT / "supplied_DV_hom0.95_1-50bp_blue.png"


def main():
    counts = Counter()

    for item in s.call_records("deepvariant", VCF):
        if item[0] == "names":
            continue

        _, sample, chrom, pos, ref, alt, group = item
        size = abs(len(alt) - len(ref))
        if 1 <= size <= 50:
            counts[size] += 1

    total = sum(counts.values())
    sizes = range(1, 51)

    fig, ax = plt.subplots(figsize=(19, 7))
    ax.bar(
        sizes,
        [counts[size] for size in sizes],
        width=0.8,
        color="#0072B2",  # Blue
    )
    ax.set_xlim(0.5, 50.5)
    ax.set_xticks([1] + list(range(5, 51, 5)))
    ax.set_yscale("symlog", linthresh=10, base=10)
    ax.set_xlabel("Indel size (bp; exact sizes 1–50)")
    ax.set_ylabel("Calls per exact size (linear through 10, log above)")
    ax.set_title(
        f"Supplied DeepVariant VCF: 1–50 bp ({total:,} calls)",
        fontweight="bold",
    )
    ax.grid(axis="y", which="both", alpha=0.18)
    ax.set_axisbelow(True)
    fig.tight_layout()

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=200)
    plt.close(fig)
    print(f"Saved: {OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
