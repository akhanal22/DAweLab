#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from bisect import bisect_right
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


# ------------------------------------------------------------
# Tentative J702 centromere coordinates supplied by the user.
# Coordinates are treated as 1-based, inclusive.
#
# For each chromosome, the plot uses the full span between
# the first start and last end, including intervening gaps.
# ------------------------------------------------------------

CENTROMERE_PIECES = [
    ("chr1_RagTag", 137822536, 137958224),
    ("chr1_RagTag", 138938038, 139057024),
    ("chr1_RagTag", 137958470, 138924181),
    ("chr2_RagTag", 96062225, 98111512),
    ("chr3_RagTag", 86349928, 88239902),
    ("chr4_RagTag", 109019215, 110699213),
    ("chr5_RagTag", 104490280, 106310257),
    ("chr6_RagTag", 61044759, 61168614),
    ("chr6_RagTag", 62540512, 63221738),
    ("chr7_RagTag", 56189237, 57071016),
    ("chr7_RagTag", 57124718, 57651903),
    ("chr7_RagTag", 57662303, 57761908),
    ("chr8_RagTag", 50558812, 52448798),
    ("chr9_RagTag", 61572438, 63207422),
    ("chr10_RagTag", 50772346, 52708768),
]


# ------------------------------------------------------------
# SV sizes supplied for the tentative centromere candidates.
#
# These seven events passed the 0.70 validation threshold
# and the other recorded spanning-read checks.
# ------------------------------------------------------------

CENTROMERE_SV_SIZES = [
    1079,  # B73D10 INS23320
    1551,  # B73D2  INS20082
    341,   # B73D2  INS24583
    161,   # B73D1  DEL1573
    6642,  # B73D1  DEL9700
    6640,  # B73D1  DEL9706
    2322,  # B73D7c1 DEL17677
]


# ------------------------------------------------------------
# knob180-associated SV sizes from the user's repeat table.
# Insertions and deletions are combined.
# ------------------------------------------------------------

KNOB180_SV_SIZES = [
    # B73D1
    3325, 360, 180, 3311, 3997,

    # B73D2
    1609, 2179, 1797, 2353, 583,

    # B73D4
    920, 2157, 762, 381, 382, 1977,

    # B73D7c1
    2902, 1460, 920, 772, 1796,

    # B73D10
    1280, 1099, 718, 212, 2377,
    1711, 2155, 2938, 2919, 1839,
]


# CentC-associated SV sizes from the user's repeat table.
CENTC_SV_SIZES = [
    1551,  # B73D2  INS20082
    1079,  # B73D10 INS23320
]


def open_text(path):
    """Open plain text or gzip-compressed text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "rt")


def normalize_chrom(chrom):
    """
    Match chr1 and chr1_RagTag naming when appropriate.
    """
    chrom = chrom.strip()

    if chrom.endswith("_RagTag"):
        return chrom

    if chrom.startswith("chr"):
        suffix = chrom[3:]

        if suffix.isdigit():
            return f"{chrom}_RagTag"

    return chrom


def merge_intervals(intervals):
    """
    Merge overlapping or adjacent zero-based half-open intervals.

    Returns:
        {
            chromosome: (list_of_starts, merged_intervals)
        }
    """

    indexed = {}

    for chrom, regions in intervals.items():
        regions = sorted(regions)
        merged = []

        for start, end in regions:
            if end <= start:
                continue

            if merged and start <= merged[-1][1]:
                old_start, old_end = merged[-1]
                merged[-1] = (old_start, max(old_end, end))
            else:
                merged.append((start, end))

        indexed[chrom] = (
            [start for start, end in merged],
            merged,
        )

    return indexed


def overlaps_index(chrom, start, end, indexed):
    """
    Test whether [start, end) overlaps an indexed interval.
    """
    chrom = normalize_chrom(chrom)

    if chrom not in indexed:
        return False

    starts, intervals = indexed[chrom]

    index = bisect_right(starts, end - 1) - 1

    if index < 0:
        return False

    interval_start, interval_end = intervals[index]

    return interval_start < end and interval_end > start


def read_repeat_intervals(gff_path):
    """
    Extract knob180 and CentC features from the J702 GFF.

    GFF uses 1-based inclusive coordinates. These are converted
    to 0-based, half-open coordinates for interval comparison.
    """

    repeat_intervals = {
        "knob180": defaultdict(list),
        "CentC": defaultdict(list),
    }

    with open_text(gff_path) as infile:
        for line in infile:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            chrom = normalize_chrom(fields[0])
            start = int(fields[3]) - 1
            end = int(fields[4])
            attributes = fields[8].lower()

            if (
                "name=knob180" in attributes
                or "classification=knob/knob180" in attributes
            ):
                repeat_intervals["knob180"][chrom].append(
                    (start, end)
                )

            elif (
                "name=centc" in attributes
                or "classification=cent/centc" in attributes
            ):
                repeat_intervals["CentC"][chrom].append(
                    (start, end)
                )

    return {
        repeat_name: merge_intervals(regions)
        for repeat_name, regions in repeat_intervals.items()
    }


def make_centromere_spans():
    """
    Make one full centromere span per chromosome.

    Includes intervening gaps between the supplied pieces.
    """

    spans = {}

    for chrom, start, end in CENTROMERE_PIECES:
        chrom = normalize_chrom(chrom)

        # Convert from 1-based inclusive to 0-based half-open.
        start0 = start - 1
        end0 = end

        if chrom not in spans:
            spans[chrom] = [start0, end0]
        else:
            spans[chrom][0] = min(spans[chrom][0], start0)
            spans[chrom][1] = max(spans[chrom][1], end0)

    return {
        chrom: tuple(region)
        for chrom, region in spans.items()
    }


def overlaps_centromere(chrom, start, end, spans):
    chrom = normalize_chrom(chrom)

    if chrom not in spans:
        return False

    region_start, region_end = spans[chrom]

    return start < region_end and end > region_start


def read_deepvariant_indels(vcf_path, repeat_index):
    """
    Read biallelic, sequence-resolved private DeepVariant indels.

    Includes sizes 1–50 bp.

    Insertions:
        Test the VCF anchor base against the intervals.

    Deletions:
        Test the deleted reference bases, excluding the anchor.
        Assumes standard left-anchored VCF representation.
    """

    centromere_spans = make_centromere_spans()

    sizes = {
        "centromere": [],
        "knob180": [],
        "CentC": [],
        "genome": [],
    }

    skipped = 0

    with open_text(vcf_path) as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 5:
                skipped += 1
                continue

            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            # Skip multiallelic and symbolic records.
            if (
                "," in alt
                or alt.startswith("<")
                or "[" in alt
                or "]" in alt
                or alt == "*"
            ):
                skipped += 1
                continue

            difference = len(alt) - len(ref)
            indel_size = abs(difference)

            if difference == 0 or not 1 <= indel_size <= 50:
                continue

            if difference > 0:
                # One-base VCF anchor interval.
                query_start = pos - 1
                query_end = pos

            else:
                # Deleted bases after the VCF anchor.
                query_start = pos
                query_end = pos + len(ref) - 1

            sizes["genome"].append(indel_size)

            if overlaps_centromere(
                chrom,
                query_start,
                query_end,
                centromere_spans,
            ):
                sizes["centromere"].append(indel_size)

            if overlaps_index(
                chrom,
                query_start,
                query_end,
                repeat_index["knob180"],
            ):
                sizes["knob180"].append(indel_size)

            if overlaps_index(
                chrom,
                query_start,
                query_end,
                repeat_index["CentC"],
            ):
                sizes["CentC"].append(indel_size)

    return sizes, skipped


def read_syri_vcf(vcf_path):
    """
    Read sequence-resolved SyRI SVs >50 bp.

    Key: (normalized chromosome, POS, ID)

    Returns the size and the called sample names.
    """

    variants = {}

    with open_text(vcf_path) as infile:
        sample_names = []

        for line in infile:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]

            if (
                "," in alt
                or alt.startswith("<")
                or "[" in alt
                or "]" in alt
                or alt == "*"
            ):
                continue

            difference = len(alt) - len(ref)
            size = abs(difference)

            if size <= 50:
                continue

            called_samples = set()

            for sample, sample_field in zip(
                sample_names,
                fields[9:],
            ):
                gt = sample_field.split(":")[0]
                alleles = gt.replace("|", "/").split("/")

                if any(
                    allele not in ("0", ".")
                    for allele in alleles
                ):
                    called_samples.add(sample)

            key = (chrom, pos, variant_id)

            if key in variants:
                raise ValueError(
                    f"Duplicate SyRI VCF key: {key}"
                )

            variants[key] = {
                "size": size,
                "called_samples": called_samples,
            }

    return variants


def read_validated_svs(
    insertion_tsv,
    deletion_tsv,
    syri_variants,
    min_depth_ratio,
):
    """
    Match validation rows to the SyRI VCF.

    Validation:
        - depth ratio >= threshold
        - all three relevant spanning/junction checks PASS
        - variant is present in the SyRI VCF
        - validation sample carries the variant in the VCF

    The TSV's original verdict is not used because it may
    have been produced with a 0.75 depth-ratio threshold.
    """

    validated_sizes = []
    validated_keys = set()

    files = [
        (insertion_tsv, "INS"),
        (deletion_tsv, "DEL"),
    ]

    for path, sv_type in files:
        with open_text(path) as infile:
            reader = csv.DictReader(infile, delimiter="\t")

            for row in reader:
                try:
                    chrom = normalize_chrom(row["chrom"])
                    pos = int(row["pos"])
                    variant_id = row["id"]
                    sample = row["sample"]

                    key = (chrom, pos, variant_id)

                    if key not in syri_variants:
                        continue

                    variant = syri_variants[key]

                    if sample not in variant["called_samples"]:
                        continue

                    if sv_type == "INS":
                        depth_ratio = float(row["own_depth_ratio"])

                        checks = [
                            row["own_start_pass"],
                            row["own_end_pass"],
                            row["j702_junction_pass"],
                        ]

                    else:
                        depth_ratio = float(row["j702_depth_ratio"])

                        checks = [
                            row["own_junction_pass"],
                            row["j702_start_pass"],
                            row["j702_end_pass"],
                        ]

                    passed = (
                        depth_ratio >= min_depth_ratio
                        and all(
                            value.upper() == "PASS"
                            for value in checks
                        )
                    )

                    if not passed:
                        continue

                    event_key = (sample, *key)

                    if event_key in validated_keys:
                        raise ValueError(
                            f"Duplicate validated event: {event_key}"
                        )

                    validated_keys.add(event_key)
                    validated_sizes.append(variant["size"])

                except (KeyError, ValueError) as error:
                    # Missing or nonnumeric validation values
                    # cannot be classified as passing.
                    if isinstance(error, ValueError) and str(
                        error
                    ).startswith("Duplicate"):
                        raise

                    continue

    return validated_sizes


def sv_x_transform(values):
    """
    Display transform:

        50–100 bp: linear
        >100 bp: logarithmic

    The numeric axis labels remain the original bp values.
    """

    values = np.asarray(values, dtype=float)

    return np.where(
        values <= 100,
        values,
        100 + 50 * np.log10(values / 100),
    )


def make_sv_bins(max_size):
    """
    Use 10-bp bins from 50–100 bp.

    Above 100 bp use five logarithmic bins per tenfold
    increase in size.
    """

    linear_edges = np.arange(50, 101, 10)

    max_size = max(max_size, 100)

    # Extend the last edge beyond the largest event.
    max_exponent = max(
        2.2,
        math.ceil(math.log10(max_size) * 5) / 5,
    )

    if 10 ** max_exponent <= max_size:
        max_exponent += 0.2

    log_exponents = np.arange(
        2.2,
        max_exponent + 0.001,
        0.2,
    )

    log_edges = 10 ** log_exponents

    return np.concatenate([linear_edges, log_edges])


def plot_deepvariant_panel(ax, sizes, title):
    """
    Plot exact sizes 1–50 bp with a gap between bars.
    """

    sizes = np.asarray(sizes, dtype=int)

    bins = np.arange(0.5, 51.5, 1)

    ax.hist(
        sizes,
        bins=bins,
        color="#0072B2",
        rwidth=0.80,
    )

    ax.set_xlim(0.5, 50.5)

    ax.set_xticks([
        1, 5, 10, 15, 20, 25,
        30, 35, 40, 45, 50,
    ])

    # Linear through 10; logarithmic above 10.
    ax.set_yscale("symlog", linthresh=10)

    ax.set_ylabel(
        "Indel count\n(linear ≤10; log above)"
    )

    ax.set_title(
        f"{title} | n = {len(sizes):,}",
        loc="left",
        fontweight="bold",
    )

    ax.grid(axis="y", alpha=0.15)
    ax.set_axisbelow(True)


def plot_syri_panel(ax, sizes, title, bins, x_max):
    """
    Plot combined SyRI insertions and deletions.
    """

    sizes = np.asarray(sizes, dtype=float)

    counts, _ = np.histogram(sizes, bins=bins)

    positions = sv_x_transform(bins)
    widths = np.diff(positions)

    ax.bar(
        positions[:-1],
        counts,
        width=widths * 0.88,
        align="edge",
        color="#D55E00",
    )

    ax.set_yscale("symlog", linthresh=1)

    candidate_ticks = [
        50, 60, 70, 80, 90, 100,
        1000, 10000, 100000, 1000000,
    ]

    ticks = [
        tick for tick in candidate_ticks
        if tick <= x_max
    ]

    ax.set_xticks(sv_x_transform(ticks))
    ax.set_xticklabels(
        [f"{tick:,}" for tick in ticks],
        fontsize=8,
    )

    ax.set_xlim(
        sv_x_transform(50),
        sv_x_transform(x_max),
    )

    ax.set_ylabel("SV count\n(symlog scale)")

    ax.set_title(
        f"{title} | n = {len(sizes):,}",
        loc="left",
        fontweight="bold",
    )

    ax.grid(axis="y", alpha=0.15)
    ax.set_axisbelow(True)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot DeepVariant indels and validated SyRI SVs "
            "for centromeres, knob180, CentC and whole genome."
        )
    )

    parser.add_argument(
        "--deepvariant",
        required=True,
        help="Private DeepVariant indel VCF, plain or .gz",
    )

    parser.add_argument(
        "--syri",
        required=True,
        help="Filtered SyRI VCF, plain or .gz",
    )

    parser.add_argument(
        "--insertions",
        required=True,
        help="SyRI insertion validation TSV",
    )

    parser.add_argument(
        "--deletions",
        required=True,
        help="SyRI deletion validation TSV",
    )

    parser.add_argument(
        "--gff",
        required=True,
        help="J702 EDTA TE annotation GFF3",
    )

    parser.add_argument(
        "--depth-ratio",
        type=float,
        default=0.70,
        help="Minimum validation depth ratio; default 0.70",
    )

    parser.add_argument(
        "--output-prefix",
        default="ma_indel_sv_distributions",
        help="Prefix for PNG and PDF output",
    )

    args = parser.parse_args()

    repeat_index = read_repeat_intervals(args.gff)

    dv_sizes, skipped = read_deepvariant_indels(
        args.deepvariant,
        repeat_index,
    )

    syri_variants = read_syri_vcf(args.syri)

    genome_sv_sizes = read_validated_svs(
        args.insertions,
        args.deletions,
        syri_variants,
        args.depth_ratio,
    )

    if not genome_sv_sizes:
        raise SystemExit(
            "No validated SyRI SVs matched the VCF. "
            "Check the VCF, sample names and TSV files."
        )

    # The three region-associated SyRI lists are fixed to the
    # supplied tables. They were checked at depth ratio >=0.70.
    if args.depth_ratio != 0.70:
        print(
            "WARNING: The centromere, knob180 and CentC SV "
            "lists were supplied and checked at 0.70. "
            "Only the whole-genome panel is recalculated "
            "for the requested threshold."
        )

    all_sv_sizes = (
        genome_sv_sizes
        + CENTROMERE_SV_SIZES
        + KNOB180_SV_SIZES
        + CENTC_SV_SIZES
    )

    max_sv_size = max(all_sv_sizes)

    # Keep at least 100,000 bp on the shared SyRI x-axis,
    # like the earlier figure.
    x_max = max(100000, max_sv_size * 1.05)

    sv_bins = make_sv_bins(x_max)

    rows = [
        (
            "Tentative centromeres",
            dv_sizes["centromere"],
            CENTROMERE_SV_SIZES,
        ),
        (
            "knob180",
            dv_sizes["knob180"],
            KNOB180_SV_SIZES,
        ),
        (
            "CentC",
            dv_sizes["CentC"],
            CENTC_SV_SIZES,
        ),
        (
            "Whole genome",
            dv_sizes["genome"],
            genome_sv_sizes,
        ),
    ]

    fig, axes = plt.subplots(
        4,
        2,
        figsize=(15, 13),
    )

    for row_index, (title, dv, sv) in enumerate(rows):
        plot_deepvariant_panel(
            axes[row_index, 0],
            dv,
            title,
        )

        plot_syri_panel(
            axes[row_index, 1],
            sv,
            title,
            sv_bins,
            x_max,
        )

    axes[0, 0].text(
        0.5,
        1.23,
        "DeepVariant private indels: 1–50 bp",
        transform=axes[0, 0].transAxes,
        ha="center",
        fontweight="bold",
        fontsize=12,
    )

    axes[0, 1].text(
        0.5,
        1.23,
        "Validated SyRI insertions and deletions: >50 bp",
        transform=axes[0, 1].transAxes,
        ha="center",
        fontweight="bold",
        fontsize=12,
    )

    axes[-1, 0].set_xlabel(
        "Indel size (bp; exact sizes 1–50)"
    )

    axes[-1, 1].set_xlabel(
        "SV size (bp; linear 50–100, logarithmic above 100)"
    )

    fig.suptitle(
        "Insertion and deletion size distributions "
        "across five maize MA lines",
        fontsize=16,
        fontweight="bold",
    )

    fig.text(
        0.04,
        0.01,
        "DeepVariant repeat panels: overlap with J702 GFF. "
        "SyRI repeat panels: supplied event classifications.\n"
        "Centromeres: full outer spans including intervening gaps. "
        f"SyRI validation: depth ratio ≥{args.depth_ratio:.2f} "
        "and spanning checks.",
        fontsize=9,
    )

    fig.tight_layout(rect=[0, 0.05, 1, 0.96])

    png_path = f"{args.output_prefix}.png"
    pdf_path = f"{args.output_prefix}.pdf"

    fig.savefig(
        png_path,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    print("\nCounts used in the figure")
    print("Region\tDeepVariant_1-50bp\tSyRI_validated_gt50bp")

    for title, dv, sv in rows:
        print(f"{title}\t{len(dv)}\t{len(sv)}")

    print(f"\nSkipped non-biallelic/symbolic DV records: {skipped}")
    print(f"PNG: {png_path}")
    print(f"PDF: {pdf_path}")


if __name__ == "__main__":
    main()
