#!/usr/bin/env python3
"""Make three indel plots and export retained SyRI/DeepVariant VCF calls."""

import gzip
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, FixedLocator

import indel_filter_steps_by_column as s


EXTRA_DV_VCF = (
    s.ROOT
    / "filtered_0.95hom_40-60her_dv"
    / "private_biallelic_indels.all_sizes.exactly_one_line.chr1-10.vcf.gz"
)
MAX_PLOT_BP = 50_000
N_RUN = re.compile(r"N{10,}", re.IGNORECASE)


def has_n_run(ref, alt):
    return bool(N_RUN.search(ref) or N_RUN.search(alt))


def count_filtered_indels():
    counts = {caller: Counter() for caller in s.INPUTS}
    syri_retained_50plus = defaultdict(set)
    syri_audit = Counter()

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

    for caller, vcf in s.INPUTS.items():
        bad_bins = {}

        for item in s.call_records(caller, vcf):
            if item[0] == "names":
                bad_bins = {
                    sample: s.read_bed(
                        s.WINDOW_ROOT
                        / caller
                        / "bad_500kb_bins"
                        / f"{sample}.bad_500kb_bins.bed"
                    )
                    for sample in item[1]
                }
                continue

            _, sample, chrom, pos, ref, alt, group = item
            if s.touches(bad_bins[sample], chrom, pos - 1, pos):
                continue

            start, end, inserted = s.event_interval(pos, ref, alt)

            if s.in_cg_core(
                cg_core, cg_runs, chrom, start, end, inserted
            ):
                continue
            if s.touches(cg_flank, chrom, start, end):
                continue

            size = abs(len(alt) - len(ref))

            # These repeat filters apply only to indels of 1–50 bp.
            if size <= s.MAX_REPEAT_INDEL_BP and s.touches(
                di, chrom, start, end
            ):
                continue
            if size <= s.MAX_REPEAT_INDEL_BP and s.touches(
                tri, chrom, start, end
            ):
                continue

            if caller == "syri" and size >= 50:
                syri_audit["before_N10_50plus"] += 1
            if has_n_run(ref, alt):
                if caller == "syri" and size >= 50:
                    syri_audit["excluded_N10_50plus"] += 1
                continue

            counts[caller][size] += 1
            if caller == "syri" and size >= 50:
                syri_retained_50plus[chrom, pos, ref, alt].add(sample)
                syri_audit["retained_50_to_50k" if size <= MAX_PLOT_BP
                           else "retained_above_50k"] += 1

    return counts, syri_retained_50plus, syri_audit


def count_extra_dv_indels():
    """Count indels in the additional VCF without applying the above filters."""
    counts = Counter()
    dv_retained_1_to_50k = defaultdict(set)

    for item in s.call_records("deepvariant", EXTRA_DV_VCF):
        if item[0] == "names":
            continue

        _, sample, chrom, pos, ref, alt, group = item
        if has_n_run(ref, alt):
            continue
        size = abs(len(alt) - len(ref))
        counts[size] += 1
        if 1 <= size <= MAX_PLOT_BP:
            dv_retained_1_to_50k[chrom, pos, ref, alt].add(sample)

    return counts, dv_retained_1_to_50k


def project_gt(gt, alt_index):
    """Keep selected ALT as 1, REF as 0, and other ALTs as missing."""
    return "".join(
        token if token in ("/", "|") else
        "1" if token == str(alt_index) else
        "0" if token == "0" else "."
        for token in re.split(r"([/|])", gt)
    )

def subset_info(info, alt_index, alt_count, number_by_id):
    """Update allele-specific INFO values after making a biallelic row."""
    if info == ".":
        return info
    result = []
    for item in info.split(";"):
        name, separator, value = item.partition("=")
        if name in {"AC", "AN", "AF", "MLEAC", "MLEAF"}:
            continue  # Original allele counts no longer describe masked samples.
        if alt_count > 1 and name == "END":
            continue  # Could describe a different ALT in a multiallelic row.
        number = number_by_id.get(name)
        if separator and number == "G":
            continue
        if separator and number == "A":
            values = value.split(",")
            if len(values) != alt_count:
                continue
            value = values[alt_index - 1]
        if separator and number == "R":
            values = value.split(",")
            if len(values) != alt_count + 1:
                continue
            value = f"{values[0]},{values[alt_index]}"
        result.append(name + ("=" + value if separator else ""))
    return ";".join(result) or "."

def write_selected_vcf(source_path, output_path, chosen):
    """Preserve source site metadata and sample columns; output GT-only calls."""
    opener = gzip.open if source_path.suffix == ".gz" else open
    written = 0
    number_by_id = {}
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    with opener(source_path, "rt") as source, temporary_path.open("w") as output:
        sample_names = None
        for line in source:
            if line.startswith("##INFO=<"):
                match = re.match(r"##INFO=<ID=([^,]+),Number=([^,]+)", line)
                if match:
                    number_by_id[match.group(1)] = match.group(2)
            if line.startswith("##"):
                output.write(line)
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                output.write(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref = fields[0], int(fields[1]), fields[3]
            alts = fields[4].split(",")
            fmt = fields[8].split(":")
            gt_index = fmt.index("GT")
            for alt_index, alt in enumerate(alts, 1):
                selected = chosen.get((chrom, pos, ref, alt), set())
                for sample in sample_names:
                    if sample not in selected:
                        continue
                    sample_index = sample_names.index(sample)
                    parts = fields[9 + sample_index].split(":")
                    original_gt = parts[gt_index]
                    gt = project_gt(original_gt, alt_index)
                    sample_values = ["./."] * len(sample_names)
                    sample_values[sample_index] = gt
                    info = subset_info(fields[7], alt_index, len(alts), number_by_id)
                    output.write("\t".join([
                        *fields[:4], alt, *fields[5:7], info, "GT",
                        *sample_values,
                    ]) + "\n")
                    written += 1
    expected = sum(len(samples) for samples in chosen.values())
    if written != expected:
        temporary_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"{output_path}: wrote {written} calls but selected {expected}"
        )
    temporary_path.replace(output_path)
    print(f"Saved VCF: {output_path} ({written:,} calls, "
          f"{output_path.stat().st_size:,} bytes)", flush=True)


def style_axes(ax, largest_size, first_size=1):
    """Match the reference figure's exact-size bars and axis scales."""
    if first_size == 1:
        ax.set_xscale("symlog", linthresh=10, base=10)
        ax.set_xlim(0.5, max(largest_size * 1.05, 1.5))
        ticks = [1, 2, 5, 10]
    else:
        ax.set_xscale("log", base=10)
        ax.set_xlim(first_size, max(largest_size * 1.05, first_size + 1))
        ticks = [first_size]
    for decade in (10 ** power for power in range(1, 12)):
        ticks.extend(mult * decade for mult in (2, 5, 10))
    ax.xaxis.set_major_locator(
        FixedLocator(sorted({tick for tick in ticks if first_size <= tick <= largest_size}))
    )
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    ax.tick_params(axis="x", labelrotation=40)
    ax.set_yscale("symlog", linthresh=10, base=10)
    ax.set_xlabel(
        "Absolute indel size (bp; linear through 10, log above)"
        if first_size == 1 else "Absolute indel size (bp; log scale)"
    )
    ax.set_ylabel("Indel calls per exact size (linear through 10, log above)")
    ax.grid(axis="y", which="both", alpha=0.18)
    ax.set_axisbelow(True)


def main():
    s.REPEAT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"Output folder: {s.REPEAT_ROOT.resolve()}", flush=True)
    print("Filtering three caller VCFs...", flush=True)
    counts, syri_selected, syri_audit = count_filtered_indels()
    print("SyRI calls >=50 bp after the original location/repeat filters:",
          f"{syri_audit['before_N10_50plus']:,}", flush=True)
    print("  Excluded because REF/ALT contains >=10 consecutive Ns:",
          f"{syri_audit['excluded_N10_50plus']:,}", flush=True)
    print("  Retained from 50 through 50,000 bp (plotted):",
          f"{syri_audit['retained_50_to_50k']:,}", flush=True)
    print("  Retained above 50,000 bp (full SyRI VCF only):",
          f"{syri_audit['retained_above_50k']:,}", flush=True)
    assert (syri_audit["before_N10_50plus"] ==
            syri_audit["excluded_N10_50plus"] +
            syri_audit["retained_50_to_50k"] +
            syri_audit["retained_above_50k"])
    summary = s.REPEAT_ROOT / "syri_50bp_count_breakdown.tsv"
    with summary.open("w") as output:
        output.write("measure\tcalls\n")
        for measure in ("before_N10_50plus", "excluded_N10_50plus",
                        "retained_50_to_50k", "retained_above_50k"):
            output.write(f"{measure}\t{syri_audit[measure]}\n")
    print(f"Saved count breakdown: {summary}", flush=True)
    print("Reading the supplied DeepVariant VCF...", flush=True)
    counts["extra_dv"], dv_selected = count_extra_dv_indels()

    # Save both the complete SyRI >=50 bp set and the plotted 50–50,000 bp
    # subset. GT-only output avoids stale FORMAT fields when a multiallelic
    # record is separated into one row per call.
    write_selected_vcf(
        s.INPUTS["syri"],
        s.REPEAT_ROOT / "syri.filtered.no_N10.50bp_and_larger.vcf",
        syri_selected,
    )
    syri_selected_to_50k = {
        key: samples for key, samples in syri_selected.items()
        if abs(len(key[3]) - len(key[2])) <= MAX_PLOT_BP
    }
    write_selected_vcf(
        s.INPUTS["syri"],
        s.REPEAT_ROOT / "syri.filtered.no_N10.50-50000bp.vcf",
        syri_selected_to_50k,
    )
    write_selected_vcf(
        EXTRA_DV_VCF,
        s.REPEAT_ROOT / "deepvariant.supplied.no_N10.1-50000bp.vcf",
        dv_selected,
    )

    series = [
        ("syri", "SyRI filtered", "#D95F02"),
        ("anchorwave", "AnchorWave filtered", "#009E73"),
        ("deepvariant", "DeepVariant filtered", "#0072B2"),
        ("extra_dv", "DeepVariant supplied VCF", "#CC79A7"),
    ]

    print("Indel calls (all sizes, excluding alleles with >=10 Ns):")
    for key, label, _ in series:
        print(f"  {label}: {sum(counts[key].values()):,}")

    # Every bar represents one exact size; this figure displays 1–50,000 bp.
    fig, ax = plt.subplots(figsize=(19, 7))
    offsets = (-0.36, -0.12, 0.12, 0.36)
    for (key, label, color), offset in zip(series, offsets):
        xs = sorted(size for size, n in counts[key].items()
                    if 1 <= size <= MAX_PLOT_BP and n)
        # Keep grouped bars visible on the logarithmic part of the x-axis.
        positions = [
            size + offset if size <= 10 else size * (1 + offset * 0.12)
            for size in xs
        ]
        widths = [0.23 if size <= 10 else 0.025 * size for size in xs]
        ax.bar(
            positions,
            [counts[key][size] for size in xs],
            width=widths, color=color,
            label=f"{label} ({sum(counts[key][size] for size in xs):,})",
        )
    style_axes(ax, MAX_PLOT_BP)
    ax.set_xlim(0.5, MAX_PLOT_BP)
    ax.set_title("Private indels, 1–50,000 bp: four VCF series (no N10 alleles)",
                 fontweight="bold")
    ax.legend(frameon=False)
    fig.tight_layout()
    output = s.REPEAT_ROOT / "filtered_callers_plus_supplied_DV_all_sizes.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved: {output}")

    # The supplied DeepVariant VCF comes from the hom-0.95/het-0.40–0.60
    # folder. Count its records as supplied, without repeating other filters.
    fig, ax = plt.subplots(figsize=(19, 7))
    sizes = range(1, 51)
    ax.bar(sizes, [counts["extra_dv"][size] for size in sizes],
           width=0.8, color="#CC79A7")
    style_axes(ax, max((size for size in sizes if counts["extra_dv"][size]),
                       default=50))
    ax.set_title(
        f"Supplied DeepVariant VCF: 1–50 bp "
        f"({sum(counts['extra_dv'][size] for size in sizes):,} calls)",
        fontweight="bold",
    )
    fig.tight_layout()
    output = s.REPEAT_ROOT / "supplied_DV_hom0.95_1-50bp.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved: {output}")

    fig, ax = plt.subplots(figsize=(19, 7))
    syri_sizes = sorted(size for size, n in counts["syri"].items()
                        if 50 <= size <= MAX_PLOT_BP and n)
    ax.bar(syri_sizes, [counts["syri"][size] for size in syri_sizes],
           width=[0.035 * size for size in syri_sizes], color="#D95F02")
    style_axes(ax, MAX_PLOT_BP, first_size=50)
    ax.set_xlim(50, MAX_PLOT_BP)
    ax.set_title(
        f"Filtered SyRI indels: 50–50,000 bp "
        f"({sum(counts['syri'][size] for size in syri_sizes):,} calls)",
        fontweight="bold",
    )
    fig.tight_layout()
    output = s.REPEAT_ROOT / "filtered_syri_50bp_and_larger.png"
    fig.savefig(output, dpi=200)
    plt.close(fig)
    print(f"Saved: {output}")
    over_limit = sum(n for size, n in counts["syri"].items()
                     if size > MAX_PLOT_BP)
    print(f"SyRI VCF also retains {over_limit:,} calls above 50,000 bp")


if __name__ == "__main__":
    main()
