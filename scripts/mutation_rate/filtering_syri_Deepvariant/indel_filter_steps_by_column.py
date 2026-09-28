#!/usr/bin/env python3
"""Count cumulative indel retention across the three caller VCFs.

Run from the directory containing vcf_to_compare. Uses reference-derived
repeat BEDs produced by 04-homopolymer_flanks_di_tri_repeat.sh.
"""

import bisect
import csv
import gzip
import re
from collections import Counter
from pathlib import Path

ROOT = Path("vcf_to_compare")
WINDOW_ROOT = ROOT / "three_callers_private_500kb"
REPEAT_ROOT = WINDOW_ROOT / "repeat_filter_comparison/04_hp_flank1_di_tri"
# Start all three caller inputs before their caller-specific 500-kb filter.
INPUTS = {
    "syri": WINDOW_ROOT / "syri/syri.private.snps_indels.vcf.gz",
    "anchorwave": WINDOW_ROOT / "anchorwave/anchorwave.private.snps_indels.vcf.gz",
    "deepvariant": WINDOW_ROOT / "deepvariant/deepvariant.private.snps_indels.vcf.gz",
}
OUTPUT = ROOT / "indel_filter_steps_by_caller.tsv"
WIDE_OUTPUT = ROOT / "indel_filter_steps_wide.tsv"
MAX_REPEAT_INDEL_BP = 50  # Match this analysis: 1-50 bp for di/tri.

STAGES = (
    "input_indels",
    "after_500kb",
    "after_CG_homopolymer_core",
    "after_CG_homopolymer_plus_1bp_flank",
    "after_dinucleotide_1_50bp",
    "after_trinucleotide_1_50bp",
)


def read_bed(path):
    if not path.is_file():
        raise SystemExit(f"Missing BED: {path}")
    intervals = {}
    with path.open() as source:
        for line in source:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, *_ = line.split()
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    result = {}
    for chrom, values in intervals.items():
        values.sort()
        starts = []
        prefix_max_end = []
        largest = -1
        for start, end in values:
            starts.append(start)
            largest = max(largest, end)
            prefix_max_end.append(largest)
        result[chrom] = (starts, prefix_max_end)
    return result


def read_core_runs(path):
    """Retain the C/G base so insertion at an edge follows the 04 script."""
    if not path.is_file():
        raise SystemExit(f"Missing BED: {path}")
    runs = {}
    with path.open() as source:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, base, *_ = line.split()
            runs.setdefault(chrom, []).append((int(start), int(end), base))
    for chrom in runs:
        runs[chrom].sort()
        runs[chrom] = ([run[0] for run in runs[chrom]], runs[chrom])
    return runs


def touches(intervals, chrom, start, end):
    """Test changed reference span; zero-length insertions touch run edges."""
    index = intervals.get(chrom)
    if index is None:
        return False
    starts, prefix_end = index
    if start == end:
        i = bisect.bisect_right(starts, start) - 1
        return i >= 0 and prefix_end[i] >= start
    i = bisect.bisect_left(starts, end) - 1
    return i >= 0 and prefix_end[i] > start


def event_interval(pos, ref, alt):
    """Trim shared prefix/suffix; return changed span in 0-based reference."""
    left = 0
    while left < min(len(ref), len(alt)) and ref[left] == alt[left]:
        left += 1
    r_end, a_end = len(ref), len(alt)
    while r_end > left and a_end > left and ref[r_end - 1] == alt[a_end - 1]:
        r_end -= 1
        a_end -= 1
    start = pos - 1 + left
    return start, start + r_end - left, alt[left:a_end]


def in_cg_core(bed, runs, chrom, start, end, inserted):
    if start != end:
        return touches(bed, chrom, start, end)
    entry = runs.get(chrom)
    if entry is None:
        return False
    starts, chrom_runs = entry
    i = bisect.bisect_right(starts, start) - 1
    for run_index in (i, i + 1):
        if 0 <= run_index < len(chrom_runs):
            core_start, core_end, base = chrom_runs[run_index]
            if core_start < start < core_end:
                return True
            if start in (core_start, core_end) and inserted and all(
                letter.upper() == base.upper() for letter in inserted
            ):
                return True
    return False


def called_alt_indices(gt, alt_count):
    return {
        int(value)
        for value in re.findall(r"\d+", gt)
        if 1 <= int(value) <= alt_count
    }


def genotype_group(gt):
    alleles = re.split(r"[|/]", gt)
    if len(alleles) == 2 and all(a.isdigit() for a in alleles):
        numbers = [int(a) for a in alleles]
        if numbers[0] != numbers[1]:
            return "HET"
        if numbers[0] > 0:
            return "HOM_ALT"
    return "OTHER_ALT"  # Haploid/partly missing SyRI genotypes stay visible.


def call_records(caller, path):
    if not path.is_file():
        raise SystemExit(f"Missing VCF: {path}")
    open_vcf = gzip.open if path.suffix == ".gz" else open
    with open_vcf(path, "rt") as source:
        sample_names = None
        for line in source:
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                yield ("names", sample_names)
                continue
            if line.startswith("#"):
                continue
            if sample_names is None:
                raise SystemExit(f"No sample header in {path}")
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref = fields[0], int(fields[1]), fields[3]
            alts = fields[4].split(",")
            fmt = fields[8].split(":")
            if "GT" not in fmt:
                raise SystemExit(f"Missing GT at {chrom}:{pos} in {path}")
            gt_index = fmt.index("GT")
            for sample, sample_value in zip(sample_names, fields[9:]):
                parts = sample_value.split(":")
                gt = parts[gt_index] if gt_index < len(parts) else "."
                indices = called_alt_indices(gt, len(alts))
                for alt_index in indices:
                    alt = alts[alt_index - 1]
                    if (not re.fullmatch(r"[ACGTNacgtn]+", ref)
                            or not re.fullmatch(r"[ACGTNacgtn]+", alt)
                            or len(ref) == len(alt)):
                        continue
                    yield ("indel", sample, chrom, pos, ref, alt,
                           genotype_group(gt))


def main():
    cg_core_path = REPEAT_ROOT / "reference.CG_homopolymers.min5.bed"
    cg_core = read_bed(cg_core_path)
    cg_core_runs = read_core_runs(cg_core_path)
    cg_flank = read_bed(REPEAT_ROOT / "reference.CG_homopolymers.min5.plus1bp.bed")
    di = read_bed(REPEAT_ROOT / "reference.perfect_dinucleotide_repeats.min3copies.bed")
    tri = read_bed(REPEAT_ROOT / "reference.perfect_trinucleotide_repeats.min3copies.bed")
    counts = Counter()
    names_by_caller = {}

    for caller, path in INPUTS.items():
        bins = {}
        for item in call_records(caller, path):
            if item[0] == "names":
                names_by_caller[caller] = item[1]
                bins = {
                    sample: read_bed(
                        WINDOW_ROOT / caller / "bad_500kb_bins"
                        / f"{sample}.bad_500kb_bins.bed"
                    )
                    for sample in item[1]
                }
                continue
            _, sample, chrom, pos, ref, alt, group = item
            key = (caller, sample, group)
            counts[key, STAGES[0]] += 1

            # Each BED was generated for the same caller and sample as input.
            bad_bin = touches(bins[sample], chrom, pos - 1, pos)
            if bad_bin:
                continue
            counts[key, STAGES[1]] += 1

            start, end, inserted = event_interval(pos, ref, alt)
            if in_cg_core(cg_core, cg_core_runs, chrom, start, end, inserted):
                continue
            counts[key, STAGES[2]] += 1

            if touches(cg_flank, chrom, start, end):
                continue
            counts[key, STAGES[3]] += 1

            size = abs(len(alt) - len(ref))
            if size <= MAX_REPEAT_INDEL_BP and touches(di, chrom, start, end):
                continue
            counts[key, STAGES[4]] += 1

            if size <= MAX_REPEAT_INDEL_BP and touches(tri, chrom, start, end):
                continue
            counts[key, STAGES[5]] += 1

    with OUTPUT.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(("caller", "sample", "genotype", *STAGES,
                         "removed_500kb", "removed_CG_core", "removed_CG_flank",
                         "removed_di_1_50bp", "removed_tri_1_50bp"))
        for caller, samples in names_by_caller.items():
            for sample in [*samples, "ALL_SAMPLES"]:
                for group in ("HET", "HOM_ALT", "OTHER_ALT", "ALL_GENOTYPES"):
                    if sample == "ALL_SAMPLES" and group == "ALL_GENOTYPES":
                        keys = [(caller, name, kind) for name in samples
                                for kind in ("HET", "HOM_ALT", "OTHER_ALT")]
                    elif sample == "ALL_SAMPLES":
                        keys = [(caller, name, group) for name in samples]
                    elif group == "ALL_GENOTYPES":
                        keys = [(caller, sample, kind) for kind in
                                ("HET", "HOM_ALT", "OTHER_ALT")]
                    else:
                        keys = [(caller, sample, group)]
                    values = [sum(counts[key, stage] for key in keys)
                              for stage in STAGES]
                    removed = [a - b for a, b in zip(values, values[1:])]
                    assert all(n >= 0 for n in removed)
                    writer.writerow((caller, sample, group, *values, *removed))
    print(f"Wrote {OUTPUT}")
    # One sheet-like table: step rows and the three callers as columns.
    with WIDE_OUTPUT.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(("genotype", "measure", "syri", "anchorwave", "deepvariant"))
        for group in ("HET", "HOM_ALT", "OTHER_ALT", "ALL_GENOTYPES"):
            values_by_caller = {}
            for caller, samples in names_by_caller.items():
                groups = ("HET", "HOM_ALT", "OTHER_ALT") if group == "ALL_GENOTYPES" else (group,)
                values_by_caller[caller] = [
                    sum(counts[(caller, sample, kind), stage]
                        for sample in samples for kind in groups)
                    for stage in STAGES
                ]
            for i, stage in enumerate(STAGES):
                writer.writerow((group, stage, *(values_by_caller[caller][i]
                                                 for caller in INPUTS)))
                if i:
                    writer.writerow((group, "removed_at_" + stage,
                                     *(values_by_caller[caller][i - 1] -
                                       values_by_caller[caller][i]
                                       for caller in INPUTS)))
    print(f"Wrote {WIDE_OUTPUT}")


if __name__ == "__main__":
    main()
