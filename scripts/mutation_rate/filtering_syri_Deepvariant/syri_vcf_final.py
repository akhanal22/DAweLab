#!/usr/bin/env python3
"""Export filtered SyRI calls of all sizes as a valid, GT-only VCF.

Run beside indel_filter_steps_by_column.py in the directory with vcf_to_compare.
Uses that module's SyRI input, 500 kb masks, reference repeat BEDs, and
insertion-boundary logic. Writes one row per retained sample and ALT allele.
"""

import gzip
import re
from collections import Counter

import indel_filter_steps_by_column as s


N_RUN = re.compile(r"N{3,}", re.IGNORECASE)
INFO_HEADER = re.compile(r"^##INFO=<ID=([^,>]+),Number=([^,>]+)")
OUTPUT = s.REPEAT_ROOT / "syri.filtered.all_variants.all_sizes.no_N3.vcf"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def selected_gt(gt, allele_number):
    """Map the chosen ALT to 1; turn other ALT alleles into missing calls."""
    return "".join(
        token if token in ("/", "|") else
        "1" if token == str(allele_number) else
        "0" if token == "0" else "."
        for token in re.split(r"([/|])", gt)
    )


def selected_info(info, allele_number, n_alts, numbers):
    """Retain safe annotations, updating allele-specific A and R fields."""
    if info == ".":
        return info
    kept = []
    for entry in info.split(";"):
        name, sep, value = entry.partition("=")
        if name in {"AC", "AN", "AF", "MLEAC", "MLEAF"}:
            continue  # Original allele counts no longer match masked samples.
        if n_alts > 1 and name == "END":
            continue  # END may describe another ALT in a mixed record.
        number = numbers.get(name)
        if number == "G":
            continue
        if sep and number == "A":
            vals = value.split(",")
            if len(vals) != n_alts:
                continue
            value = vals[allele_number - 1]
        elif sep and number == "R":
            vals = value.split(",")
            if len(vals) != n_alts + 1:
                continue
            value = vals[0] + "," + vals[allele_number]
        kept.append(name + ("=" + value if sep else ""))
    return ";".join(kept) or "."


def main():
    source = s.INPUTS["syri"]
    cg_path = s.REPEAT_ROOT / "reference.CG_homopolymers.min5.bed"
    paths = [
        source,
        cg_path,
        s.REPEAT_ROOT / "reference.CG_homopolymers.min5.plus1bp.bed",
        s.REPEAT_ROOT / "reference.perfect_dinucleotide_repeats.min3copies.bed",
        s.REPEAT_ROOT / "reference.perfect_trinucleotide_repeats.min3copies.bed",
    ]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")
    cg_core = s.read_bed(cg_path)
    cg_runs = s.read_core_runs(cg_path)
    cg_flank = s.read_bed(paths[2])
    di = s.read_bed(paths[3])
    tri = s.read_bed(paths[4])
    s.REPEAT_ROOT.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_name(OUTPUT.name + ".tmp")
    counts = Counter()
    numbers = {}
    try:
        with open_text(source) as src, temporary.open("w") as dst:
            names = None
            masks = None
            for line in src:
                if line.startswith("##"):
                    match = INFO_HEADER.match(line)
                    if match:
                        numbers[match.group(1)] = match.group(2)
                    dst.write(line)
                    continue
                if line.startswith("#CHROM"):
                    names = line.rstrip("\n").split("\t")[9:]
                    if not names:
                        raise ValueError("SyRI input has no sample columns")
                    masks = {}
                    for sample in names:
                        bed = (s.WINDOW_ROOT / "syri" / "bad_500kb_bins"
                               / f"{sample}.bad_500kb_bins.bed")
                        if not bed.is_file():
                            raise FileNotFoundError(f"Required sample mask missing: {bed}")
                        masks[sample] = s.read_bed(bed)
                    dst.write(line)
                    continue
                if not line.strip():
                    continue
                if names is None:
                    raise ValueError("Missing VCF #CHROM header")
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 + len(names):
                    raise ValueError("VCF record has an unexpected number of columns")
                chrom, pos, ref = fields[0], int(fields[1]), fields[3]
                alts = fields[4].split(",")
                fmt = fields[8].split(":")
                if "GT" not in fmt:
                    raise ValueError(f"Missing GT at {chrom}:{pos}")
                gt_index = fmt.index("GT")
                gts = []
                for sample_field in fields[9:]:
                    parts = sample_field.split(":")
                    gts.append(parts[gt_index] if gt_index < len(parts) else ".")

                for allele_number, alt in enumerate(alts, 1):
                    if alt in {".", "*", "<NON_REF>"}:
                        counts["skipped_nonvariant_allele"] += 1
                        continue
                    symbolic = alt.startswith("<") or "[" in alt or "]" in alt
                    has_n3 = bool(N_RUN.search(ref) or N_RUN.search(alt))
                    if symbolic:
                        # Symbolic ALT has no sequence length. Use INFO/END as
                        # its reference span where present; repeat-size rules
                        # cannot be evaluated for it.
                        match = re.search(r"(?:^|;)END=(\d+)(?:;|$)", fields[7])
                        start = pos - 1
                        end = max(pos, int(match.group(1))) if match else pos
                        inserted = False
                        size = None
                    else:
                        start, end, inserted = s.event_interval(pos, ref, alt)
                        size = abs(len(alt) - len(ref))

                    for idx, (sample, gt) in enumerate(zip(names, gts)):
                        if str(allele_number) not in re.split(r"[/|]", gt):
                            continue
                        counts["input_calls"] += 1
                        if has_n3:
                            counts["N3"] += 1
                            continue
                        if s.touches(masks[sample], chrom, pos - 1, pos):
                            counts["bad_500kb_bin"] += 1
                            continue
                        if s.in_cg_core(cg_core, cg_runs, chrom, start, end, inserted):
                            counts["CG_core"] += 1
                            continue
                        if s.touches(cg_flank, chrom, start, end):
                            counts["CG_plus1bp"] += 1
                            continue
                        if size is not None and 1 <= size <= s.MAX_REPEAT_INDEL_BP:
                            if s.touches(di, chrom, start, end):
                                counts["dinucleotide"] += 1
                                continue
                            if s.touches(tri, chrom, start, end):
                                counts["trinucleotide"] += 1
                                continue
                        sample_values = ["./."] * len(names)
                        sample_values[idx] = selected_gt(gt, allele_number)
                        info = selected_info(fields[7], allele_number, len(alts), numbers)
                        dst.write("\t".join([
                            *fields[:4], alt, fields[5], fields[6], info,
                            "GT", *sample_values,
                        ]) + "\n")
                        counts["retained_calls"] += 1
                        counts["symbolic" if symbolic else
                               "indel" if size else
                               "SNP" if len(ref) == len(alt) == 1 else
                               "other_sequence_variant"] += 1
        temporary.replace(OUTPUT)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Source: {source.resolve()}")
    print(f"Filtered VCF: {OUTPUT.resolve()}")
    for label in ("input_calls", "retained_calls", "SNP", "indel",
                  "other_sequence_variant", "symbolic", "bad_500kb_bin",
                  "CG_core", "CG_plus1bp", "dinucleotide",
                  "trinucleotide", "N3", "skipped_nonvariant_allele"):
        print(f"{label}: {counts[label]:,}")
    print("One VCF row per retained sample and ALT; other samples are ./.")


if __name__ == "__main__":
    main()
