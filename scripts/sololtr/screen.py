#!/usr/bin/env python3

import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path("/workspace/scratch/48118ffd80f5")
GFF = ROOT / "upload/J702_flip.fa.mod.EDTA.intact.gff3"
DELETIONS = (
    ROOT / "upload/deletion_validation_withcigar_for500kbfilteredvcfagainstj702.tsv"
)
OUTDIR = ROOT / "solo_LTR_screen"
MAX_DISTANCE_BP = 5

OUTDIR.mkdir(parents=True, exist_ok=True)

# Group LTR annotations by parent element. Convert GFF3 coordinates
# from 1-based inclusive to 0-based half-open intervals.
parents = {}
ltrs_by_parent = defaultdict(list)

with GFF.open() as handle:
    for line in handle:
        if line.startswith("##FASTA"):
            break
        if line.startswith("#"):
            continue

        columns = line.rstrip("\n").split("\t")
        if len(columns) != 9:
            continue

        chrom, _, feature, start, end, _, _, _, attributes = columns
        if feature not in {"repeat_region", "long_terminal_repeat"}:
            continue

        attrs = dict(
            item.split("=", 1) for item in attributes.split(";") if "=" in item
        )

        if feature == "repeat_region":
            parents[attrs["ID"]] = attrs
        else:
            ltrs_by_parent[attrs["Parent"]].append(
                (chrom, int(start) - 1, int(end))
            )

elements_by_chrom = defaultdict(list)

for parent_id, ltrs in ltrs_by_parent.items():
    if len(ltrs) != 2:
        raise ValueError(f"{parent_id} has {len(ltrs)} LTRs. Expected two.")

    left, right = sorted(ltrs, key=lambda ltr: ltr[1])
    chrom, left_start, left_end = left
    right_chrom, right_start, right_end = right

    if chrom != right_chrom or left_end > right_start:
        raise ValueError(f"Invalid LTR pair for {parent_id}: {left}, {right}")

    elements_by_chrom[chrom].append(
        (parent_id, left_start, left_end, right_start, right_end)
    )

required_pass = ("own_junction_pass", "j702_start_pass", "j702_end_pass")

fields = [
    "sample", "chrom", "vcf_pos_1based", "deleted_start_1based",
    "deleted_end_1based", "deletion_length", "id", "element_id",
    "classification", "left_LTR_start_1based", "left_LTR_end_1based",
    "right_LTR_start_1based", "right_LTR_end_1based",
    "breakpoint_distance_left_bp", "breakpoint_distance_right_bp",
    "tier", "LTR_identity", "remaining_element_bp",
    "left_LTR_length", "right_LTR_length", "depth_ratio",
]

validated_count = 0
overlapping_deletions = set()
candidate_deletions = set()
candidate_pairs = 0
out_tsv = OUTDIR / "solo_LTR_candidates_within5bp.tsv"

with DELETIONS.open(newline="") as source, out_tsv.open("w", newline="") as output:
    reader = csv.DictReader(source, delimiter="\t")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for row in reader:
        length = int(row["length"])
        depth_ratio = float(row["j702_depth_ratio"])

        if length <= 50 or not (depth_ratio >= 0.70):
            continue
        if any(row[column] != "PASS" for column in required_pass):
            continue

        validated_count += 1
        chrom = row["chrom"]
        pos = int(row["pos"])
        end = int(row["end"])

        # For these simple VCF deletions, POS is the anchor base.
        # The deleted bases are POS+1 through END, inclusive.
        # The corresponding 0-based cut coordinates are POS and END.
        if end - pos != length:
            raise ValueError(
                f"{row['sample']} {row['id']}: END - POS does not match "
                f"the reported deletion length ({length})."
            )

        deletion_key = (row["sample"], chrom, pos, end)

        for parent_id, a, b, c, d in elements_by_chrom[chrom]:
            if not (pos < d and end > a):
                continue

            overlapping_deletions.add(deletion_key)

            # Distance to each LTR interval, measured between cut coordinates.
            # A cut inside the LTR or on either edge has distance zero.
            left_distance = max(a - pos, 0, pos - b)
            right_distance = max(c - end, 0, end - d)

            if left_distance > MAX_DISTANCE_BP or right_distance > MAX_DISTANCE_BP:
                continue

            attrs = parents[parent_id]
            candidate_deletions.add(deletion_key)
            candidate_pairs += 1

            writer.writerow({
                "sample": row["sample"],
                "chrom": chrom,
                "vcf_pos_1based": pos,
                "deleted_start_1based": pos + 1,
                "deleted_end_1based": end,
                "deletion_length": length,
                "id": row["id"],
                "element_id": parent_id,
                "classification": attrs.get("classification", ""),
                "left_LTR_start_1based": a + 1,
                "left_LTR_end_1based": b,
                "right_LTR_start_1based": c + 1,
                "right_LTR_end_1based": d,
                "breakpoint_distance_left_bp": left_distance,
                "breakpoint_distance_right_bp": right_distance,
                "tier": "paired_LTR_within_5bp",
                "LTR_identity": attrs.get("ltr_identity", ""),
                "remaining_element_bp": (d - a) - (min(end, d) - max(pos, a)),
                "left_LTR_length": b - a,
                "right_LTR_length": d - c,
                "depth_ratio": depth_ratio,
            })

summary = {
    "breakpoint_tolerance_bp": MAX_DISTANCE_BP,
    "intact_LTR_elements": len(ltrs_by_parent),
    "validated_deletions": validated_count,
    "deletions_overlapping_intact_LTR_elements": len(overlapping_deletions),
    "candidate_deletion_element_pairs": candidate_pairs,
    "within5bp_unique_deletions": len(candidate_deletions),
}

summary_text = json.dumps(summary, indent=2)
(OUTDIR / "summary.json").write_text(summary_text + "\n")
print(summary_text)
print(f"Candidate table: {out_tsv}")
