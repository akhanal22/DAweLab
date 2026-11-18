#!/usr/bin/env python3

import argparse
import bisect
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def normalize_chromosome(name):
    # chr1_RagTag_R / chr1_RagTag_ / chr1_RagTag -> chr1
    return re.sub(r"_RagTag.*$", "", name.strip())


def classify_annotation(classification):
    if classification.startswith(("DNA/", "MITE/")):
        return "DNA transposon"

    if classification.startswith(("LTR/", "LINE/", "SINE/")):
        return "RNA transposon"

    return "Other repeat"


# Full RNA element features.
# LTR and TSD subfeatures are excluded.
RNA_FEATURES = {
    "LTR_retrotransposon",
    "Gypsy_LTR_retrotransposon",
    "Copia_LTR_retrotransposon",
    "LINE_element",
    "L1_LINE_retrotransposon",
    "RTE_LINE_retrotransposon",
    "SINE_element",
    "retrotransposon",
}

CHILD_FEATURES = {
    "long_terminal_repeat",
    "target_site_duplication",
}


def read_insertions(path, convention):
    rows = []

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                raise ValueError(
                    f"Line {line_number}: expected chromosome start end"
                )

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                # Allow an optional column header.
                if not rows and fields[1].lower().startswith("start"):
                    continue
                raise ValueError(
                    f"Line {line_number}: coordinates must be integers"
                )

            if start > end:
                raise ValueError(
                    f"Line {line_number}: start is greater than end"
                )

            # SyRI:
            # start is the anchor before the insertion.
            # Inserted bases are start+1 through end.
            #
            # BED:
            # start is 0-based and end is exclusive.
            # Converted inserted bases are start+1 through end.
            #
            # Inclusive:
            # start and end are both inserted bases.
            inserted_start = (
                start if convention == "inclusive" else start + 1
            )

            if inserted_start < 1 or end < inserted_start:
                raise ValueError(
                    f"Line {line_number}: invalid insertion interval"
                )

            rows.append({
                "insertion_id": f"INS_{len(rows) + 1:03d}",
                "input_line": line_number,
                "chromosome": normalize_chromosome(fields[0]),
                "input_start": start,
                "input_end": end,
                "inserted_start_1based": inserted_start,
                "inserted_end_1based": end,
                "insertion_bp": end - inserted_start + 1,
                "matches": [],
                "best": {},
                "alternative_matches": [],
            })

    return rows


def scan_gff(
    path,
    source,
    insertions,
    min_coverage,
    boundary_tolerance,
    convention,
):
    """
    TEanno:
        Screen DNA transposons, MITEs and other repeats.

    intact:
        Screen full RNA transposon features only.
        Exclude outer repeat_region wrappers and LTR/TSD children.
    """

    by_chromosome = defaultdict(list)

    for insertion in insertions:
        by_chromosome[insertion["chromosome"]].append(insertion)

    # Index the insertion intervals for efficient GFF overlap searches.
    indices = {}

    for chromosome, chromosome_rows in by_chromosome.items():
        chromosome_rows.sort(
            key=lambda row: row["inserted_start_1based"]
        )

        starts = [
            row["inserted_start_1based"]
            for row in chromosome_rows
        ]

        prefix_max_ends = []
        maximum_end = 0

        for row in chromosome_rows:
            maximum_end = max(
                maximum_end,
                row["inserted_end_1based"],
            )
            prefix_max_ends.append(maximum_end)

        indices[chromosome] = (
            chromosome_rows,
            starts,
            prefix_max_ends,
        )

    seen = set()
    annotation_chromosomes = set()

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) != 9:
                raise ValueError(
                    f"{path}, line {line_number}: expected 9 GFF columns"
                )

            chromosome = normalize_chromosome(fields[0])
            annotation_chromosomes.add(chromosome)

            attributes = {}

            for item in fields[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = unquote(value)

            classification = attributes.get("classification", "")
            category = classify_annotation(classification)
            feature_type = fields[2]

            if feature_type in CHILD_FEATURES:
                continue

            if source == "intact":
                if category != "RNA transposon":
                    continue

                if feature_type not in RNA_FEATURES:
                    continue

            elif category == "RNA transposon":
                # RNA candidates must come from the intact GFF.
                continue

            if chromosome not in indices:
                continue

            annotation_start = int(fields[3])
            annotation_end = int(fields[4])

            annotation_key = (
                fields[0],
                annotation_start,
                annotation_end,
                feature_type,
                attributes.get("ID", ""),
                classification,
            )

            if annotation_key in seen:
                continue

            seen.add(annotation_key)

            chromosome_rows, starts, prefix_max_ends = indices[
                chromosome
            ]

            upper = bisect.bisect_right(starts, annotation_end)

            lower = bisect.bisect_left(
                prefix_max_ends,
                annotation_start,
                0,
                upper,
            )

            for insertion in chromosome_rows[lower:upper]:
                inserted_start = insertion["inserted_start_1based"]
                inserted_end = insertion["inserted_end_1based"]

                overlap_bp = max(
                    0,
                    min(inserted_end, annotation_end)
                    - max(inserted_start, annotation_start)
                    + 1,
                )

                if overlap_bp == 0:
                    continue

                coverage = overlap_bp / insertion["insertion_bp"]

                start_difference = (
                    annotation_start - inserted_start
                )

                end_difference = (
                    annotation_end - inserted_end
                )

                qualifies = (
                    coverage >= min_coverage
                    and abs(start_difference) <= boundary_tolerance
                    and abs(end_difference) <= boundary_tolerance
                )

                annotation = {
                    "category": category,
                    "classification": classification,
                    "annotation_source": source,
                    "annotation_chromosome": fields[0],
                    "annotation_start": annotation_start,
                    "annotation_end": annotation_end,
                    "annotation_bp": (
                        annotation_end - annotation_start + 1
                    ),
                    "feature_type": feature_type,
                    "annotation_id": attributes.get("ID", ""),
                    "annotation_name": attributes.get("Name", ""),
                    "parent": attributes.get("Parent", ""),
                    "strand": fields[6],
                    "method": attributes.get("method", ""),
                    "identity": attributes.get("identity", ""),
                    "ltr_identity": attributes.get("ltr_identity", ""),
                    "TSD": attributes.get(
                        "tsd", attributes.get("TSD", "")
                    ),
                    "motif": attributes.get("motif", ""),
                    "overlap_bp": overlap_bp,
                    "insertion_coverage": coverage,
                    "start_difference_bp": start_difference,
                    "end_difference_bp": end_difference,
                    "qualifies": qualifies,
                    "gff_line": line_number,
                }

                # Best annotation: greatest coverage, then closest boundaries.
                ranking = (
                    coverage,
                    -max(
                        abs(start_difference),
                        abs(end_difference),
                    ),
                    -abs(start_difference) - abs(end_difference),
                )

                previous = insertion["best"].get(category)

                if previous is None or ranking > previous[0]:
                    insertion["best"][category] = (
                        ranking,
                        annotation,
                    )

                if qualifies:
                    insertion["matches"].append(annotation)

                # Check whether the one-base coordinate interpretation
                # changes the result. This does not relax the tolerance.
                alternative_start = (
                    insertion["input_start"]
                    if convention != "inclusive"
                    else insertion["input_start"] + 1
                )

                alternative_length = (
                    inserted_end - alternative_start + 1
                )

                alternative_overlap = max(
                    0,
                    min(inserted_end, annotation_end)
                    - max(alternative_start, annotation_start)
                    + 1,
                )

                if alternative_length > 0:
                    alternative_pass = (
                        alternative_overlap / alternative_length
                        >= min_coverage
                        and abs(annotation_start - alternative_start)
                        <= boundary_tolerance
                        and abs(end_difference)
                        <= boundary_tolerance
                    )

                    if alternative_pass:
                        insertion["alternative_matches"].append(
                            annotation
                        )

    missing = (
        set(by_chromosome) - annotation_chromosomes
    )

    if missing:
        raise ValueError(
            f"{path}: chromosome names not found after normalization: "
            f"{sorted(missing)}"
        )


INSERTION_COLUMNS = [
    "insertion_id",
    "input_line",
    "chromosome",
    "input_start",
    "input_end",
    "inserted_start_1based",
    "inserted_end_1based",
    "insertion_bp",
]

ANNOTATION_COLUMNS = [
    "category",
    "classification",
    "annotation_source",
    "annotation_chromosome",
    "annotation_start",
    "annotation_end",
    "annotation_bp",
    "feature_type",
    "annotation_id",
    "annotation_name",
    "parent",
    "strand",
    "method",
    "identity",
    "ltr_identity",
    "TSD",
    "motif",
    "overlap_bp",
    "insertion_coverage",
    "start_difference_bp",
    "end_difference_bp",
    "qualifies",
    "gff_line",
]


def combine_record(insertion, annotation):
    record = {
        key: insertion[key]
        for key in INSERTION_COLUMNS
    }
    record.update(annotation)
    return record


def write_csv(path, columns, records):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(records)


def is_pif_harbinger(annotation):
    return (
        annotation["classification"] in {"MITE/DTH", "DNA/DTH"}
        or annotation["feature_type"]
        == "PIF_Harbinger_TIR_transposon"
    )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Match insertions to EDTA DNA/repeat annotations "
            "and intact RNA transposons."
        )
    )

    parser.add_argument("--insertions", required=True)
    parser.add_argument("--teanno", required=True)
    parser.add_argument("--intact", required=True)

    parser.add_argument(
        "--coordinates",
        choices=["syri", "inclusive", "bed"],
        default="syri",
    )

    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.90,
    )

    parser.add_argument(
        "--boundary-tolerance",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--min-len",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--max-len",
        type=int,
        default=None,
    )

    parser.add_argument("--outdir", required=True)

    args = parser.parse_args()

    if not 0 < args.min_coverage <= 1:
        parser.error("--min-coverage must be between 0 and 1")

    if args.boundary_tolerance < 0:
        parser.error("--boundary-tolerance cannot be negative")

    original = read_insertions(
        args.insertions,
        args.coordinates,
    )

    insertions = [
        row
        for row in original
        if row["insertion_bp"] >= args.min_len
        and (
            args.max_len is None
            or row["insertion_bp"] <= args.max_len
        )
    ]

    if not insertions:
        parser.error("No insertion intervals pass the size filter")

    scan_gff(
        args.teanno,
        "TEanno",
        insertions,
        args.min_coverage,
        args.boundary_tolerance,
        args.coordinates,
    )

    scan_gff(
        args.intact,
        "intact",
        insertions,
        args.min_coverage,
        args.boundary_tolerance,
        args.coordinates,
    )

    combined_matches = [
        combine_record(insertion, annotation)
        for insertion in insertions
        for annotation in insertion["matches"]
    ]

    all_records = []
    audit_141bp = []

    for insertion in insertions:
        categories = {
            annotation["category"]
            for annotation in insertion["matches"]
        }

        if insertion["matches"]:
            selected_annotation = insertion["matches"][0]
        elif insertion["best"]:
            selected_annotation = max(
                insertion["best"].values(),
                key=lambda item: item[0],
            )[1]
        else:
            selected_annotation = {}

        coordinate_sensitive = (
            bool(insertion["matches"])
            != bool(insertion["alternative_matches"])
        )

        record = combine_record(
            insertion,
            selected_annotation,
        )

        record.update({
            "result": (
                "PASS" if insertion["matches"] else "FAIL"
            ),
            "dna_match": "DNA transposon" in categories,
            "intact_rna_match": "RNA transposon" in categories,
            "other_repeat_match": "Other repeat" in categories,
            "qualifying_annotations": len(insertion["matches"]),
            "alternative_coordinate_pass": bool(
                insertion["alternative_matches"]
            ),
            "coordinate_sensitive": coordinate_sensitive,
        })

        all_records.append(record)

        if insertion["insertion_bp"] == 141:
            qualifying_pif = [
                annotation
                for annotation in insertion["matches"]
                if is_pif_harbinger(annotation)
            ]

            if qualifying_pif:
                audit_annotation = qualifying_pif[0]
            else:
                audit_annotation = insertion["best"].get(
                    "DNA transposon",
                    (None, {}),
                )[1]

            audit_record = combine_record(
                insertion,
                audit_annotation,
            )

            audit_record.update({
                "pif_harbinger_pass": bool(qualifying_pif),
                "coordinate_sensitive": coordinate_sensitive,
            })

            audit_141bp.append(audit_record)

    summary = {
        "input_intervals": len(original),
        "tested_intervals": len(insertions),
        "coordinate_convention": args.coordinates,
        "min_coverage": args.min_coverage,
        "boundary_tolerance_bp": args.boundary_tolerance,
        "min_length": args.min_len,
        "max_length": args.max_len,
        "qualifying_insertions": sum(
            bool(row["matches"])
            for row in insertions
        ),
        "qualifying_annotation_rows": len(combined_matches),
        "unique_insertions_by_category": {
            category: sum(
                any(
                    annotation["category"] == category
                    for annotation in row["matches"]
                )
                for row in insertions
            )
            for category in [
                "DNA transposon",
                "RNA transposon",
                "Other repeat",
            ]
        },
        "dna_families": dict(Counter(
            annotation["classification"]
            for row in insertions
            for annotation in row["matches"]
            if annotation["category"] == "DNA transposon"
        )),
        "rna_families": dict(Counter(
            annotation["classification"]
            for row in insertions
            for annotation in row["matches"]
            if annotation["category"] == "RNA transposon"
        )),
        "insertions_141bp": len(audit_141bp),
        "pif_harbinger_matches_141bp": sum(
            row["pif_harbinger_pass"]
            for row in audit_141bp
        ),
        "coordinate_sensitive_intervals": sum(
            row["coordinate_sensitive"]
            for row in all_records
        ),
        "sources": {
            "insertions": args.insertions,
            "TEanno": args.teanno,
            "intact": args.intact,
        },
    }

    output_directory = Path(args.outdir)
    output_directory.mkdir(parents=True, exist_ok=True)

    match_columns = (
        INSERTION_COLUMNS + ANNOTATION_COLUMNS
    )

    write_csv(
        output_directory / "combined_matches.csv",
        match_columns,
        combined_matches,
    )

    for category, filename in [
        ("DNA transposon", "dna_matches.csv"),
        ("RNA transposon", "rna_intact_matches.csv"),
        ("Other repeat", "other_repeat_matches.csv"),
    ]:
        write_csv(
            output_directory / filename,
            match_columns,
            [
                row
                for row in combined_matches
                if row["category"] == category
            ],
        )

    write_csv(
        output_directory / "all_insertions.csv",
        INSERTION_COLUMNS + [
            "result",
            "dna_match",
            "intact_rna_match",
            "other_repeat_match",
            "qualifying_annotations",
            "alternative_coordinate_pass",
            "coordinate_sensitive",
        ] + ANNOTATION_COLUMNS,
        all_records,
    )

    write_csv(
        output_directory / "audit_141bp.csv",
        INSERTION_COLUMNS + [
            "pif_harbinger_pass",
            "coordinate_sensitive",
        ] + ANNOTATION_COLUMNS,
        audit_141bp,
    )

    with open(
        output_directory / "summary.json",
        "w",
    ) as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
