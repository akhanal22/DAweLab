#!/usr/bin/env python3
"""Compare CDS CG methylation using identical sites in each sample pair.

Input BED columns (pb-CpG-tools):
chrom, begin, end, mod_score, type, cov,
est_mod_count, est_unmod_count

Only rows with type == Total are used.
BED coordinates must refer to the same J702 reference in every sample.
"""

import argparse
import bisect
import csv
import gzip
import math
from collections import defaultdict, Counter
from pathlib import Path
from urllib.parse import unquote


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def write_table(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_cds_index(path):
    """Resolve CDS parents to genes and merge CDS intervals per gene."""
    parents = {}
    genes = {}
    cds = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attributes = fields
            attrs = {}
            for item in attributes.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = value

            feature_id = unquote(attrs.get("ID", ""))
            parent_ids = [
                unquote(p)
                for p in attrs.get("Parent", "").split(",")
                if p
            ]

            if feature_id:
                parents[feature_id] = parent_ids

            if feature.lower() == "gene" and feature_id:
                genes[feature_id] = (chrom, strand)

            if feature == "CDS":
                # Convert GFF 1-based inclusive to BED coordinates.
                cds.append(
                    (chrom, int(start) - 1, int(end), parent_ids)
                )

    cache = {}

    def gene_ancestors(node, visiting=None):
        if node in cache:
            return cache[node]
        if node in genes:
            return {node}

        visiting = set() if visiting is None else visiting
        if node in visiting:
            raise ValueError(f"Cycle in GFF Parent relationships: {node}")

        result = set()
        for parent in parents.get(node, []):
            result.update(gene_ancestors(parent, visiting | {node}))

        cache[node] = result
        return result

    intervals = defaultdict(list)

    for chrom, start, end, parent_ids in cds:
        resolved = set()
        for parent in parent_ids:
            resolved.update(gene_ancestors(parent))

        for gene in resolved:
            intervals[(chrom, gene)].append((start, end))

    if not intervals:
        raise ValueError("No CDS intervals could be linked to gene IDs.")

    by_chrom = defaultdict(list)
    cds_genes = set()

    for (chrom, gene), spans in intervals.items():
        merged = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))

        for start, end in merged:
            by_chrom[chrom].append((start, end, gene))
        cds_genes.add(gene)

    index = {}
    for chrom, spans in by_chrom.items():
        spans.sort()
        starts = [span[0] for span in spans]
        prefix_max_end = []
        maximum = -1
        for _, end, _ in spans:
            maximum = max(maximum, end)
            prefix_max_end.append(maximum)
        index[chrom] = (spans, starts, prefix_max_end)

    return index, cds_genes


def overlapping_genes(index, chrom, start, end):
    if chrom not in index:
        return set()

    spans, starts, prefix_max_end = index[chrom]
    i = bisect.bisect_left(starts, end) - 1
    result = set()

    while i >= 0 and prefix_max_end[i] > start:
        left, right, gene = spans[i]
        if left < end and right > start:
            result.add(gene)
        i -= 1

    return result


def bed_rows(path, min_cov, max_cov):
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if (
                not line.strip()
                or line.startswith(("#", "track", "browser"))
            ):
                continue

            fields = line.split()
            if len(fields) < 8:
                raise ValueError(
                    f"{path}:{line_number}: expected at least 8 columns"
                )

            if fields[4] != "Total":
                continue

            try:
                chrom = fields[0]
                start, end = int(fields[1]), int(fields[2])
                coverage = float(fields[5])
                modified = float(fields[6])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid numeric field"
                ) from exc

            if not (
                math.isfinite(coverage)
                and math.isfinite(modified)
                and end > start
                and 0 <= modified <= coverage
            ):
                raise ValueError(
                    f"{path}:{line_number}: invalid counts or interval"
                )

            if min_cov <= coverage <= max_cov:
                yield (chrom, start, end), modified, coverage


def status(methylation, um_max, methylated_min):
    if methylation <= um_max:
        return "UM"
    if methylation >= methylated_min:
        return "CG_methylated"
    return "ambiguous"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--founder", required=True)
    parser.add_argument(
        "--sample", action="append", required=True,
        help="Repeat as --sample NAME=path/to/sample.bed.gz"
    )
    parser.add_argument("--min-cov", type=float, default=40)
    parser.add_argument("--max-cov", type=float, default=100)
    parser.add_argument("--min-sites", type=int, default=20)
    parser.add_argument("--um-max", type=float, required=True)
    parser.add_argument("--methylated-min", type=float, required=True)
    parser.add_argument("--out", default="shared_CDS_CG_vs_J702")
    args = parser.parse_args()

    if not 0 <= args.um_max < args.methylated_min <= 1:
        parser.error("Require 0 <= um-max < methylated-min <= 1")
    if not 0 < args.min_cov <= args.max_cov or args.min_sites < 1:
        parser.error("Invalid coverage or minimum-site settings")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    index, cds_genes = read_cds_index(args.gff)
    print(f"Genes with annotated CDS: {len(cds_genes):,}", flush=True)

    # Store only qualifying founder sites overlapping CDS.
    founder = {}
    for key, modified, coverage in bed_rows(
        args.founder, args.min_cov, args.max_cov
    ):
        gene_ids = overlapping_genes(index, *key)
        if not gene_ids:
            continue
        if key in founder:
            raise ValueError(f"Duplicate founder Total site: {key}")
        founder[key] = (modified, coverage, gene_ids)

    if not founder:
        raise ValueError(
            "No qualifying founder sites overlap CDS. "
            "Check reference contig names and coverage filters."
        )

    print(f"Qualifying J702 CDS CG sites: {len(founder):,}", flush=True)

    detail_fields = [
        "gene_id", "shared_CG_sites",
        "J702_modified_sum", "J702_coverage_sum",
        "sample_modified_sum", "sample_coverage_sum",
        "J702_CG_fraction", "sample_CG_fraction",
        "difference_sample_minus_J702",
        "J702_status", "sample_status", "change",
    ]

    summary_fields = [
        "sample", "CDS_genes",
        "genes_with_20_or_more_shared_sites",
        "confident_in_both", "gain", "loss",
        "retained_methylated", "retained_UM",
        "ambiguous_in_either", "insufficient_shared_sites",
    ]
    # Keep the heading accurate if --min-sites is changed.
    summary_fields[2] = "genes_passing_min_shared_sites"
    summaries = []

    for specification in args.sample:
        if "=" not in specification:
            parser.error("--sample must be NAME=path")

        name, path = specification.split("=", 1)
        if not name or "/" in name or "\\" in name:
            parser.error("Sample name must be a nonempty filename-safe name")

        print(f"Comparing {name} with J702...", flush=True)

        # n_sites, founder modified, founder coverage,
        # sample modified, sample coverage
        totals = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0])
        seen = set()

        for key, modified, coverage in bed_rows(
            path, args.min_cov, args.max_cov
        ):
            ref = founder.get(key)
            if ref is None:
                continue
            if key in seen:
                raise ValueError(f"Duplicate {name} Total site: {key}")
            seen.add(key)

            ref_modified, ref_coverage, gene_ids = ref
            for gene in gene_ids:
                values = totals[gene]
                values[0] += 1
                values[1] += ref_modified
                values[2] += ref_coverage
                values[3] += modified
                values[4] += coverage

        rows = []
        counts = Counter()

        for gene in sorted(cds_genes):
            n, rm, rc, sm, sc = totals.get(
                gene, [0, 0.0, 0.0, 0.0, 0.0]
            )

            ref_fraction = rm / rc if rc else None
            sample_fraction = sm / sc if sc else None

            if n < args.min_sites:
                ref_status = sample_status = "insufficient_shared_sites"
                change = "insufficient_shared_sites"
            else:
                ref_status = status(
                    ref_fraction, args.um_max, args.methylated_min
                )
                sample_status = status(
                    sample_fraction, args.um_max, args.methylated_min
                )

                if "ambiguous" in (ref_status, sample_status):
                    change = "ambiguous_in_either"
                elif ref_status == "UM" and sample_status == "CG_methylated":
                    change = "gain"
                elif ref_status == "CG_methylated" and sample_status == "UM":
                    change = "loss"
                elif ref_status == "CG_methylated":
                    change = "retained_methylated"
                else:
                    change = "retained_UM"

            counts[change] += 1
            rows.append({
                "gene_id": gene,
                "shared_CG_sites": n,
                "J702_modified_sum": rm,
                "J702_coverage_sum": rc,
                "sample_modified_sum": sm,
                "sample_coverage_sum": sc,
                "J702_CG_fraction":
                    ref_fraction if ref_fraction is not None else "",
                "sample_CG_fraction":
                    sample_fraction if sample_fraction is not None else "",
                "difference_sample_minus_J702":
                    sample_fraction - ref_fraction
                    if ref_fraction is not None and sample_fraction is not None
                    else "",
                "J702_status": ref_status,
                "sample_status": sample_status,
                "change": change,
            })

        write_table(out / f"{name}.all_genes.tsv", rows, detail_fields)

        for category in ("gain", "loss", "retained_methylated"):
            write_table(
                out / f"{name}.{category}.tsv",
                [r for r in rows if r["change"] == category],
                detail_fields,
            )

        confident = sum(
            counts[c] for c in
            ("gain", "loss", "retained_methylated", "retained_UM")
        )

        summaries.append({
            "sample": name,
            "CDS_genes": len(cds_genes),
            "genes_passing_min_shared_sites":
                len(cds_genes) - counts["insufficient_shared_sites"],
            "confident_in_both": confident,
            **{
                category: counts[category]
                for category in summary_fields[4:]
            },
        })

        print(
            f"{name}: confident={confident:,}; "
            f"gain={counts['gain']:,}; loss={counts['loss']:,}; "
            f"common methylated={counts['retained_methylated']:,}",
            flush=True,
        )

    write_table(out / "summary.tsv", summaries, summary_fields)
    print(f"Saved results: {out.resolve()}")


if __name__ == "__main__":
    main()
