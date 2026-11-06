#!/usr/bin/env python3

from pathlib import Path
from itertools import combinations
import pandas as pd

# Lines to compare. Remove a line here if you do not want it included.
samples = [
    "B73D1",
    "B73D2",
    "B73D4",
    "B73D7c",
    "B73D10",
    "B73D7B",
    "B73all",
]

input_dir = Path(".")
output_dir = Path("shared_gain_loss_comparisons")
output_dir.mkdir(exist_ok=True)


def read_gene_ids(path):
    """Read unique gene IDs from a gain or loss TSV."""
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    df = pd.read_csv(path, sep="\t")

    if "gene_id" not in df.columns:
        raise ValueError(f"'gene_id' column missing in {path}")

    return set(df["gene_id"].dropna().astype(str))


def compare(change):
    """Compare gain or loss gene lists across all samples."""

    gene_sets = {}

    for sample in samples:
        path = input_dir / f"{sample}.{change}.tsv"
        gene_sets[sample] = read_gene_ids(path)

    # Every gene found in at least one sample.
    all_genes = sorted(set().union(*gene_sets.values()))

    rows = []

    for gene in all_genes:
        present = [sample for sample in samples if gene in gene_sets[sample]]

        row = {
            "gene_id": gene,
            "number_of_lines": len(present),
            "lines": ",".join(present),
        }

        for sample in samples:
            row[sample] = gene in gene_sets[sample]

        rows.append(row)

    columns = ["gene_id", "number_of_lines", "lines"] + samples
    matrix = pd.DataFrame(rows, columns=columns)

    if not matrix.empty:
        matrix = matrix.sort_values(
            ["number_of_lines", "gene_id"],
            ascending=[False, True],
        )

    # All gains/losses and the lines where each appears.
    matrix.to_csv(
        output_dir / f"{change}_presence_matrix.tsv",
        sep="\t",
        index=False,
    )

    # Genes shared by at least two lines.
    matrix[matrix["number_of_lines"] >= 2].to_csv(
        output_dir / f"{change}_shared_at_least_2_lines.tsv",
        sep="\t",
        index=False,
    )

    # Genes shared by every line in the samples list.
    matrix[matrix["number_of_lines"] == len(samples)].to_csv(
        output_dir / f"{change}_shared_all_lines.tsv",
        sep="\t",
        index=False,
    )

    # Every pair of lines.
    pair_rows = []

    for sample1, sample2 in combinations(samples, 2):
        common = sorted(gene_sets[sample1] & gene_sets[sample2])

        pair_rows.append({
            "line_1": sample1,
            "line_2": sample2,
            "line_1_total": len(gene_sets[sample1]),
            "line_2_total": len(gene_sets[sample2]),
            f"shared_{change}_genes": len(common),
            "shared_gene_ids": ",".join(common),
        })

    pairwise = pd.DataFrame(pair_rows)

    pairwise.to_csv(
        output_dir / f"{change}_pairwise_comparisons.tsv",
        sep="\t",
        index=False,
    )

    print(f"\n{change.upper()}")
    for sample in samples:
        print(f"  {sample}: {len(gene_sets[sample])} genes")

    print(
        "  Shared by at least 2 lines:",
        int((matrix["number_of_lines"] >= 2).sum())
    )
    print(
        "  Shared by all lines:",
        int((matrix["number_of_lines"] == len(samples)).sum())
    )


compare("gain")
compare("loss")

print(f"\nResults saved in: {output_dir.resolve()}")
