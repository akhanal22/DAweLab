#!/usr/bin/env python3

import sys
import argparse
import matplotlib.pyplot as plt

def parse_paf(paf_file, query_chr=None, target_chr=None, min_len=10000):
    """
    Parse PAF and return alignment segments.
    Filters by chromosome name and minimum alignment length.
    """
    segments = []

    with open(paf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")

            qname = cols[0]
            qstart = int(cols[2])
            qend   = int(cols[3])
            strand = cols[4]
            tname = cols[5]
            tstart = int(cols[7])
            tend   = int(cols[8])
            aln_len = int(cols[10])

            if aln_len < min_len:
                continue

            if query_chr and qname != query_chr:
                continue

            if target_chr and tname != target_chr:
                continue

            segments.append((qstart, qend, tstart, tend, strand))

    return segments


def make_dotplot(segments, query_label, target_label, output):
    fig, ax = plt.subplots(figsize=(8, 8))

    for qstart, qend, tstart, tend, strand in segments:
        if strand == "+":
            ax.plot([qstart, qend], [tstart, tend], linewidth=0.5)
        else:
            ax.plot([qstart, qend], [tend, tstart], linewidth=0.5)

    ax.set_xlabel(query_label)
    ax.set_ylabel(target_label)
    ax.set_title("Dotplot")

    ax.ticklabel_format(style='plain')
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    print(f"Dotplot saved to {output}")


def main():
    parser = argparse.ArgumentParser(description="Generate dotplot from PAF file")
    parser.add_argument("--paf", required=True, help="Input PAF file")
    parser.add_argument("--query", help="Query chromosome name (optional)")
    parser.add_argument("--target", help="Target chromosome name (optional)")
    parser.add_argument("--min_len", type=int, default=10000,
                        help="Minimum alignment length to plot (default=10kb)")
    parser.add_argument("--out", required=True, help="Output PNG file")

    args = parser.parse_args()

    segments = parse_paf(args.paf, args.query, args.target, args.min_len)

    if not segments:
        print("No alignments found after filtering.")
        sys.exit(1)

    make_dotplot(
        segments,
        args.query if args.query else "Query",
        args.target if args.target else "Target",
        args.out
    )


if __name__ == "__main__":
    main()
