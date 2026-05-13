#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle


def clean_chr_name(chrom):
    return chrom.replace("_RagTag_RagTag", "")


def parse_window_name(name):
    m = re.match(r"(.+):(\d+)-(\d+)$", name)
    if m:
        return m.group(1), int(m.group(2))
    return name, 0


def read_windowed_paf(paf, ref_chr, qry_chr, start, end):
    cols = [
        "qry", "qry_len", "qry_start", "qry_end", "strand",
        "ref", "ref_len", "ref_start", "ref_end",
        "matches", "block_len", "mapq"
    ]

    df = pd.read_csv(paf, sep="\t", header=None, usecols=range(12), names=cols)

    qry_parsed = df["qry"].apply(parse_window_name)
    df["qry_chr_real"] = qry_parsed.apply(lambda x: x[0])
    df["qry_offset"] = qry_parsed.apply(lambda x: x[1])

    ref_parsed = df["ref"].apply(parse_window_name)
    df["ref_chr_real"] = ref_parsed.apply(lambda x: x[0])
    df["ref_offset"] = ref_parsed.apply(lambda x: x[1])

    df["qry_start_genome"] = df["qry_offset"] + df["qry_start"]
    df["qry_end_genome"] = df["qry_offset"] + df["qry_end"]

    df["ref_start_genome"] = df["ref_offset"] + df["ref_start"]
    df["ref_end_genome"] = df["ref_offset"] + df["ref_end"]

    df = df[
        (df["ref_chr_real"] == ref_chr) &
        (df["qry_chr_real"] == qry_chr) &
        (df["ref_end_genome"] >= start) &
        (df["ref_start_genome"] <= end) &
        (df["qry_end_genome"] >= start) &
        (df["qry_start_genome"] <= end)
    ].copy()

    return df


def read_gff(gff, chrom, start, end):
    records = []

    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            seqid, source, feature, s, e, score, strand, phase, attrs = parts
            s = int(s)
            e = int(e)

            if seqid != chrom:
                continue
            if e < start or s > end:
                continue

            classification = "Other"
            for item in attrs.split(";"):
                if item.lower().startswith("classification="):
                    classification = item.split("=", 1)[1]

            records.append({
                "chrom": seqid,
                "start": max(s, start),
                "end": min(e, end),
                "feature": feature,
                "strand": strand,
                "classification": classification
            })

    return pd.DataFrame(records)


def simplify_class(x):
    x = str(x)

    if "LTR/Gypsy" in x:
        return "LTR/Gypsy"
    if "LTR/Copia" in x:
        return "LTR/Copia"
    if "LTR" in x:
        return "LTR/other"
    if "DNA" in x or "TIR" in x or "DTC" in x or "DTA" in x or "DTH" in x:
        return "DNA/TIR"
    if "Helitron" in x or "DHH" in x:
        return "Helitron"
    if "LINE" in x:
        return "LINE"
    if "SINE" in x:
        return "SINE"

    return "Other"


def plot_dotplot_with_te(
    df, ref_te, qry_te,
    start, end, output,
    ref_chr, qry_chr,
    ref_name, qry_name,
    x_genome
):
    fig = plt.figure(figsize=(12, 9))

    ax_dot = fig.add_axes([0.24, 0.38, 0.68, 0.55])
    ax_xte = fig.add_axes([0.24, 0.14, 0.68, 0.14], sharex=ax_dot)
    ax_yte = fig.add_axes([0.04, 0.38, 0.12, 0.55], sharey=ax_dot)

    classes = [
        "LTR/Gypsy",
        "LTR/Copia",
        "LTR/other",
        "DNA/TIR",
        "Helitron",
        "LINE",
        "SINE",
        "Other"
    ]

    class_y = {c: i for i, c in enumerate(classes)}

    if x_genome == "query":
        x_start_col = "qry_start_genome"
        x_end_col = "qry_end_genome"
        y_start_col = "ref_start_genome"
        y_end_col = "ref_end_genome"

        x_te = qry_te
        y_te = ref_te

        x_name = qry_name
        y_name = ref_name

        x_chr = qry_chr
        y_chr = ref_chr

    else:
        x_start_col = "ref_start_genome"
        x_end_col = "ref_end_genome"
        y_start_col = "qry_start_genome"
        y_end_col = "qry_end_genome"

        x_te = ref_te
        y_te = qry_te

        x_name = ref_name
        y_name = qry_name

        x_chr = ref_chr
        y_chr = qry_chr

    segments = []

    for _, row in df.iterrows():
        x1 = row[x_start_col] / 1e6
        x2 = row[x_end_col] / 1e6

        if row["strand"] == "+":
            y1 = row[y_start_col] / 1e6
            y2 = row[y_end_col] / 1e6
        else:
            y1 = row[y_end_col] / 1e6
            y2 = row[y_start_col] / 1e6

        segments.append([(x1, y1), (x2, y2)])

    lc = LineCollection(segments, linewidths=0.35)
    ax_dot.add_collection(lc)

    ax_dot.set_xlim(start / 1e6, end / 1e6)
    ax_dot.set_ylim(start / 1e6, end / 1e6)

    ax_dot.tick_params(axis="x", labelbottom=True)
    ax_dot.tick_params(axis="y", labelleft=True)

    short_chr = clean_chr_name(x_chr)

    ax_dot.set_xlabel(f"{x_name} {short_chr} position (Mb)")
    ax_dot.set_ylabel(f"{y_name} {clean_chr_name(y_chr)} position (Mb)", labelpad=18)

    ax_dot.set_title(
        f"{short_chr}:{start/1e6:.1f}-{end/1e6:.1f} Mb dotplot with TE annotation"
    )

    # X-axis TE track
    for _, row in x_te.iterrows():
        c = simplify_class(row["classification"])
        y = class_y.get(c, class_y["Other"])
        x = row["start"] / 1e6
        w = (row["end"] - row["start"]) / 1e6
        ax_xte.add_patch(Rectangle((x, y - 0.35), w, 0.7))

    ax_xte.set_xlim(start / 1e6, end / 1e6)
    ax_xte.set_ylim(-1, len(classes))
    ax_xte.set_yticks(range(len(classes)))
    ax_xte.set_yticklabels(classes, fontsize=8)
    ax_xte.set_ylabel(f"{x_name}\nTEs", fontsize=9)
    ax_xte.set_xlabel(f"{x_name} {clean_chr_name(x_chr)} position (Mb)")

    # Y-axis TE track
    for _, row in y_te.iterrows():
        c = simplify_class(row["classification"])
        x = class_y.get(c, class_y["Other"])
        y = row["start"] / 1e6
        h = (row["end"] - row["start"]) / 1e6
        ax_yte.add_patch(Rectangle((x - 0.35, y), 0.7, h))

    ax_yte.set_xlim(-1, len(classes))
    ax_yte.set_ylim(start / 1e6, end / 1e6)

    ax_yte.set_xticks(range(len(classes)))
    ax_yte.set_xticklabels(classes, rotation=90, fontsize=8)

    ax_yte.tick_params(axis="y", labelleft=True)
    ax_yte.set_ylabel(f"{y_name} {clean_chr_name(y_chr)} position (Mb)", fontsize=9)
    ax_yte.set_xlabel(f"{y_name} TEs", fontsize=9)
    ax_yte.set_title(f"{y_name}\nTEs", fontsize=9)

    plt.savefig(output, dpi=600, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot windowed PAF dotplot with EDTA TE annotation tracks."
    )

    parser.add_argument("--paf", required=True)
    parser.add_argument("--ref_gff", required=True)
    parser.add_argument("--qry_gff", required=True)

    parser.add_argument("--ref_chr", required=True)
    parser.add_argument("--qry_chr", required=True)

    parser.add_argument("--ref_name", default="J702")
    parser.add_argument("--qry_name", default="B73all")

    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)

    parser.add_argument(
        "--x_genome",
        choices=["ref", "query"],
        default="query",
        help="Which genome to put on x-axis. Use query for B73all/B73D1 on x-axis and J702 on y-axis."
    )

    parser.add_argument("--output", default="dotplot_TE.png")

    args = parser.parse_args()

    paf_df = read_windowed_paf(
        args.paf,
        args.ref_chr,
        args.qry_chr,
        args.start,
        args.end
    )

    ref_te = read_gff(
        args.ref_gff,
        args.ref_chr,
        args.start,
        args.end
    )

    qry_te = read_gff(
        args.qry_gff,
        args.qry_chr,
        args.start,
        args.end
    )

    print(f"PAF alignments retained: {len(paf_df)}")
    print(f"Reference TE records retained: {len(ref_te)}")
    print(f"Query TE records retained: {len(qry_te)}")

    plot_dotplot_with_te(
        paf_df,
        ref_te,
        qry_te,
        args.start,
        args.end,
        args.output,
        args.ref_chr,
        args.qry_chr,
        args.ref_name,
        args.qry_name,
        args.x_genome
    )


if __name__ == "__main__":
    main()
