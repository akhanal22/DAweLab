#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle


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


def plot_dotplot_with_te(df, ref_te, qry_te, start, end, output):
    fig = plt.figure(figsize=(12, 9))

    # Layout:
    # left TE panel | main dotplot
    #                bottom TE panel
    ax_dot = fig.add_axes([0.24, 0.38, 0.68, 0.55])
    ax_ref = fig.add_axes([0.24, 0.14, 0.68, 0.14], sharex=ax_dot)
    ax_qry = fig.add_axes([0.04, 0.38, 0.12, 0.55])

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

    # -----------------------
    # Main dotplot
    # -----------------------
    segments = []

    for _, row in df.iterrows():
        x1 = row["ref_start_genome"] / 1e6
        x2 = row["ref_end_genome"] / 1e6

        if row["strand"] == "+":
            y1 = row["qry_start_genome"] / 1e6
            y2 = row["qry_end_genome"] / 1e6
        else:
            y1 = row["qry_end_genome"] / 1e6
            y2 = row["qry_start_genome"] / 1e6

        segments.append([(x1, y1), (x2, y2)])

    lc = LineCollection(segments, linewidths=0.35)
    ax_dot.add_collection(lc)

    ax_dot.set_xlim(start / 1e6, end / 1e6)
    ax_dot.set_ylim(start / 1e6, end / 1e6)

    # Show chromosome coordinate tick labels on BOTH axes
    ax_dot.tick_params(axis="x", labelbottom=True)
    ax_dot.tick_params(axis="y", labelleft=True)

    ax_dot.set_xlabel("J702 chr1 position (Mb)")
    ax_dot.set_ylabel("B73all chr1 position (Mb)", labelpad=18)

    ax_dot.set_title(
        f"chr1:{start/1e6:.1f}-{end/1e6:.1f} Mb dotplot with TE annotation"
    )

    # -----------------------
    # J702 TE track below x-axis
    # -----------------------
    for _, row in ref_te.iterrows():
        c = simplify_class(row["classification"])
        y = class_y.get(c, class_y["Other"])
        x = row["start"] / 1e6
        w = (row["end"] - row["start"]) / 1e6

        ax_ref.add_patch(Rectangle((x, y - 0.35), w, 0.7))

    ax_ref.set_xlim(start / 1e6, end / 1e6)
    ax_ref.set_ylim(-1, len(classes))

    ax_ref.set_yticks(range(len(classes)))
    ax_ref.set_yticklabels(classes, fontsize=8)

    ax_ref.set_ylabel("J702\nTEs", fontsize=9)
    ax_ref.set_xlabel("J702 chr1 position (Mb)")

    # -----------------------
    # B73all TE track on left side
    # -----------------------
    for _, row in qry_te.iterrows():
        c = simplify_class(row["classification"])
        x = class_y.get(c, class_y["Other"])
        y = row["start"] / 1e6
        h = (row["end"] - row["start"]) / 1e6

        ax_qry.add_patch(Rectangle((x - 0.35, y), 0.7, h))

    ax_qry.set_xlim(-1, len(classes))
    ax_qry.set_ylim(start / 1e6, end / 1e6)

    # TE names for left panel
    ax_qry.set_xticks(range(len(classes)))
    ax_qry.set_xticklabels(classes, rotation=90, fontsize=8)

    # Y-axis coordinates for B73all TE panel
    ax_qry.tick_params(axis="y", labelleft=True)
    ax_qry.set_ylabel("B73all chr1 position (Mb)", fontsize=9)

    ax_qry.set_xlabel("B73all TEs", fontsize=9)
    ax_qry.set_title("B73all\nTEs", fontsize=9)

    plt.savefig(output, dpi=600, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot windowed PAF dotplot with EDTA TE annotation tracks."
    )

    parser.add_argument("--paf", required=True)
    parser.add_argument("--ref_gff", required=True)
    parser.add_argument("--qry_gff", required=True)

    parser.add_argument("--ref_chr", default="chr1_RagTag_RagTag")
    parser.add_argument("--qry_chr", default="chr1_RagTag_RagTag")

    parser.add_argument("--start", type=int, default=227000000)
    parser.add_argument("--end", type=int, default=250000000)

    parser.add_argument("--output", default="chr1_dotplot_TE.png")

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
        args.output
    )


if __name__ == "__main__":
    main()
