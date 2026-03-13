This script is for generating the chromosome level dotplot. To begin with, I scaffolded the raw hifi assmebly with J702 reference genome. Then, genome to genome alignment is done to
generate the paf file.
1. genome-genome-alignment.sh
```
#!/bin/bash
#SBATCH --job-name=genome_genome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=400G
#SBATCH --array=1-7
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# --------------------
# Load environment
# --------------------
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

set -euo pipefail

# --------------------
# 1. Get sample name
# --------------------
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linename.lst | tr -d '\r' | xargs)

if [[ -z "$file" ]]; then
  echo "ERROR: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

echo "[$(date)] Processing sample: $file"

# --------------------
# 2. Input files
# --------------------
ref="J702_v3_modified_ragtag_re/ragtag.scaffold.fasta"
qry="${file}_v3_modified_ragtag_re/ragtag.scaffold.fasta"

outdir="paf_mapped_against_J702_modified"
outfile="${outdir}/${file}_vs_J702_modified.paf"

# --------------------
# 3. Sanity checks
# --------------------
[[ -s "$ref" ]] || { echo "ERROR: Missing reference $ref"; exit 1; }
[[ -s "$qry" ]] || { echo "ERROR: Missing query $qry"; exit 1; }

mkdir -p "$outdir"
# --------------------
# 4. Create samtools FASTA index (.fai)
# --------------------
#if [[ ! -s "${qry}.fai" ]]; then
 #   echo "[$(date)] Creating FASTA index (.fai)"
#    samtools faidx "$qry"
#fi

# --------------------
# 5. Create minimap2 index (.mmi)
# --------------------
#ref_index="${ref}.mmi"

#if [[ ! -s "$ref_index" ]]; then
 #   echo "[$(date)] Creating minimap2 index (.mmi)"
  #  minimap2 -d "$ref_index" "$ref"
#fi

# --------------------
# 6. Run minimap2 alignment
# --------------------
echo "[$(date)] Running minimap2..."

minimap2 \
  -x asm5 \
  -t "${SLURM_CPUS_PER_TASK}" \
  "$ref" "$qry" > "$outfile"

echo "[$(date)] Finished: $outfile"

```

After generating the paf file, i use following python script (grid_dotplot_by_chr.py) to generate chromosome wide dotplot for all individual and all chromosome.

```
#!/usr/bin/env python3

import os
import re
import math
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

PAF_COLS = [
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend",
    "nmatch", "alen", "mapq"
]

def read_paf(paf_file):
    return pd.read_csv(
        paf_file,
        sep="\t",
        header=None,
        usecols=range(12),
        names=PAF_COLS,
        low_memory=False
    )

def sample_name_from_file(paf_file):
    base = os.path.basename(paf_file)
    base = re.sub(r"_vs_.*\.paf$", "", base)
    return base

def chr_num(name):
    """
    Extract chromosome number from names like:
    chr10_RagTag_RagTag
    chr10_RagTag
    chr1
    """
    m = re.match(r"chr(\d+)", str(name))
    return int(m.group(1)) if m else None

def filter_for_chromosome(df, chrom_num, min_alen=0, min_mapq=0, primary_only=False):
    """
    Keep rows where query chromosome number == chrom_num
    and target chromosome number == chrom_num.

    This gives same-chromosome dotplots (e.g. chr10 vs chr10),
    which is usually what you want for clean per-chromosome plots.

    If you want off-chromosome hits too, I can give that variant too.
    """
    qchr = df["qname"].apply(chr_num)
    tchr = df["tname"].apply(chr_num)

    keep = (qchr == chrom_num) & (tchr == chrom_num)
    keep &= (df["alen"] >= min_alen)
    keep &= (df["mapq"] >= min_mapq)

    if primary_only:
        # optional: keep only primary alignments if tag exists in extra columns
        pass

    return df.loc[keep].copy()

def make_segments(df):
    fwd = []
    rev = []

    for _, row in df.iterrows():
        q1, q2 = row["qstart"], row["qend"]
        t1, t2 = row["tstart"], row["tend"]

        if row["strand"] == "+":
            fwd.append([(q1, t1), (q2, t2)])
        else:
            rev.append([(q1, t2), (q2, t1)])

    return fwd, rev

def plot_one_chromosome(paf_files, chrom_num, outdir, min_alen, min_mapq,
                        ncols=2, panel_size=6.5, linewidth=0.18,
                        same_scale=True, output_format="pdf"):
    samples = []
    dfs = []
    qmax_all = 0
    tmax_all = 0

    for paf in paf_files:
        sample = sample_name_from_file(paf)
        df = read_paf(paf)
        sub = filter_for_chromosome(df, chrom_num, min_alen=min_alen, min_mapq=min_mapq)

        samples.append(sample)
        dfs.append(sub)

        if not sub.empty:
            qmax_all = max(qmax_all, int(sub["qlen"].max()))
            tmax_all = max(tmax_all, int(sub["tlen"].max()))

    n = len(samples)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(panel_size * ncols, panel_size * nrows),
        squeeze=False
    )
    axes = axes.flatten()

    for ax, sample, df in zip(axes, samples, dfs):
        if df.empty:
            ax.set_title(f"{sample}\n(no chr{chrom_num} vs chr{chrom_num} hits)", fontsize=10)
            ax.set_xlabel(f"{sample} query chr{chrom_num} position")
            ax.set_ylabel(f"J702 target chr{chrom_num} position")
            continue

        fwd, rev = make_segments(df)

        if fwd:
            lc1 = LineCollection(fwd, linewidths=linewidth, alpha=0.9, rasterized=True)
            ax.add_collection(lc1)

        if rev:
            lc2 = LineCollection(rev, linewidths=linewidth, alpha=0.9, rasterized=True)
            ax.add_collection(lc2)

        if same_scale:
            ax.set_xlim(0, qmax_all)
            ax.set_ylim(0, tmax_all)
        else:
            ax.set_xlim(0, int(df["qlen"].max()))
            ax.set_ylim(0, int(df["tlen"].max()))

        ax.set_title(sample, fontsize=11)
        ax.set_xlabel(f"{sample} chr{chrom_num} query position (bp)")
        ax.set_ylabel(f"J702 chr{chrom_num} target position (bp)")
        ax.ticklabel_format(style="plain", axis="both")
        ax.grid(False)

    for i in range(len(samples), len(axes)):
        axes[i].axis("off")

    fig.suptitle(
        f"Chromosome {chrom_num} dotplot grid\n"
        f"Filter: chr{chrom_num} query vs chr{chrom_num} target | "
        f"min alignment length = {min_alen:,} bp | min MAPQ = {min_mapq}",
        fontsize=14
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"chr{chrom_num}_grid_dotplot.{output_format}")

    if output_format.lower() == "pdf":
        fig.savefig(outfile, bbox_inches="tight")
    else:
        fig.savefig(outfile, dpi=800, bbox_inches="tight")

    plt.close(fig)
    print(f"Saved: {outfile}")

def main():
    parser = argparse.ArgumentParser(description="Grid dotplot per chromosome from multiple PAF files")
    parser.add_argument("--paf_dir", required=True, help="Directory containing *.paf files")
    parser.add_argument("--pattern", default="*.paf", help="Glob pattern, default *.paf")
    parser.add_argument("--chromosomes", nargs="+", type=int, required=True,
                        help="Chromosome numbers to plot, e.g. 1 2 7 8 10")
    parser.add_argument("--outdir", default="grid_dotplots", help="Output directory")
    parser.add_argument("--min_alen", type=int, default=5000, help="Minimum alignment block length")
    parser.add_argument("--min_mapq", type=int, default=0, help="Minimum MAPQ")
    parser.add_argument("--ncols", type=int, default=2, help="Grid columns")
    parser.add_argument("--panel_size", type=float, default=6.5, help="Panel size")
    parser.add_argument("--linewidth", type=float, default=0.18, help="Line width")
    parser.add_argument("--same_scale", action="store_true", help="Use same axes across all individuals")
    parser.add_argument("--format", choices=["pdf", "png"], default="pdf", help="Output format")
    args = parser.parse_args()

    paf_files = sorted(glob.glob(os.path.join(args.paf_dir, args.pattern)))
    if not paf_files:
        raise SystemExit(f"No PAF files found in {args.paf_dir} with pattern {args.pattern}")

    print("Found PAF files:")
    for p in paf_files:
        print(" ", p)

    for chrom_num in args.chromosomes:
        plot_one_chromosome(
            paf_files=paf_files,
            chrom_num=chrom_num,
            outdir=args.outdir,
            min_alen=args.min_alen,
            min_mapq=args.min_mapq,
            ncols=args.ncols,
            panel_size=args.panel_size,
            linewidth=args.linewidth,
            same_scale=args.same_scale,
            output_format=args.format
        )

if __name__ == "__main__":
    main()
```

The submission scirpt (dotplot.sh)
```
#!/bin/bash
#SBATCH --job-name=chrgrid
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=chrgrid.%j.out
#SBATCH --error=chrgrid.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@uga.edu

set -euo pipefail

source ~/ENTER/etc/profile.d/conda.sh
conda activate base   # or your plotting/python env

python grid_dotplot_by_chr.py \
  --paf_dir . \
  --chromosomes 1 2 3 4 5 6 7 8 9 10 \
  --outdir chr_grid_plots \
  --min_alen 0 \
  --min_mapq 0 \
  --ncols 3 \
  --panel_size 7 \
  --linewidth 2 \
  --same_scale \
  --format pdf
```
