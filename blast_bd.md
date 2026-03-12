I used the consensus.dimer.fasta to generate the blast output for TR1knob, CentC, Knob180 and generate the bed file for each of them as follows:
```
#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --partition=batch
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --array=1-8

#set -euo pipefail

# --- conda / tools ---
module purge
export CONDA_NO_PLUGINS=true
source ~/ENTER/etc/profile.d/conda.sh
conda activate blast_env   # must contain bedtools + awk etc

LIST="allname.lst"
REPEATS=("CentC" "Cent4" "knob180" "TR1")
MIN_ALN=100   # same as your screenshot ($4>=30)

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r' | xargs)
[[ -n "$SAMPLE" ]] || { echo "ERROR: empty sample for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

FOLDER="${SAMPLE}_v3_modified_ragtag_re"
OUTDIR="${FOLDER}/repeat_blast"
BLASTOUT="${OUTDIR}/${SAMPLE}.blast.tsv"

#[[ -s "$BLASTOUT" ]] || { echo "ERROR: missing/empty BLAST table: $BLASTOUT" >&2; exit 1; }

mkdir -p "$OUTDIR"

command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
# ----------------------------
# Step 1: BLAST repeats to genome
# ----------------------------

QUERY="consensus.dimer.fasta"
GENOME="${SAMPLE}_v3_modified_ragtag_re/ragtag.scaffold.fasta"

# make blast db if it doesn't exist
if [[ ! -f "${GENOME}.nin" && ! -f "${GENOME}.nhr" ]]; then
    echo "Building BLAST database for $GENOME"
    makeblastdb -in "$GENOME" -dbtype nucl
fi

echo "Running BLAST..."

blastn \
  -query "$QUERY" \
  -db "$GENOME" \
  -out "$BLASTOUT" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads $SLURM_CPUS_PER_TASK \
  -evalue 1e-10
# ----------------------------
# Step 2: intervals + spans
# ----------------------------
for rep in "${REPEATS[@]}"; do
  intervals="${OUTDIR}/intervals_${SAMPLE}_${rep}.bed"
  summary="${OUTDIR}/span_${SAMPLE}_${rep}.tsv"

  # Avoid killing the job if a repeat isn't present
  if ! grep -Fq "$rep" "$BLASTOUT"; then
    : > "$intervals"
    : > "$summary"
    continue
  fi

  # 2A) merged intervals BED4: contig start end rep
  grep -F "$rep" "$BLASTOUT" \
    | awk '$3>=90 && $4>=100' \
    | cut -f2,9,10 \
    | awk 'BEGIN{OFS="\t"}{ if($3<$2) print $1,$3,$2; else print $1,$2,$3 }' \
    | bedtools sort -i - \
    | bedtools merge -i - \
    | awk -v r="$rep" 'BEGIN{OFS="\t"}{print $1,$2,$3,r}' \
    > "$intervals"

  # 2B) per-contig span + merged_bp
  # contig  repeat  min_start  max_end  span_bp  merged_bp
  awk -v r="$rep" 'BEGIN{OFS="\t"}
    {
      c=$1; s=$2; e=$3;
      if(!(c in min) || s<min[c]) min[c]=s;
      if(!(c in max) || e>max[c]) max[c]=e;
      merged[c] += (e-s);
    }
    END{
      for(c in min){
        span = max[c]-min[c];
        print c, r, min[c], max[c], span, merged[c];
      }
    }' "$intervals" | sort -k1,1 > "$summary"
done

# Combined outputs
cat "${OUTDIR}/intervals_${SAMPLE}_"*.bed | sort -k1,1 -k2,2n > "${OUTDIR}/all_intervals_${SAMPLE}.bed"

{
  echo -e "contig\trepeat\tmin_start\tmax_end\tspan_bp\tmerged_bp"
  cat "${OUTDIR}/span_${SAMPLE}_"*.tsv
} > "${OUTDIR}/all_span_${SAMPLE}.tsv"

echo "DONE:"
echo "  ${OUTDIR}/all_intervals_${SAMPLE}.bed"
echo "  ${OUTDIR}/all_span_${SAMPLE}.tsv"
```
After this to generate stacked plot for chr1, chr2, chr7 and chr8. I use following python script named as lineplot.py. You could also do it trhough submission script
```
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

# =========================================================
# USER SETTINGS
# =========================================================
SAMPLES = [
    "B73D1_v3_modified_ragtag_re",
    "B73D2_v3_modified_ragtag_re",
    "B73D4_v3_modified_ragtag_re",
    "B73D7B_v3_modified_ragtag_re",
    "B73D7c1_v3_modified_ragtag_re",
    "B73D10_v3_modified_ragtag_re",
    "B73all_v3_modified_ragtag_re",
    "J702_v3_modified_ragtag_re",
]

SUBDIR = "repeat_blast"
CHROMOSOMES = ["chr1", "chr2", "chr7", "chr8"]

REPEAT_COLORS = {
    "TR1": "#f03b91",
    "knob180": "#32cd32",
    "CentC": "#ffcc33",
    "Cent4": "#ff8c00",
}

TRACK_COLOR = "#1f77a8"      # blue chromosome bar
AGP_FACE = "#d9d9d9"         # gray AGP blocks
AGP_EDGE = "#333333"

OUTDIR = "repeat_landscape_by_chr"
DPI = 300

# =========================================================
# HELPERS
# =========================================================
def sample_prefix(sample):
    return re.sub(r"_v3_modified_ragtag_re$", "", sample)

def sample_label(sample):
    p = sample_prefix(sample)
    p = p.replace("c1", "c")
    return f"{p}_v3"

def normalize_chr(name):
    """
    Normalize:
      chr8
      chr8_RagTag
      chr8_RagTag_RagTag
    all to:
      chr8
    """
    name = str(name).strip()
    name = re.sub(r'(_ragtag)+$', '', name, flags=re.IGNORECASE)
    return name.lower()

def load_repeat_bed(sample):
    prefix = sample_prefix(sample)
    bedfile = os.path.join(sample, SUBDIR, f"all_intervals_{prefix}.bed")

    if not os.path.exists(bedfile):
        print(f"WARNING: missing BED file: {bedfile}")
        return pd.DataFrame(columns=["chrom", "start", "end", "repeat", "chrom_norm"])

    df = pd.read_csv(
        bedfile,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "repeat"],
        dtype={"chrom": str, "start": int, "end": int, "repeat": str}
    )

    df["chrom_norm"] = df["chrom"].apply(normalize_chr)
    return df

def load_agp(sample):
    agpfile = os.path.join(sample, "ragtag.scaffold.agp")

    if not os.path.exists(agpfile):
        print(f"WARNING: missing AGP file: {agpfile}")
        return pd.DataFrame(columns=[
            "object", "object_beg", "object_end", "part_number", "component_type",
            "component_id", "component_beg", "component_end", "orientation", "object_norm"
        ])

    rows = []
    with open(agpfile) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue

            obj = parts[0]
            obj_beg = int(parts[1])
            obj_end = int(parts[2])
            part_number = parts[3]
            ctype = parts[4]

            component_id = parts[5] if len(parts) > 5 else ""
            component_beg = parts[6] if len(parts) > 6 and parts[6].isdigit() else None
            component_end = parts[7] if len(parts) > 7 and parts[7].isdigit() else None
            orientation = parts[8] if len(parts) > 8 else ""

            rows.append([
                obj, obj_beg, obj_end, part_number, ctype,
                component_id, component_beg, component_end, orientation
            ])

    df = pd.DataFrame(rows, columns=[
        "object", "object_beg", "object_end", "part_number", "component_type",
        "component_id", "component_beg", "component_end", "orientation"
    ])
    df["object_norm"] = df["object"].apply(normalize_chr)
    return df

def chromosome_max_length(chrom, repeat_data, agp_data):
    mx = 1
    for sample in SAMPLES:
        rdf = repeat_data[sample]
        adf = agp_data[sample]

        subr = rdf[rdf["chrom_norm"] == chrom]
        if not subr.empty:
            mx = max(mx, int(subr["end"].max()))

        suba = adf[adf["object_norm"] == chrom]
        if not suba.empty:
            mx = max(mx, int(suba["object_end"].max()))

    return mx

def draw_repeat_track(ax, df, chrom, y0, bar_h=0.34):
    sub = df[df["chrom_norm"] == chrom].copy()
    if sub.empty:
        return

    max_end = sub["end"].max()
    ax.add_patch(Rectangle((0, y0), max_end, bar_h, facecolor=TRACK_COLOR, edgecolor=TRACK_COLOR, lw=0.8))

    for _, row in sub.iterrows():
        rep = row["repeat"]
        color = REPEAT_COLORS.get(rep)
        if color is None:
            continue
        width = max(1, row["end"] - row["start"])
        ax.add_patch(Rectangle((row["start"], y0), width, bar_h, facecolor=color, edgecolor=color, lw=0))

def draw_agp_track(ax, agp_df, chrom, y0, bar_h=0.16):
    sub = agp_df[agp_df["object_norm"] == chrom].copy()
    if sub.empty:
        return

    # plot only W components for contig structure
    sub = sub[sub["component_type"] == "W"]
    if sub.empty:
        return

    for _, row in sub.iterrows():
        x = row["object_beg"] - 1
        w = row["object_end"] - row["object_beg"] + 1
        ax.add_patch(Rectangle((x, y0), w, bar_h, facecolor=AGP_FACE, edgecolor=AGP_EDGE, lw=0.8))

# =========================================================
# LOAD DATA
# =========================================================
os.makedirs(OUTDIR, exist_ok=True)

repeat_data = {sample: load_repeat_bed(sample) for sample in SAMPLES}
agp_data = {sample: load_agp(sample) for sample in SAMPLES}

# =========================================================
# MAKE ONE FIGURE PER CHROMOSOME
# =========================================================
for chrom in CHROMOSOMES:
    chr_max = chromosome_max_length(chrom, repeat_data, agp_data)

    fig_h = max(8, len(SAMPLES) * 1.25)
    fig, ax = plt.subplots(figsize=(16, fig_h))

    y_gap = 1.2
    repeat_h = 0.38
    agp_h = 0.18

    top_y = len(SAMPLES) * y_gap

    for i, sample in enumerate(SAMPLES):
        y = top_y - (i + 1) * y_gap

        rdf = repeat_data[sample]
        adf = agp_data[sample]

        # repeat bar
        draw_repeat_track(ax, rdf, chrom, y0=y + 0.28, bar_h=repeat_h)

        # agp bar underneath
        draw_agp_track(ax, adf, chrom, y0=y + 0.08, bar_h=agp_h)

        # sample label
        ax.text(
            -chr_max * 0.01,
            y + 0.46,
            sample_label(sample),
            ha="right",
            va="center",
            fontsize=14
        )

    ax.set_xlim(0, chr_max * 1.03)
    ax.set_ylim(0, top_y + 0.9)

    ax.set_yticks([])
    ax.set_xlabel(f"Position on {chrom} (Mb)", fontsize=14)
    ax.set_title(f"{chrom} repeat landscape (stacked)", fontsize=17)

    xticks = list(range(0, int(chr_max / 1_000_000) + 20, 10))
    ax.set_xticks([x * 1_000_000 for x in xticks])
    ax.set_xticklabels([str(x) for x in xticks], fontsize=12)

    legend_handles = [
        Patch(facecolor=REPEAT_COLORS["TR1"], edgecolor=REPEAT_COLORS["TR1"], label="TR1"),
        Patch(facecolor=REPEAT_COLORS["knob180"], edgecolor=REPEAT_COLORS["knob180"], label="knob180"),
        Patch(facecolor=REPEAT_COLORS["CentC"], edgecolor=REPEAT_COLORS["CentC"], label="CentC"),
        Patch(facecolor=AGP_FACE, edgecolor=AGP_EDGE, label="AGP contig"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=13)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    plt.tight_layout()

    outfile = os.path.join(OUTDIR, f"{chrom}_repeat_landscape_stacked.png")
    fig.savefig(outfile, dpi=DPI, bbox_inches="tight")
    print("Saved:", outfile)
    plt.close(fig)
```






