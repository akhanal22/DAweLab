#!/bin/bash
#SBATCH --job-name=AW_proali
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=120:00:00
#SBATCH --mem=100G
#SBATCH --array=1-3
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# --------------------------------------------------
# Activate AnchorWave environment
# --------------------------------------------------

#source ~/ENTER/etc/profile.d/conda.sh
#conda activate anchorwave

# --------------------------------------------------
# Input configuration
# --------------------------------------------------

LIST="threename.lst"

REF_DIR="./J702_v3_modified_ragtag_re"

REF_FASTA="${REF_DIR}/J702_flip.fa"
REF_GFF="${REF_DIR}/J702.liftoff.gff3"
REF_CDS="${REF_DIR}/J702.cds.fa"

REF_SAM="anchorwave_alignment/reference_mapping/J702_cds_to_J702.sam"

QUERY_SAM_DIR="anchorwave_alignment/query_mapping"
OUTDIR="anchorwave_alignment/maf_re"

mkdir -p "$OUTDIR"

# --------------------------------------------------
# Check required programs
# --------------------------------------------------

command -v anchorwave >/dev/null 2>&1 || {
    echo "ERROR: anchorwave was not found in the current environment."
    exit 1
}

# --------------------------------------------------
# Read sample
# --------------------------------------------------

[[ -s "$LIST" ]] || {
    echo "ERROR: Sample list not found or empty: $LIST"
    exit 1
}

sample=$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" |
    tr -d '\r' |
    xargs
)

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: $sample"

# J702 is the reference, not a query
if [[ "$sample" == "J702" ]]; then
    echo "J702 is the reference genome. Skipping this task."
    exit 0
fi

# --------------------------------------------------
# Query genome and SAM paths
# --------------------------------------------------

QUERY_DIR="${sample}_v3_modified_ragtag_re_after_telomere_rearrangement"
QUERY_FASTA="${QUERY_DIR}/${sample}_ragtag.scaffold.fasta"

QUERY_SAM="${QUERY_SAM_DIR}/J702_cds_to_${sample}.sam"

ANCHORS="${OUTDIR}/J702_vs_${sample}.anchors"
MAF="${OUTDIR}/J702_vs_${sample}.maf"
FRAGMENT_MAF="${OUTDIR}/J702_vs_${sample}.fragment.maf"
METHOD_BED="${OUTDIR}/J702_vs_${sample}.alignment_method.bed"
LOG="${OUTDIR}/J702_vs_${sample}.anchorwave.log"

# Temporary output files
TMP_ANCHORS="${ANCHORS}.tmp"
TMP_MAF="${MAF}.tmp"
TMP_FRAGMENT_MAF="${FRAGMENT_MAF}.tmp"
TMP_METHOD_BED="${METHOD_BED}.tmp"

# --------------------------------------------------
# Check all inputs
# --------------------------------------------------

for file in \
    "$REF_FASTA" \
    "$REF_GFF" \
    "$REF_CDS" \
    "$REF_SAM" \
    "$QUERY_FASTA" \
    "$QUERY_SAM"
do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Required input is missing or empty:"
        echo "$file"
        exit 1
    fi
done

# --------------------------------------------------
# Report inputs
# --------------------------------------------------

echo "Reference FASTA: $REF_FASTA"
echo "Reference GFF:   $REF_GFF"
echo "Reference CDS:   $REF_CDS"
echo "Reference SAM:   $REF_SAM"
echo "Query FASTA:     $QUERY_FASTA"
echo "Query SAM:       $QUERY_SAM"
echo "Output MAF:      $MAF"
echo "Threads:         $SLURM_CPUS_PER_TASK"
echo "Start time:      $(date)"

anchorwave --version 2>/dev/null || true

# --------------------------------------------------
# Run whole-genome alignment
# --------------------------------------------------
#
# -R 1 : maximum alignment depth on reference
# -Q 1 : maximum alignment depth on query
#
# For maize versus maize, 1-to-1 is the appropriate
# initial setting.
# --------------------------------------------------

anchorwave proali \
    -i "$REF_GFF" \
    -r "$REF_FASTA" \
    -as "$REF_CDS" \
    -ar "$REF_SAM" \
    -a "$QUERY_SAM" \
    -s "$QUERY_FASTA" \
    -n "$TMP_ANCHORS" \
    -R 1 \
    -Q 1 \
    -o "$TMP_MAF" \
    -f "$TMP_FRAGMENT_MAF" \
    -b "$TMP_METHOD_BED" \
    -t "$SLURM_CPUS_PER_TASK" \
    > "$LOG" 2>&1

# --------------------------------------------------
# Confirm that AnchorWave produced output
# --------------------------------------------------

if [[ ! -s "$TMP_MAF" ]]; then
    echo "ERROR: AnchorWave did not produce a nonempty MAF."
    echo "Check log: $LOG"
    exit 1
fi

if [[ ! -s "$TMP_ANCHORS" ]]; then
    echo "ERROR: AnchorWave did not produce a nonempty anchor file."
    echo "Check log: $LOG"
    exit 1
fi

# Rename outputs only after successful completion
mv "$TMP_ANCHORS" "$ANCHORS"
mv "$TMP_MAF" "$MAF"

if [[ -s "$TMP_FRAGMENT_MAF" ]]; then
    mv "$TMP_FRAGMENT_MAF" "$FRAGMENT_MAF"
fi

if [[ -s "$TMP_METHOD_BED" ]]; then
    mv "$TMP_METHOD_BED" "$METHOD_BED"
fi

# --------------------------------------------------
# Basic validation
# --------------------------------------------------

block_count=$(grep -c '^a' "$MAF" || true)
sequence_line_count=$(grep -c '^s' "$MAF" || true)
anchor_count=$(grep -cv '^#' "$ANCHORS" || true)

echo
echo "AnchorWave completed for $sample"
echo "MAF alignment blocks: $block_count"
echo "MAF sequence lines:   $sequence_line_count"
echo "Anchor-file rows:     $anchor_count"
echo
echo "Output files:"
ls -lh "$ANCHORS" "$MAF" "$LOG"

[[ -s "$FRAGMENT_MAF" ]] && ls -lh "$FRAGMENT_MAF"
[[ -s "$METHOD_BED" ]] && ls -lh "$METHOD_BED"

echo
echo "End time: $(date)"
