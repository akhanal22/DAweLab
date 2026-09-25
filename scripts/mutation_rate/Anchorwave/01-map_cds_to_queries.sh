#!/bin/bash
#SBATCH --job-name=AW_query_CDS
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=60G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# --------------------------------------------------
# Load minimap2
# --------------------------------------------------

# Option 1: module
# module purge
# module load Minimap2

# --------------------------------------------------
# Input configuration
# --------------------------------------------------

LIST="allname.lst"

REF_DIR="J702_v3_modified_ragtag_re"
REF_CDS="${REF_DIR}/J702.cds.fa"

OUTDIR="anchorwave_alignment/query_mapping"
mkdir -p "$OUTDIR"

# --------------------------------------------------
# Check list and reference CDS
# --------------------------------------------------

[[ -s "$LIST" ]] || {
    echo "ERROR: Sample list not found or empty: $LIST"
    exit 1
}

[[ -s "$REF_CDS" ]] || {
    echo "ERROR: Reference CDS FASTA not found: $REF_CDS"
    exit 1
}

command -v minimap2 >/dev/null 2>&1 || {
    echo "ERROR: minimap2 was not found."
    exit 1
}

# --------------------------------------------------
# Read sample name
# --------------------------------------------------

sample=$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" |
    tr -d '\r' |
    xargs
)

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Sample: $sample"

# J702 is the reference and was handled in Step 3
if [[ "$sample" == "J702" ]]; then
    echo "J702 is the reference genome."
    echo "Skipping it in the query-genome mapping step."
    exit 0
fi

# --------------------------------------------------
# Construct query genome path
# --------------------------------------------------

QUERY_DIR="${sample}_v3_modified_ragtag_re_after_telomere_rearrangement"
QUERY_FASTA="${QUERY_DIR}/ragtag.scaffold.fasta"

QUERY_SAM="${OUTDIR}/J702_cds_to_${sample}.sam"

# --------------------------------------------------
# Check query assembly
# --------------------------------------------------

[[ -d "$QUERY_DIR" ]] || {
    echo "ERROR: Query directory not found:"
    echo "$QUERY_DIR"
    exit 1
}

[[ -s "$QUERY_FASTA" ]] || {
    echo "ERROR: Query genome not found or empty:"
    echo "$QUERY_FASTA"
    exit 1
}

# --------------------------------------------------
# Report job
# --------------------------------------------------

echo "Reference CDS: $REF_CDS"
echo "Query genome:  $QUERY_FASTA"
echo "Output SAM:    $QUERY_SAM"
echo "Threads:       $SLURM_CPUS_PER_TASK"
echo "Start time:    $(date)"

minimap2 --version

# --------------------------------------------------
# Map J702 CDS sequences to query genome
# --------------------------------------------------

minimap2 \
    -x splice \
    -k 12 \
    -a \
    -t "$SLURM_CPUS_PER_TASK" \
    "$QUERY_FASTA" \
    "$REF_CDS" \
    > "${QUERY_SAM}.tmp"

# Keep final name only after successful completion
mv "${QUERY_SAM}.tmp" "$QUERY_SAM"

# --------------------------------------------------
# Validate output
# --------------------------------------------------

[[ -s "$QUERY_SAM" ]] || {
    echo "ERROR: SAM output is empty: $QUERY_SAM"
    exit 1
}

alignment_count=$(awk '$1 !~ /^@/ && $3 != "*" {count++} END {print count+0}' "$QUERY_SAM")

echo "Mapped SAM records: $alignment_count"
echo "Output size:"
ls -lh "$QUERY_SAM"

if [[ "$alignment_count" -eq 0 ]]; then
    echo "ERROR: No J702 CDS sequences mapped to $sample."
    exit 1
fi

echo "Query CDS mapping completed successfully for $sample."
echo "End time: $(date)"
