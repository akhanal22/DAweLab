#!/bin/bash
#SBATCH --job-name=nucplot_chr
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --array=1-10
#SBATCH --output=nucplot.%A_%a.out
#SBATCH --error=nucplot.%A_%a.err

set -euo pipefail

# ----------------------------
# USER SETTINGS (edit these)
# ----------------------------

# File listing chromosomes (one per line). SLURM_ARRAY_TASK_ID chooses the line.
CHR_LIST="chr.list"

# File listing samples (one per line), e.g.:
# B73D1
# B73D2
# ...
SAMPLES_LIST="linename.lst"

# Where the full-genome BAMs live.
# If BAMs are in the submit directory as ${sample}.sorted.bam, leave as ".".
BAM_DIR="."

# Name pattern for full-genome BAMs (coordinate-sorted):
# If they are "${sample}.sorted.bam", leave as-is.
BAM_PATTERN="{sample}.sorted.bam"

# Path to NucPlot.py (adjust to your location)
NUCPLOT="../NucFreq/NucPlot.py"

# NucPlot options (optional)
FONT_SIZE=16

# Montage layout: 2 columns x 4 rows for 8 samples
GRID_COLS=2
GRID_ROWS=4

# ----------------------------
# Derived vars / setup
# ----------------------------

CHR="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHR_LIST" | tr -d '\r' | xargs)"
if [[ -z "${CHR}" ]]; then
  echo "[ERROR] Could not read chromosome from ${CHR_LIST} at line ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

if [[ ! -s "$SAMPLES_LIST" ]]; then
  echo "[ERROR] Missing samples list: $SAMPLES_LIST"
  exit 1
fi

if [[ ! -x "$NUCPLOT" && ! -f "$NUCPLOT" ]]; then
  echo "[ERROR] NucPlot.py not found at: $NUCPLOT"
  exit 1
fi

# Output structure per chromosome
OUTDIR="nucplot_${CHR}"
CHR_BAM_DIR="${OUTDIR}/chr_bams"
PNG_DIR="${OUTDIR}/pngs"
GRID_DIR="${OUTDIR}/grids"
mkdir -p "$CHR_BAM_DIR" "$PNG_DIR" "$GRID_DIR"

echo "[INFO] Chromosome: $CHR"
echo "[INFO] Samples file: $SAMPLES_LIST"
echo "[INFO] Output dir: $OUTDIR"

# ----------------------------
# Helper: resolve BAM path
# ----------------------------
resolve_bam () {
  local sample="$1"
  local fn="${BAM_PATTERN/\{sample\}/$sample}"
  echo "${BAM_DIR%/}/${fn}"
}
# -------------------
# Modules (avoid conda/module mixing)
# -------------------
log "Loading modules"
module purge
unset PYTHONPATH

# Adjust to what exists on your cluster; the key is: ONE python stack only.
module load SAMtools/1.21-GCC-13.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load SciPy-bundle/2024.05-gfbf-2024a
module load Pysam/0.22.1-GCC-13.3.0
module load matplotlib
module load pandas/2.3.0
module load Seaborn/0.13.2-gfbf-2024a

# sanity
log "Python: $(which python)"
python -c "import sys; print(sys.version)"
python -c "import numpy, pysam, pandas, matplotlib; print('Imports OK')"


# ----------------------------
# Main loop: per sample
# ----------------------------
pngs_to_montage=()

while IFS= read -r sample || [[ -n "$sample" ]]; do
  sample="$(echo "$sample" | tr -d '\r' | xargs)"
  [[ -z "$sample" ]] && continue

  BAM="$(resolve_bam "$sample")"
  if [[ ! -s "$BAM" ]]; then
    echo "[ERROR] Missing BAM for sample $sample: $BAM"
    exit 1
  fi

  # Ensure BAM is indexed (full-genome BAM)
  if [[ ! -s "${BAM}.bai" ]]; then
    echo "[INFO] Indexing full BAM (missing .bai): $BAM"
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$BAM"
  fi

  CHR_BAM="${CHR_BAM_DIR}/${sample}.${CHR}.bam"
  CHR_PNG="${PNG_DIR}/${sample}.${CHR}.png"

  # 1) Create chr-level BAM only if missing/empty (and index)
  if [[ -s "$CHR_BAM" && -s "${CHR_BAM}.bai" ]]; then
    echo "[INFO] Chr BAM exists — skipping extraction: $CHR_BAM"
  else
    echo "[INFO] Extracting $CHR for $sample -> $CHR_BAM"
    tmp="${CHR_BAM}.tmp"
    samtools view -@ "$SLURM_CPUS_PER_TASK" -b "$BAM" "$CHR" > "$tmp"
    mv "$tmp" "$CHR_BAM"
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$CHR_BAM"
  fi

  # 2) Run NucPlot only if PNG missing/empty
  if [[ -s "$CHR_PNG" ]]; then
    echo "[INFO] PNG exists — skipping NucPlot: $CHR_PNG"
  else
    echo "[INFO] Running NucPlot for $sample on $CHR"
    # IMPORTANT: pass an explicit region to avoid "read whole bam" behavior
    # We ask NucPlot to use the whole chromosome but as an explicit region.
    # NucPlot's own log usually prints chr length if it parses header.
    "$NUCPLOT" "$CHR_BAM" "$CHR_PNG" --regions "${CHR}"
  fi

  pngs_to_montage+=("$CHR_PNG")

done < "$SAMPLES_LIST"

# ----------------------------
# Montage into a 2x4 grid (or as many as you have)
# ----------------------------
GRID_OUT="${GRID_DIR}/${CHR}.grid.${GRID_COLS}x${GRID_ROWS}.png"

if [[ -s "$GRID_OUT" ]]; then
  echo "[INFO] Grid already exists — skipping montage: $GRID_OUT"
else
  echo "[INFO] Creating grid: $GRID_OUT"

  # Check ImageMagick
  if ! command -v montage >/dev/null 2>&1; then
    echo "[ERROR] ImageMagick 'montage' not found. Load module or install it."
    exit 1
  fi

  # montage expects row-major ordering; we keep file order from linename.lst
  montage "${pngs_to_montage[@]}" \
    -tile "${GRID_COLS}x${GRID_ROWS}" \
    -geometry +10+10 \
    "$GRID_OUT"
fi

echo "[DONE] Chromosome $CHR complete."
echo "[DONE] Grid: $GRID_OUT"
