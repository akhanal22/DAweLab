#!/bin/bash
#SBATCH --job-name=hifiasm_re
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=100:00:00
#SBATCH --mem=500G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# 1. Pull the name directly from your list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linename.lst)

# --- config ---
OUTDIR_v3="/scratch/ak18787/B73_project/HiFi_assembly_v3"
READDIR="/scratch/ak18787/B73_project/fastX"
THREADS=20

cd /scratch/ak18787/B73_project

ml SAMtools/1.21-GCC-13.3.0
module load hifiasm/0.24.0-GCCcore-12.3.0

# --- check & prepare directories ---
if [ ! -d "$OUTDIR_v3" ]; then
  echo "[*] Creating main output directory: $OUTDIR_v3"
  mkdir -p "$OUTDIR_v3"
else
  echo "[✓] Main output directory already exists: $OUTDIR_v3"
fi

if [ ! -d "$OUTDIR_v4" ]; then
  echo "[*] Creating main output directory: $OUTDIR_v4"
  mkdir -p "$OUTDIR_v4"
else
  echo "[✓] Main output directory already exists: $OUTDIR_v4"
fi

# --- run hifiasm ---
  FQ="$READDIR/${file}.hifireads.fastq.gz"

  echo "=== Running hifiasm v3 for $LAB ==="
  hifiasm --write-ec --write-paf -u 0 \
    -o "$OUTDIR_v3/${file}_v3.asm" \
    -l0 -t "$THREADS" \
    "$FQ"
done
