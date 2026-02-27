#!/bin/bash
#SBATCH --job-name=b73all_bam
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --array=1-1
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

#set -euo pipefail

# ---- activate conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

# 1) One BAM per array task
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" b73all.lst)

if [[ -z "${file}" ]]; then
  echo "ERROR: empty sample name for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

bam="${file}.sorted.bam"
if [[ ! -s "${bam}" ]]; then
  echo "ERROR: BAM not found or empty: ${bam}"
  exit 1
fi

# 2) Make per-sample output folders
outdir="./${file}"
mkdir -p "${outdir}/depth_primary" "${outdir}/depth_secondary" "${outdir}/plots"

# 3) Index BAM if needed (recommended)
if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
  samtools index -@ "${SLURM_CPUS_PER_TASK}" "${bam}"
fi

# 4) One depth file per sample (across ALL contigs)

# PRIMARY = exclude secondary (0x100) + supplementary (0x800) => 0x900
primary_out="${outdir}/depth_primary/${file}.primary.depth"
samtools view -@ "${SLURM_CPUS_PER_TASK}" -b -F 0x900 "${bam}" \
  | samtools depth -aa -@ "${SLURM_CPUS_PER_TASK}" - \
  > "${primary_out}"

# SECONDARY = include secondary (0x100), exclude supplementary (0x800)
secondary_out="${outdir}/depth_secondary/${file}.secondary.depth"
samtools view -@ "${SLURM_CPUS_PER_TASK}" -b -f 0x100 -F 0x800 "${bam}" \
  | samtools depth -aa -@ "${SLURM_CPUS_PER_TASK}" - \
  > "${secondary_out}"

echo "Done:"
echo "  PRIMARY:   ${primary_out}"
echo "  SECONDARY: ${secondary_out}"
