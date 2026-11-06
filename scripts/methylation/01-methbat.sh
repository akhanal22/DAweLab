#!/bin/bash
#SBATCH --job-name=pbcpg
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=100:00:00
#SBATCH --mem=400G
#SBATCH --output=%x.%J.out
#SBATCH --error=%x.%J.err
#SBATCH --array=1-8
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL


# ---- conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools



# ---- sample name for array submission----
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst | tr -d '\r' | xargs)

if [[ -z "${file}" ]]; then
  echo "ERROR: empty sample name for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

#----Load genome file for indexing and outdir----
#genome="../J702_v3_modified_ragtag_re/J702_flip.fa"

#1 lets index the genome file now-----
#pbmm2 index $genome $genome.mmi

#2 after genome index, we should align the fastq pacbio reads to the genome
#pbmm2 align $genome /scratch/ak18787/B73_project/CCS/CCS/${file}.merged.bam ${file}.hifireads.bam

#3, In this step, I index the bam file
#samtools sort -@ 8 -m 2G -o methylation_after_manual_curation_t2t/${file}.hifireads.sorted.bam methylation_after_manual_curation_t2t/${file}.hifireads.bam
#samtools index methylation_after_manual_curation_t2t/${file}.hifireads.sorted.bam

#4 before running the next step lets activate methylation conda environment
conda activate methylation

#4 In this step, the bam file is converted to .pbmm2 using following code. This essentially summarize site probabilities for 5mC methylation.

#aligned_bam_to_cpg_scores \
 #--bam methylation_after_manual_curation_t2t/${file}.hifireads.sorted.bam \
 #--output-prefix methylation_after_manual_curation_t2t/${file}.hifireads.pbmm2 \
  #--threads 16

#5 In this step, i do the segmentation step with methbat to get the blocks that are methylated
methbat segment \
  --input-prefix ${file}.hifireads.pbmm2 \
  --output-prefix ${file}.hifireads.methbat
