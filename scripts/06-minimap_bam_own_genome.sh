#!/bin/bash
#SBATCH --job-name=J702_minimap
#SBATCH --partition=highmem_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=400G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# 1. Pull the name directly from your list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linename.lst)

genome="../../raw_genome/${file}_v3.contig.fa"

# ---- 1) activate conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

# 2. I want to map the raw hifi reads against the raw reference genome for each lines, converted it to bam files and sorted it.
# -t 8: Uses the 8 CPUs you requested
# -o: Creates a unique folder for this specific file to avoid conflicts

minimap2 -t 8 -ax map-hifi "${genome}" ../../fastX/${file}.hifireads.fastq.gz \
| samtools view -@ 8 -b \
| samtools sort -@ 8 -o ${file}.sorted.bam


 samtools index ${file}.sorted.bam
