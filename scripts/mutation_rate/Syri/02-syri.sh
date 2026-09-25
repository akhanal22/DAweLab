#!/bin/bash
#SBATCH --job-name=syri_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --array=1-7
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL


set -euo pipefail

# load the module first
module load SyRI/1.7.1-foss-2023a
#module load plotsr/1.1.0
# --------------------
# 1. Get sample name
# --------------------
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linenameall.lst | tr -d '\r' | xargs)

# 2 Run Syri

syri -c ${file}_vs_J702_modified.sorted.paf -r ../J702_v3_modified_ragtag_re/flipped/J702_chr1-10.fa -q ../${file}_v3_modified_ragtag_re_after_telomere_rearrangement/${file}_chr1-10.fa -F P --prefix ${file}
#plotsr --sr B73D1_chr1_251_252Mb.zoom.syri.out --genomes genomes.txt -o B73D2_v2.png
#plotsr --sr B73D1syri.out --genomes genomes.txt --reg J702:chr1_RagTag_RagTag:250000000-260000000  -o b73d1_rough.png
#plotsr --sr B73D1syri.out --genomes genomes.txt -s 500 --reg J702:chr1_RagTag_RagTag:250000000-250500000  -o b73d1_nosyn.png
