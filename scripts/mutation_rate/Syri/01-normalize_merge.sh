#!/bin/bash
#SBATCH --job-name=merge_syri
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=100:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL


module load SAMtools/1.23.1-GCC-13.3.0
module load BCFtools/1.23.1-GCC-13.3.0

for i in {1..10}; do
    printf 'chr%s\tchr%s_RagTag\n' "$i" "$i"
done > syri_to_J702_chromosomes.tsv

bash merge_syri_normalize.sh \
  /scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa \
  syri_samples.tsv \
  syri_normalized_results \
  syri_to_J702_chromosomes.tsv
