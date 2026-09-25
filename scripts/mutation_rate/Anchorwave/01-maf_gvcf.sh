#!/bin/bash
#SBATCH --job-name=maftovcf
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=100:00:00
#SBATCH --mem=500G
#SBATCH --array=1-7
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# 1. Pull the name directly from your list
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst)


# 2. maf to gvcf

 ./tassel-5-standalone/run_pipeline.pl -Xmx100G \
 -debug \
 -MAFToGVCFPlugin  \
 -referenceFasta ../../J702_v3_modified_ragtag_re/J702_flip.fa \
 -mafFile J702_vs_${line}.maf \
 -sampleName ${line} \
 -gvcfOutput ./${line}.gvcf \
 -fillGaps false > ./${line}GVCF.txt
