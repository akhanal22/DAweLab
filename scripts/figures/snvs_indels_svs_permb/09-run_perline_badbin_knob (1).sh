#!/bin/bash
#SBATCH --job-name=perline_knob
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# 1000 shuffles per MA line, each line using only its own SyRI bad bins
cd /scratch/ak18787/B73_project/final_vcf/permutation
ENV=$HOME/ENTER/envs/py310_edta/bin
$ENV/python 09-region_enrichment_perline_badbin_knob.py \
    --bedtools $ENV/bedtools \
    --n-shuffle 1000 \
    --threads $SLURM_CPUS_PER_TASK
