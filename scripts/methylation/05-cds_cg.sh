#!/bin/bash
#SBATCH --job-name=pbg_input
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=400G
#SBATCH --output=TE_CpG_profiles.%j.out
#SBATCH --error=TE_CpG_profiles.%j.err

python compare_shared_CDS.py \
    --gff J702.liftoff.gff3 \
    --founder J702.hifireads.pbmm2.combined.bed.gz \
    --sample B73D1=B73D1.hifireads.pbmm2.combined.bed.gz \
    --sample B73D2=B73D2.hifireads.pbmm2.combined.bed.gz \
    --sample B73D4=B73D4.hifireads.pbmm2.combined.bed.gz \
    --sample B73D7c=B73D7c.hifireads.pbmm2.combined.bed.gz \
    --sample B73D10=B73D10.hifireads.pbmm2.combined.bed.gz \
    --sample B73D7B=B73D7B.hifireads.pbmm2.combined.bed.gz \
    --sample B73all=B73all.hifireads.pbmm2.combined.bed.gz \
    --um-max 0.05 \
    --methylated-min 0.20
