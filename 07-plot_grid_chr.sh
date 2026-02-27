#!/bin/bash
#SBATCH --job-name=chr_plot
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --array=1-10
#SBATCH --output=depth_chr.%A_%a.out
#SBATCH --error=depth_chr.%A_%a.err

chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chr.txt | tr -d '\r' | xargs)
mkdir -p chr_plots

python3 plot_chr_grid.py \
  --samples samples8.lst \
  --chr "$chr" \
  --depth_pattern "./{sample}/depth/{sample}_q20_2308primary.depth" \
  --bin_bp 10000 \
  --sharey \
  --out_png "chr_plots/${chr}.8samples_primaryq20_2308.grid.png" \
  --title "${chr} depth (8 individuals; 2x4 panels)"
