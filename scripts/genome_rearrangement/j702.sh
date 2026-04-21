#!/bin/bash
#SBATCH --job-name=j702
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# input and output file names
FA=J702_v3.contig.fa
OUT=J702_v3.contig_modified.fa

# 1. Split ptg000016l at 23 700 225 bp. note that cutoff starts at 122094  23700222. I just did three bp extra on right end.
samtools faidx "$FA" ptg000016l:0-23700222 \
  | sed 's/^>ptg000016l/>ptg000016la/' > "$OUT"
samtools faidx "$FA" ptg000016l:23700223-71528896 \
  | sed 's/^>ptg000016l/>ptg000016lb/' >> "$OUT"

# 2. Split ptg000016l at 5 377 852 bp (adjust contig name if needed)
samtools faidx "$FA" ptg000008l:1-5377840 \
  | sed 's/^>ptg000008l/>ptg000008la/' >> "$OUT"
samtools faidx "$FA" ptg000008l:5377841-200930298 \
  | sed 's/^>ptg000008l/>ptg000008lb/' >> "$OUT"

# 3. Extract all other contigs except the original ptg000016l
#    (build a list of contigs to keep from the FAI index)
cut -f1 "$FA.fai" | grep -Ev '^(ptg000016l|ptg000008l)$' > keep.list
samtools faidx -r keep.list "$FA" >> "$OUT"

# Clean up the temporary list
rm keep.list
