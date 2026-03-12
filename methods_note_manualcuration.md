I need to manually curate the genome because the size of chromosome 1 is shorter and chromosome 2is longer for J702 and B73D4. When aligned against B73D5 reference, some parts of chromosom0e 
2 are perfectly syntenic to chromosome 1 of reference. Also, the extended part is at the junction of raw contigs arrangement. The coverega is also low in the junction region.

1. For J702:
```
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
```

2. For B73D2. Chromosome 1 is longer (310Mb) , i foundthat contig ptg000025l has some parts that align to chr2 on refernce genome and junction have reduction in read depth so i am breaking up contgs. 
```
#SBATCH --job-name=genome_genome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

FA=B73D2_v3.contig.fa
OUT=B73D2_v3.contig_modified.fa

set -euo pipefail

# ensure index exists
samtools faidx "$FA"

# get contig length
LEN025=$(awk '$1=="ptg000025l"{print $2}' "$FA.fai")

# split contig
samtools faidx "$FA" ptg000025l:1-1659726 \
 | sed 's/^>ptg000025l/>ptg000025la/' > "$OUT"

samtools faidx "$FA" ptg000025l:1659727-"$LEN025" \
 | sed 's/^>ptg000025l/>ptg000025lb/' >> "$OUT"

# keep all other contigs
cut -f1 "$FA.fai" | grep -Ev '^ptg000025l$' > keep.list
samtools faidx -r keep.list "$FA" >> "$OUT"

rm keep.list

echo "Modified FASTA written to $OUT"
```
3) For B73d7c
```
#!/bin/bash
#SBATCH --job-name=genome_genome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# ---- conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

FA=B73D7c_v3.contig.fa
OUT=B73D7c_v3.contig_modified.fa

set -euo pipefail
# get contig length automatically
LEN=$(awk '$1=="ptg000002l"{print $2}' "$FA.fai")

# split contig
samtools faidx "$FA" ptg000002l:1-5352613 \
  | sed 's/^>ptg000002l/>ptg000002la/' > "$OUT"

samtools faidx "$FA" ptg000002l:5352614-$LEN \
  | sed 's/^>ptg000002l/>ptg000002lb/' >> "$OUT"

# keep all other contigs
cut -f1 "$FA.fai" | grep -v '^ptg000002l$' > keep.list
samtools faidx -r keep.list "$FA" >> "$OUT"

rm keep.list
echo "Modified FASTA written to $OUT"
```
4) For B73d1
```
#!/bin/bash
#SBATCH --job-name=b73d1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# ---- conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools
set -euo pipefail
FA="B73D1_v3.contig.fa"
OUT="B73D1_v3_modified.contig.fa"

# index fasta if needed
samtools faidx "$FA"

# get contig lengths from .fai
LEN020=$(awk '$1=="ptg000020l"{print $2}' "$FA.fai")
LEN014=$(awk '$1=="ptg000014l"{print $2}' "$FA.fai")

echo "ptg000020l length = $LEN020"
echo "ptg000014l length = $LEN014"

# 1) keep all contigs except the ones being replaced
awk '
BEGIN{skip=0}
/^>/{
    name=$1
    sub(/^>/,"",name)
    if(name=="ptg000020l" || name=="ptg000014l") skip=1
    else skip=0
}
!skip
' "$FA" > kept_contigs.fa

# 2) split ptg000020l into a, b, c
samtools faidx "$FA" ptg000020l:0-1022450 \
  | sed '1s/.*/>ptg000020la/' > ptg000020la.fa

samtools faidx "$FA" ptg000020l:1022451-2524178 \
  | sed '1s/.*/>ptg000020lb/' > ptg000020lb.fa

samtools faidx "$FA" ptg000020l:2524179-$LEN020 \
  | sed '1s/.*/>ptg000020lc/' > ptg000020lc.fa

# 3) split ptg000014l into a, b
samtools faidx "$FA" ptg000014l:0-9274519 \
  | sed '1s/.*/>ptg000014la/' > ptg000014la.fa

samtools faidx "$FA" ptg000014l:9274520-$LEN014 \
  | sed '1s/.*/>ptg000014lb/' > ptg000014lb.fa

# 4) combine: unchanged contigs + replacement pieces
cat kept_contigs.fa \
    ptg000020la.fa ptg000020lb.fa ptg000020lc.fa \
    ptg000014la.fa ptg000014lb.fa > "$OUT"

# 5) index final fasta
samtools faidx "$OUT"

echo "Done. Output fasta: $OUT"
```
5) For B73D4
```
#!/bin/bash
#SBATCH --job-name=genome_genome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

FA=B73D4_v3.contig.fa
OUT=B73D4_v3.contig_modified.fa

set -euo pipefail

# ensure index exists
samtools faidx "$FA"

# get contig length
LEN025=$(awk '$1=="ptg000004l"{print $2}' "$FA.fai")

# split contig
samtools faidx "$FA" ptg000004l:1-195746274 \
 | sed 's/^>ptg000004l/>ptg000004la/' > "$OUT"

samtools faidx "$FA" ptg000004l:195746275-"$LEN025" \
 | sed 's/^>ptg000004l/>ptg000004lb/' >> "$OUT"

# keep all other contigs
cut -f1 "$FA.fai" | grep -Ev '^ptg000004l$' > keep.list
samtools faidx -r keep.list "$FA" >> "$OUT"

rm keep.list
samtools faidx "$OUT"
echo "Modified FASTA written to $OUT"
```
After breaking up the contig, I reran ragtag using B73D5 as reference. However, with that still B73D1 had not resolved issue for chromsome 7 longer. This could be because the repeat 
length is shorter in NAM version. So, I first did scaffolding for J792 using B73_V5 genome and used J702 for scaffolding others
```
#!/bin/bash
#SBATCH --job-name=d10
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=500G
#SBATCH --array=1-1
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

source ~/ENTER/etc/profile.d/conda.sh
conda activate base

set -euo pipefail

# sample name
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" d10.lst | tr -d '\r' | xargs)
[[ -n "$file" ]] || { echo "ERROR: empty sample name for task $SLURM_ARRAY_TASK_ID"; exit 1; }

# output directory
outdir="${file}_v3_modified_ragtag_re"
mkdir -p "$outdir"

# run ragtag
ragtag.py scaffold J702_v3_modified_ragtag/ragtag.scaffold.fasta raw_genome/${file}_v3.contig_modified.fa -o "$outdir"
```



