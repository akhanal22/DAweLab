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

# 1. Split ptg000016l at 23 700 222 bp. 
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
Even after manual curationb, for B63D2, 251-256mb of chromosome 1 is inverted
```
<img width="1097" height="1089" alt="image" src="https://github.com/user-attachments/assets/b9c55315-c4aa-49b6-9371-17e0686268eb" />

```
so to check whether this inversion is real or assembly artifact, i first looked at agp files and found that that region actually belongs to the contig ptg000025lb:1659727-29128499 that I broke earlier for manual curation of B73D2. the bam on both end 251mb and 256mb showed the coverage dip which led me conclude that this might be assembly artificat. While it is common to have dip on the edges of raw contig the coverage dip on the middle where we expect uniform coverage (251mb in our case) suggest assembly artifacts.   
chr1_RagTag_RagTag      1       138479420       1       W       ptg000018l      1       138479420       -
chr1_RagTag_RagTag      138479421       138479520       2       U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      138479521       138505638       3       W       ptg000760l      1       26118   -
chr1_RagTag_RagTag      138505639       138505738       4       U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      138505739       227754277       5       W       ptg000017l      1       89248539        -
chr1_RagTag_RagTag      227754278       227754377       6       U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      227754378       227782828       7       W       ptg000655l      1       28451   -
chr1_RagTag_RagTag      227782829       227782928       8       U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      227782929       227843963       9       W       ptg000142l      1       61035   -
chr1_RagTag_RagTag      227843964       227844063       10      U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      227844064       229351416       11      W       ptg000048l      1       1507353 -
chr1_RagTag_RagTag      229351417       229351516       12      U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      229351517       256820289       13      W       ptg000025lb:1659727-29128499    1       27468773        -
chr1_RagTag_RagTag      256820290       256820389       14      U       100     scaffold        yes     align_genus
chr1_RagTag_RagTag      256820390       309279454       15      W       ptg000007l      1       52459065        -
To exactly break the contig based on where i see inverted orintation, I mapped raw fasta for ptg000025lb:1659727-29128499 against J702 and the result for that contig is as follows
ptg000025lb:1659727-29128499	27468773	5393322	27464957	-	chr1_RagTag	309250393	229464846	251536282	16992826	22072226	60	tp:A:P	cm:i:1603102	s1:i:16992540	s2:i:14046	dv:f:0.0000	rl:i:9035536
ptg000025lb:1659727-29128499	27468773	13082	5374650	+	chr1_RagTag	309250393	251543717	256905353	4062584	5361711	60	tp:A:P	cm:i:382330	s1:i:4062522	s2:i:46912	dv:f:0.0000	rl:i:9035536. Based on this pattern, following submission script named 13-b73d2_rearrangement_with_j702.sh was submitted for manual curation.
```
#!/bin/bash
#SBATCH --job-name=j702
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# input and output file names
FA=B73D2_v3.contig_modified.fa
OUT=B73D2_v3.contig_modified_re.fa

# make sure fasta is indexed
samtools faidx "$FA"

# split ptg000025lb into two contigs
samtools faidx "$FA" ptg000025lb:1-5374650 \
  | sed 's/^>ptg000025lb:1-5374650$/>ptg000025lba/' > "$OUT"

samtools faidx "$FA" ptg000025lb:5374651-27468773 \
  | sed 's/^>ptg000025lb:5374651-27468773$/>ptg000025lbb/' >> "$OUT"

# extract all other contigs except the original ptg000025lb
cut -f1 "$FA.fai" | grep -Ev '^(ptg000025lb)$' > keep.list
samtools faidx -r keep.list "$FA" >> "$OUT"

# clean up
rm keep.list
```
I had to rearrange B73D7B genome too because based on syri vcf it looks likes there is inversion "chr1_RagTag_RagTag	227949567	INVAL33	N	<INVAL>	.	PASS	END=229457231;ChrB=chr1_RagTag_RagTag;StartB=229572573;EndB=228064848;Parent=INV12;VarType=.;DupType=.	GT	1" where a is J702 and B is B73D2. But, the inversion breakposint has no coverage when i reads of B73D7b against its own genome. So, I am manually breaking the contig at the exact inversion breakpoint to redo scaffolding
```
#!/bin/bash
#SBATCH --job-name=b73d7b
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# input and output file names
FA=B73D7B_v3.contig_modified.fa
OUT=B73D7B_v3.contig_modified_re.fa

# make sure fasta is indexed
samtools faidx "$FA"

# start fresh
: > "$OUT"

# split ptg000024l into three contigs
samtools faidx "$FA" ptg000024l:1-134805 \
  | sed 's/^>ptg000024l:1-134805$/>ptg000024la/' > "$OUT"

samtools faidx "$FA" ptg000024l:134806-1631659 \
  | sed 's/^>ptg000024l:134806-1631659$/>ptg000024lb/' >> "$OUT"

samtools faidx "$FA" ptg000024l:1631660-23733356 \
  | sed 's/^>ptg000024l:1631660-23733356$/>ptg000024lc/' >> "$OUT"

# extract all other contigs except the original ptg000024l
cut -f1 "$FA.fai" | grep -Ev '^(ptg000024l)$' > keep.list
samtools faidx -r keep.list "$FA" >> "$OUT"

# index new fasta
samtools faidx "$OUT"

# clean up
rm keep.list
```

