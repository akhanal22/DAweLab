When creating the dotplot for whole genome and whole chromosomes, it doesnot show the small insertions and deletion. So, we need to separate the fasta sequence and here we separate it in
1kb bins. To do that, first I create the bed file for each chromosome of each individual.
```
while read chr
do
    echo "Processing $chr"

    len=$(grep -w "$chr" ../../B73D2_v3_modified_ragtag_re_v2/chr1_10_only.fa.fai | cut -f2)

    echo -e "${chr}\t${len}" > ${chr}.genome

    bedtools makewindows \
        -g ${chr}.genome \
        -w 1000 \
        > B73D2/${chr}_1kb_b73d2.bed

done < ./chr1_10.list
```
Do this for each individual. We now use bedtools get fasta to break the fasta sequence for each individual. Ideally, you could do this is array job however naming should be consistent. I, unfortunately didnot have consistent naming for genome file so did it separately
```
#!/bin/bash
#SBATCH --job-name=b73d7c_scaffold
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

#set -euo pipefail

#source ~/ENTER/etc/profile.d/conda.sh
conda activate base

#--- 1)--- sample name
#file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" names.txt | tr -d '\r' | xargs)

#--- 2)--- get fasta for each chromosomes
while read chr
do
    echo "Extracting FASTA windows for $chr"

    bedtools getfasta \
        -fi ../../B73D4_v3_modified_ragtag_re/chromosome.fa \
        -bed B73D4/${chr}_1kb_b73d4.bed \
        -fo B73D4/${chr}_1kb_b73d4.fa

   bedtools getfasta \
        -fi ../../B73D7c1_v3_modified_ragtag_re/chromosome.fa \
        -bed B73D7c/${chr}_1kb_b73d7c.bed \
        -fo B73D7c/${chr}_1kb_b73d7c.fa

   bedtools getfasta \
        -fi ../../B73D10_v3_modified_ragtag_re/chromosome.fa \
        -bed B73D10/${chr}_1kb_b73d10.bed \
        -fo B73D10/${chr}_1kb_b73d10.fa

done < ./chr1_10.list

```
Note that I have shown it for three individuals only here, because i did for others separately. 
