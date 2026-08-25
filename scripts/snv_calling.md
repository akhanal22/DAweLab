
This is the code used for calling snvs using long read pacbio hifi data, I already had the bam file that was made using minimap2 by mapping raw reads against the J702 modified reference genomes.
So I used longcall (Version: 0.0.8) to run following script 
```
#!/bin/bash
#SBATCH --job-name=longcall
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/ENTER/etc/profile.d/conda.sh
conda activate base

file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linename.lst | tr -d '\r' | xargs)

ref="../J702_v3_modified_ragtag_re/J702_flip.fa"
bam="${file}.sorted.bam"
outdir="vcf_against_j702_snvs_longcall"

mkdir -p "$outdir"

echo "Sample: $file"
echo "Reference: $ref"
echo "BAM: $bam"
echo "Threads: $SLURM_CPUS_PER_TASK"
date

[[ -s "$ref" ]] || { echo "ERROR: missing reference $ref"; exit 1; }
[[ -s "$bam" ]] || { echo "ERROR: missing BAM $bam"; exit 1; }

if [[ ! -s "${bam}.bai" ]]; then
    samtools index "$bam"
fi

longcallD call \
  -t ${SLURM_CPUS_PER_TASK} \
  "$ref" \
  "$bam" \
  > "${outdir}/${file}.longcallD.vcf"
  ```
  After this i only select for snps
  
  ```
  bcftools view -v snps B73D1.longcallD.vcf -Oz -o B73D1.longcallD.snps.vcf
```

Normalize according to J702 positions
```
#!/bin/bash
#SBATCH --job-name=longcall_J702
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --array=1-6
#SBATCH --output=deepvariant_%A_%a.out
#SBATCH --error=deepvariant_%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL

source ~/ENTER/etc/profile.d/conda.sh
conda activate base
module load BCFtools/1.23.1-GCC-13.3.0
# reference used for mapping/calling
ref="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"

# sample name from list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst | tr -d '\r' | xargs)


gzvcf="${file}.longcallD.snps.vcf.gz"
normvcf="${file}.longcallD.snps.norm.vcf.gz"

echo "Sample: $file"
#echo "Input VCF: $invcf"
echo "Reference: $ref"
date

# checks
[[ -s "$ref" ]] || { echo "ERROR: missing reference $ref"; exit 1; }
[[ -s "${ref}.fai" ]] || samtools faidx "$ref"

# compress if needed
#if [[ ! -s "$gzvcf" ]]; then
 #   [[ -s "$invcf" ]] || { echo "ERROR: missing input VCF $invcf"; exit 1; }
  #  bgzip -c "$invcf" > "$gzvcf"
#fi

# index input
if [[ ! -s "${gzvcf}.csi" ]]; then
    bcftools index -c "$gzvcf"
fi

# normalize
bcftools norm \
  -f "$ref" \
  -Oz \
  -o "$normvcf" \
"$gzvcf"

# index normalized output
bcftools index -c "$normvcf"

echo "Finished normalization for $file"
date
```

  Now, the rigourous filtering only for confident snps

  ```
   samples="B73D1 B73D2 B73D4 B73D7c B73D10 J702"

for sample in $samples
do
    echo "Processing ${sample}..."

    bcftools filter \
        -e 'FILTER!="PASS" || FORMAT/DP<40 || FORMAT/DP>100 || FORMAT/GQ<30' \
        "${sample}.longcallD.snps.norm.vcf.gz" \
        -Ou |
    bcftools view \
        -i '(
              GT="1|1" &&
              FORMAT/AD[0:1]/FORMAT/DP >= 0.90
            ) || (
              (GT="0|1" || GT="1|0") &&
              FORMAT/AD[0:1]/FORMAT/DP >= 0.30 &&
              FORMAT/AD[0:1]/FORMAT/DP <= 0.70
            )' \
        -Oz \
        -o "${sample}.longcallD.snps.DP40-100.GQ20.rigorous.vcf.gz"

    bcftools index -f -t \
        "${sample}.longcallD.snps.DP40-100.GQ20.rigorous.vcf.gz"

    echo -n "${sample} retained SNPs: "
    bcftools view -H \
        "${sample}.longcallD.snps.DP40-100.GQ20.rigorous.vcf.gz" |
        wc -l

done


  ```

 Remove the indels that is shared with J702
```
while read file; do
    invcf="${file}.longcallD.snps.DP40-100.GQ20.rigorous.vcf.gz"
    outvcf="${file}.longcallD.snps.DP40-100.GQ20.rigorous.noJ702.vcf.gz"

    echo "Removing J702-shared SNPs from $file"
    echo "Input:  $invcf"
    echo "Output: $outvcf"

    bcftools isec \
      -n=1 -w1 \
      "$invcf" \
      J702.longcallD.snps.norm.vcf.gz \
      -Oz -o "$outvcf"

    bcftools index -c "$outvcf"

done < name.lst

```
MErge all and select unique snps
```
Merge all and select the unique indels only

module load BCFtools/1.21-GCC-13.3.0
bcftools merge   -l vcf.list   -Oz   -o five_lines.merged.vcf.gz # vcflist is the list of vcf files after rigorous filtering
bcftools index -c five_lines.merged.vcf.gz
bcftools view   -m2 -M2 -v indels   -i 'COUNT(GT="alt")=1'   -Oz   -o private_snps.only_one_line.vcf.gz   five_lines.merged.vcf.gz
bcftools index -c private_snps.only_one_line.vcf.gz
 awk '
BEGIN{OFS="\t"}
/^#CHROM/{
  for(i=10;i<=NF;i++) sample[i]=$i
  next
}
!/^#/{
  for(i=10;i<=NF;i++){
    gt=$i
    split(gt,a,":")
    if(a[1]=="0/1" || a[1]=="1/0" || a[1]=="1/1" || a[1]=="0|1" || a[1]=="1|0" || a[1]=="1|1"){
      count[sample[i]]++
    }
  }
}
END{
  for(s in count) print s, count[s]
}' private_snps.only_one_line.vcf > private_snps_counts.tsv
```
But i counted from chr1-10 only so i removed the vcf called for unscaffolded region
```
grep -v "ptg" private_snps.only_one_line.vcf | grep -v "scaf" > private_snps.only_chr1-10.vcf
```
then used the same approach above to count
```
awk '
BEGIN{OFS="\t"}
/^#CHROM/{
  for(i=10;i<=NF;i++) sample[i]=$i
  next
}
!/^#/{
  for(i=10;i<=NF;i++){
    gt=$i
    split(gt,a,":")
    if(a[1]=="0/1" || a[1]=="1/0" || a[1]=="1/1" || a[1]=="0|1" || a[1]=="1|0" || a[1]=="1|1"){
      count[sample[i]]++
    }
  }
}
END{
  for(s in count) print s, count[s]
}' private_snps.only_chr1-10.vcf> private_snps_counts_ch1-10.tsv
```
now, for genome size estimation, i use all the read with depth between 40-100 and only primary reads with mapq of 20 from chromosome 1-10 only.
```
#!/bin/bash
#SBATCH --job-name=array_minimap
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# 1. Pull the name directly from your list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" all.lst)

# ---- 1) activate conda ----
conda activate samtools



samtools depth \
    -@ "$SLURM_CPUS_PER_TASK" \
    -d 0 \
    -b J702_chr1_chr10.bed \
    "${file}.sorted.primary_mapq20.bam" |
awk '$3 > 40 && $3 < 100 {callable++}
     END {
         print "Callable genome size:", callable+0
     }' \
> "${file}.chr1_chr10.DP41_99.callable_size.txt"
```
This gives following genome size and mutation
Ind	 SNVs	Gen	genome_size	Mutation
	 (unique)			rate
B73D1	3437	33	2060981176	5.05E-08
B73D2	3565	33	2113955231	5.11E-08
B73D7c	3526	33	2084179119	5.13E-08
B73D10	3360	33	2107497801	4.83E-08
B73D4	3733	33	2121867275	5.33E-08

Now, same method can be used for indels too
