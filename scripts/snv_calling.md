<img width="919" height="551" alt="image" src="https://github.com/user-attachments/assets/69bb1c9d-7129-42c0-8b11-5e1c122efef8" />This is the code used for calling snvs using long read pacbio hifi data, I already had the bam file that was made using minimap2 by mapping raw reads against the J702 modified reference genomes.
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

ref="../J702_v3_modified_ragtag_re/ragtag.scaffold.fasta"
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
  Now, the rigourous filtering only for confident snps
  ```
   bcftools filter   -e 'QUAL<50 || FILTER!="PASS" || FORMAT/DP<40 || FORMAT/GQ<30'   B73D1.longcallD.snps.vcf -Ou | bcftools view   -i '(
        GT="1|1" && FORMAT/AD[0:1]/FORMAT/DP >= 0.90
      ) || (
        (GT="0|1" || GT="1|0") &&
        FORMAT/AD[0:1]/FORMAT/DP >= 0.30 &&
        FORMAT/AD[0:1]/FORMAT/DP <= 0.70
      )'   -Oz -o B73D1.longcallD.snps.rigorous.vcf

  ```

Same filtering was applied for all other individuals
Now to be sure that the variant are accurate and not generated due to sequencing error, i filtered out the snvs that were present in J702 using following code but first i normalized the vcf for each of them
```
#SBATCH --array=1-7
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/ENTER/etc/profile.d/conda.sh
conda activate base

# reference used for mapping/calling
ref="../../J702_v3_modified_ragtag_re/ragtag.scaffold.fasta"

# sample name from list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linename.lst | tr -d '\r' | xargs)


gzvcf="${file}.longcallD.snps.rigorous.vcf.gz"
normvcf="${file}.longcallD.snps.rigorous.norm.vcf.gz"

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

After normalization, I finally do the filtering
```
while read file; do
    invcf="${file}.longcallD.snps.rigorous.norm.vcf.gz"
    outvcf="${file}.longcallD.snps.rigorous.norm.noJ702.vcf.gz"

    echo "Removing J702-shared SNPs from $file"
    echo "Input:  $invcf"
    echo "Output: $outvcf"

    bcftools isec \
      -n=1 -w1 \
      "$invcf" \
      J702.longcallD.snps.norm.vcf.gz \
      -Oz -o "$outvcf"

    bcftools index -c "$outvcf"

done < linename.lst

```
For estimating mutation rate, we usually look for prive snvs, because we want to focus on new mutation that arose independently per generation, If it is common across individual, it could easily just be background mutation present in the founder line i.e. J702. So, I generate the merge vcf file first using bcftools and count the private snvs as below:
```
module load BCFtools/1.21-GCC-13.3.0
bcftools merge   -l vcf.list   -Oz   -o seven_lines.merged.vcf.gz # vcflist is the list of vcf files after rigorous filtering
bcftools index -c seven_lines.merged.vcf.gz
bcftools view   -m2 -M2 -v snps   -i 'COUNT(GT="alt")=1'   -Oz   -o private_snvs.only_one_line.vcf.gz   seven_lines.merged.vcf.gz
bcftools index -c private_snvs.only_one_line.vcf.gz
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
    if(a[=="0/1" || a[1]=="1/0" || a[1]=="1/1" || a[1]=="0|1" || a[1]=="1|0" || a[1]=="1|1"){
      count[sample[i]]++
    }
  }
}
END{
  for(s in count) print s, count[s]
}' private_snvs.only_one_line.vcf > private_snv_counts.tsv
```
Based on this following table represents private snv counts for each

```
| Line | Private SNVs | Generations |
|------|---------------|-------------|
| B73D1 | 3520 | 33 |
| B73D2 | 3580 | 33 |
| B73D4 | 398 | 33 |
| B73D7B | 2142 | 33 |
| B73D7c | 2289 | 33 |
| B73D10 | 3478 | 33 |
| B73all | 770 | 35 |
| **Total** | **16177** | **233** |
```
Now we calculate the mutation rate as below:

For estimating mutation  rate:
Total unique SNVs (sum all) = 16177
Total generation = 6 *33 + 35 = 233
Genome size = 2.1gb
Mutation rate = Total snvs/(2 * genome size * total generations) 
                             = 16177/(2 * 2100000000 * 233)  
	          = 1.65* 10 ^-8/ site/generation 















