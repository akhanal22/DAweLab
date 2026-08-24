
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
   bcftools filter \
  -e 'FILTER!="PASS" || FORMAT/DP<40 || FORMAT/DP>100 || FORMAT/GQ<20' \
  B73D1.longcallD.snps.vcf \
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
  -o B73D1.longcallD.snps.DP40-100.GQ20.rigorous.vcf.gz

bcftools index -t B73D1.longcallD.snps.DP40-100.GQ20.rigorous.vcf.gz

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
    if(a[1]=="0/1" || a[1]=="1/0" || a[1]=="1/1" || a[1]=="0|1" || a[1]=="1|0" || a[1]=="1|1"){
      count[sample[i]]++
    }
  }
}
END{
  for(s in count) print s, count[s]
}' private_snvs.only_one_line.vcf > private_snv_counts.tsv
```
However, since, B73D4 and B73all were same for 33 generations and B73D7b and B73D7c were same until 15 generation , I found the common snps across them using following code
```
bcftools view -m2 -M2 -v snps seven_lines.merged.vcf.gz \
| awk '
BEGIN{OFS="\t"}

function is_alt(gt_field){
  split(gt_field,a,":")
  gt=a[1]
  return (gt=="0/1" || gt=="1/0" || gt=="1/1" || gt=="0|1" || gt=="1|0" || gt=="1|1")
}

/^#CHROM/{
  for(i=10;i<=NF;i++){
    sample[$i]=i
  }
  next
}

!/^#/{
  d4     = is_alt($sample["B73D4.sorted.bam"])
  all    = is_alt($sample["B73all.sorted.bam"])
  d7b    = is_alt($sample["B73D7B.sorted.bam"])
  d7c    = is_alt($sample["B73D7C.sorted.bam"])

  if(d4 && all){
    count["B73D4_B73all_shared"]++
  }

  if(d7b && d7c){
    count["B73D7B_B73D7C_shared"]++
  }
}

END{
  print "Category","Count"
  print "B73D4_B73all_shared", count["B73D4_B73all_shared"]+0
  print "B73D7B_B73D7C_shared", count["B73D7B_B73D7C_shared"]+0
}
' > shared_snv_counts.tsv
```
Based on this following table represents private snv counts for each (except B73D4) and shared snvs between B73D4 and B73all and B7D7B and B73D7C

| Line | Private SNVs | Generations |
|------|---------------|-------------|
| B73D1 | 3520 | 33 |
| B73D2 | 3580 | 33 |
| B73D4+B73all | 3878 | 33 |
| B73D7B | 2142 | 18 |
| B73D7c | 2289 | 18 |
| B73D10 | 3478 | 33 |
| B73all | 770 | 2 |
| B73D7B+B73D7C | 1787 | 15 | 
| **Total** | **21444** | **185** |

Now we calculate the mutation rate as below:

For estimating mutation  rate:
Total unique SNVs (sum all) = 21444
Total generation = 185
Genome size = 2.1gb
Mutation rate = Total snvs/(2 * genome size * total generations)  
	          = 2.76* 10 ^-8/ site/generation 



Following the approach from arabidopsis paper, lets count private heterozygous and homozygous mutation first:
```
bcftools view -m2 -M2 -v snps seven_lines.merged.vcf.gz \
| awk '
BEGIN{OFS="\t"}

function gt_type(gt_field){
  split(gt_field,a,":")
  gt=a[1]

  if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") return "het"
  if(gt=="1/1" || gt=="1|1") return "hom"

  return "ref"
}

/^#CHROM/{
  for(i=10;i<=NF;i++){
    sample[$i]=i
  }
  next
}

!/^#/{
  d4  = gt_type($sample["B73D4.sorted.bam"])
  all = gt_type($sample["B73all.sorted.bam"])
  d7b = gt_type($sample["B73D7B.sorted.bam"])
  d7c = gt_type($sample["B73D7C.sorted.bam"])

  if(d4!="ref" && all!="ref"){
    count["B73D4_B73all_shared_total"]++

    if(d4=="het" && all=="het") count["B73D4_B73all_shared_het_het"]++
    else if(d4=="hom" && all=="hom") count["B73D4_B73all_shared_hom_hom"]++
    else count["B73D4_B73all_shared_mixed_het_hom"]++
  }

  if(d7b!="ref" && d7c!="ref"){
    count["B73D7B_B73D7C_shared_total"]++

    if(d7b=="het" && d7c=="het") count["B73D7B_B73D7C_shared_het_het"]++
    else if(d7b=="hom" && d7c=="hom") count["B73D7B_B73D7C_shared_hom_hom"]++
    else count["B73D7B_B73D7C_shared_mixed_het_hom"]++
  }
}

END{
  print "Category","Count"

  print "B73D4_B73all_shared_total", count["B73D4_B73all_shared_total"]+0
  print "B73D4_B73all_shared_het_het", count["B73D4_B73all_shared_het_het"]+0
  print "B73D4_B73all_shared_hom_hom", count["B73D4_B73all_shared_hom_hom"]+0
  print "B73D4_B73all_shared_mixed_het_hom", count["B73D4_B73all_shared_mixed_het_hom"]+0

  print "B73D7B_B73D7C_shared_total", count["B73D7B_B73D7C_shared_total"]+0
  print "B73D7B_B73D7C_shared_het_het", count["B73D7B_B73D7C_shared_het_het"]+0
  print "B73D7B_B73D7C_shared_hom_hom", count["B73D7B_B73D7C_shared_hom_hom"]+0
  print "B73D7B_B73D7C_shared_mixed_het_hom", count["B73D7B_B73D7C_shared_mixed_het_hom"]+0
}
' > shared_snv_counts_het_hom.tsv
```
Now shared count
```
bcftools view -m2 -M2 -v snps seven_lines.merged.vcf.gz \
| awk '
BEGIN{OFS="\t"}

function gt_type(gt_field){
  split(gt_field,a,":")
  gt=a[1]

  if(gt=="0/1" || gt=="1/0" || gt=="0|1" || gt=="1|0") return "het"
  if(gt=="1/1" || gt=="1|1") return "hom"

  return "ref"
}

/^#CHROM/{
  for(i=10;i<=NF;i++){
    sample[$i]=i
  }
  next
}

!/^#/{
  d4  = gt_type($sample["B73D4.sorted.bam"])
  all = gt_type($sample["B73all.sorted.bam"])
  d7b = gt_type($sample["B73D7B.sorted.bam"])
  d7c = gt_type($sample["B73D7C.sorted.bam"])

  if(d4!="ref" && all!="ref"){
    count["B73D4_B73all_shared_total"]++

    if(d4=="het" && all=="het") count["B73D4_B73all_shared_het_het"]++
    else if(d4=="hom" && all=="hom") count["B73D4_B73all_shared_hom_hom"]++
    else count["B73D4_B73all_shared_mixed_het_hom"]++
  }

  if(d7b!="ref" && d7c!="ref"){
    count["B73D7B_B73D7C_shared_total"]++

    if(d7b=="het" && d7c=="het") count["B73D7B_B73D7C_shared_het_het"]++
    else if(d7b=="hom" && d7c=="hom") count["B73D7B_B73D7C_shared_hom_hom"]++
    else count["B73D7B_B73D7C_shared_mixed_het_hom"]++
  }
}

END{
  print "Category","Count"

  print "B73D4_B73all_shared_total", count["B73D4_B73all_shared_total"]+0
  print "B73D4_B73all_shared_het_het", count["B73D4_B73all_shared_het_het"]+0
  print "B73D4_B73all_shared_hom_hom", count["B73D4_B73all_shared_hom_hom"]+0
  print "B73D4_B73all_shared_mixed_het_hom", count["B73D4_B73all_shared_mixed_het_hom"]+0

  print "B73D7B_B73D7C_shared_total", count["B73D7B_B73D7C_shared_total"]+0
  print "B73D7B_B73D7C_shared_het_het", count["B73D7B_B73D7C_shared_het_het"]+0
  print "B73D7B_B73D7C_shared_hom_hom", count["B73D7B_B73D7C_shared_hom_hom"]+0
  print "B73D7B_B73D7C_shared_mixed_het_hom", count["B73D7B_B73D7C_shared_mixed_het_hom"]+0
}
' > shared_snv_counts_het_hom.tsv
```


okAy now to estimate genome size based on deepvariant
```
#!/bin/bash

MIN_GQ=30
MIN_DP=40

regions=$(printf "chr%s_RagTag," {1..10} | sed 's/,$//')

for gvcf in \
B73D1.deepvariant.g.vcf.gz \
B73D2.deepvariant.g.vcf.gz \
B73D4.deepvariant.g.vcf.gz \
B73D10.deepvariant.g.vcf.gz \
B73D7c.deepvariant.g.vcf.gz
do

    sample=$(basename "$gvcf" .deepvariant.g.vcf.gz)

    echo "Processing $sample ..." >&2

    bcftools query \
        -r "$regions" \
        -f '%CHROM\t%POS\t%INFO/END[\t%GT\t%GQ\t%MIN_DP\t%DP]\n' \
        "$gvcf" |
    awk -v sample="$sample" \
        -v mingq="$MIN_GQ" \
        -v mindp="$MIN_DP" '
    {
        start=$2
        end=$3

        # If END is absent, this record represents one base
        if(end=="." || end=="")
            end=start

        gt=$4
        gq=$5
        min_dp=$6
        dp=$7

        # DeepVariant reference blocks generally use MIN_DP.
        # Variant records generally use DP.
        if(min_dp!="." && min_dp!="")
            depth=min_dp
        else
            depth=dp

        # Skip missing genotypes
        if(gt=="./." || gt==".|." || gt==".")
            next

        # GQ < 30 fails
        if(gq=="." || gq < mingq)
            next

        # DP <= 40 fails
        if(depth=="." || depth <= mindp)
            next

        # Count all bases represented by this passing gVCF record
        callable += end-start+1
    }

    END {
        print sample "\t" callable
    }'

done
```
```
bcftools view \
    -f PASS \
    -v snps \
    -i 'QUAL>=50' \
    five_lines.merged.vcf.gz \
    -Oz \
    -o five_lines.merged.pass.snps.vcf.gz

```
get private snps only. each vcf were already filtered for dp<40 GQ<30 and qual>50 and normalized and variant shared with J702 were removed.
```
bcftools view \
    -m2 -M2 \
    -v snps \
    -i 'COUNT(GT="alt")=1' \
    five_lines.merged.filtered.vcf.gz \
    -Oz \
    -o private_snvs.only_one_line.vcf.gz

bcftools index -c private_snvs.only_one_line.vcf.gz
```
count private snps
```
bcftools query \
    -f '[%SAMPLE\t%GT\n]' \
    private_snvs.only_one_line.vcf.gz |
awk '
BEGIN{OFS="\t"}
$2=="0/1" || $2=="1/0" || $2=="1/1" ||
$2=="0|1" || $2=="1|0" || $2=="1|1" {
    count[$1]++
}
END{
    print "Sample","Private_SNVs"
    for(s in count)
        print s,count[s]
}' > private_snv_counts.tsv
```
