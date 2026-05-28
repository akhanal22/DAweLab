I begin with the merge vcf for all individuals  (except b73all and b73d7b, although the file says nob73d7b only it actually doesnot have b73all too) but only private or unique vcf for 
each individual are present. The following module were loaded for this analysis
module load BCFtools/1.23.1-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0

1) First lets creat bed file for genome gfff file
```
awk 'BEGIN{OFS="\t"} $3=="gene"{
    print $1,$4-1,$5,$9
}' ../../J702_v3_modified_ragtag_re/J702.liftoff.gff3 > ./j702_genes.bed
```
2) lets make directory to store our information
```
mkdir -p per_Sample_genic_intergenic
```
3)Now, extract the information for the count of snps that fall in genic vs intergenic region for which i use the for loop
```
for sample in $(bcftools query -l private_snvs.only_one_linenob73d7b.vcf); do
    echo "Processing $sample"

    bcftools view -s "$sample" -m2 -M2 -v snps private_snvs.only_one_linenob73d7b.vcf \
    | bcftools view -i 'GT!="0/0" && GT!="./."' \
    | bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' \
    > per_sample_genic_intergenic/${sample}.private_snps.bed

    bedtools intersect \
        -a per_sample_genic_intergenic/${sample}.private_snps.bed \
        -b j702_genes.bed \
        -u \
    > per_sample_genic_intergenic/${sample}.genic_private_snps.bed

    bedtools intersect \
        -a per_sample_genic_intergenic/${sample}.private_snps.bed \
        -b j702_genes.bed \
        -v \
    > per_sample_genic_intergenic/${sample}.intergenic_private_snps.bed

    total=$(wc -l < per_sample_genic_intergenic/${sample}.private_snps.bed)
    genic=$(wc -l < per_sample_genic_intergenic/${sample}.genic_private_snps.bed)
    intergenic=$(wc -l < per_sample_genic_intergenic/${sample}.intergenic_private_snps.bed)
    ```
    

    echo -e "$sample\t$total\t$genic\t$intergenic"

done > per_sample_genic_intergenic/private_snp_genic_intergenic_counts.tsv
