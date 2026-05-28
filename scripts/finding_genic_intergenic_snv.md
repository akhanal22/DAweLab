I begin with the merge vcf for all individuals  (except b73all and b73d7b, although the file says nob73d7b only it actually doesnot have b73all too) but only private or unique vcf for 
each individual are present.This vcf was created from seven_lines.merged.vcf.gz file by selecting only private snvsthat i am interested in. The code can be found in "snv_calling" scripts in the same directory.  The following module were loaded for this analysis
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
echo -e "sample\tgt_type\ttotal\tgenic\tintergenic" \
> per_sample_genic_intergenic/private_snp_genic_intergenic_counts.tsv

for sample in $(bcftools query -l private_snvs.only_one_linenob73d7b.vcf); do
    echo "Processing $sample"

    # Make BED with genotype
    bcftools view -s "$sample" -m2 -M2 -v snps private_snvs.only_one_linenob73d7b.vcf \
    | bcftools view -i 'GT!="0/0" && GT!="./."' \
    | bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\t[%GT]\n' \
    > per_sample_genic_intergenic/${sample}.private_snps.withGT.bed

    # Separate heterozygous SNPs
    awk '$5=="0/1" || $5=="1/0" || $5=="0|1" || $5=="1|0"' \
    per_sample_genic_intergenic/${sample}.private_snps.withGT.bed \
    > per_sample_genic_intergenic/${sample}.het_private_snps.bed

    # Separate homozygous ALT SNPs
    awk '$5=="1/1" || $5=="1|1"' \
    per_sample_genic_intergenic/${sample}.private_snps.withGT.bed \
    > per_sample_genic_intergenic/${sample}.hom_private_snps.bed

    for type in het hom; do

        bedtools intersect \
            -a per_sample_genic_intergenic/${sample}.${type}_private_snps.bed \
            -b j702_genes.bed \
            -u \
        > per_sample_genic_intergenic/${sample}.${type}_genic_private_snps.bed

        bedtools intersect \
            -a per_sample_genic_intergenic/${sample}.${type}_private_snps.bed \
            -b j702_genes.bed \
            -v \
        > per_sample_genic_intergenic/${sample}.${type}_intergenic_private_snps.bed

        total=$(wc -l < per_sample_genic_intergenic/${sample}.${type}_private_snps.bed)
        genic=$(wc -l < per_sample_genic_intergenic/${sample}.${type}_genic_private_snps.bed)
        intergenic=$(wc -l < per_sample_genic_intergenic/${sample}.${type}_intergenic_private_snps.bed)

        echo -e "$sample\t$type\t$total\t$genic\t$intergenic" \
        >> per_sample_genic_intergenic/private_snp_genic_intergenic_counts.tsv

    done
done
    ```
    

    echo -e "$sample\t$total\t$genic\t$intergenic"

done > per_sample_genic_intergenic/private_snp_genic_intergenic_counts.tsv

Now to make the figure for snp distribution, First make one combined BED with SNP type:

awk -v s="$sample" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"het"}' \
per_sample_genic_intergenic/${sample}.het_private_snps.bed \
> ${sample}.private_snps.plot.bed

awk -v s="$sample" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"hom"}' \
per_sample_genic_intergenic/${sample}.hom_private_snps.bed \
>> ${sample}.private_snps.plot.bed

Then use this R code to plot SNPs across chromosomes:

library(ggplot2)
library(dplyr)

sample <- "B73D1.sorted.bam"

snps <- read.table(
  paste0(sample, ".private_snps.plot.bed"),
  header = FALSE,
  sep = "\t"
)

colnames(snps) <- c("chr", "start", "end", "id", "gt_type")

snps <- snps %>%
  mutate(
    pos_mb = end / 1e6,
    chr = factor(chr, levels = unique(chr))
  )

ggplot(snps, aes(x = pos_mb, y = chr, color = gt_type)) +
  geom_point(size = 1.6, alpha = 0.7) +
  theme_bw() +
  labs(
    x = "Position on chromosome (Mb)",
    y = "Chromosome",
    color = "Genotype",
    title = paste("Private SNP distribution:", sample)
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  paste0(sample, "_private_snps_chromosome_distribution.png"),
  width = 10,
  height = 6,
  dpi = 300
)
