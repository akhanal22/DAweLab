This is for intersecting vcf insertions with reference to get flanking region 
``
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' B73D7B_ragtag.vcf | awk 'BEGIN{OFS="\t"} $4=="INS" {s=$2-2001; if(s<0)s=0; e=$2+2000; print $1,s,e,$3,$5}' > b73d7Binsertion_flanks_2kb.bed
``

So:

awk 'BEGIN{OFS="\t"} {print $1,$6,$8,$9}' insertion_flanks_vs_species_b73d7b.paf

This gives:

insertion/flank ID
species-reference chromosome
mapped start
mapped end

That mapped interval is the species-own-reference coordinate of the insertion locus.

If you want a BED of those mapped loci
awk 'BEGIN{OFS="\t"} {print $6,$8,$9,$1}' insertion_flanks_vs_species_b73d7b.paf > B73D7B.insertion_loci_on_own_reference.bed
Then overlap with EDTA

If your EDTA annotation is already in BED:

bedtools intersect -wa -wb \
  -a B73D7B.insertion_loci_on_own_reference.bed \
  -b EDTA.bed > B73D7B.insertion_loci_vs_EDTA.tsv

If EDTA is GFF3 first convert it:
