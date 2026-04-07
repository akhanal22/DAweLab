This is for intersecting vcf insertions with reference to get flanking region 
``
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' B73D7B_ragtag.vcf | awk 'BEGIN{OFS="\t"} $4=="INS" {s=$2-2001; if(s<0)s=0; e=$2+2000; print $1,s,e,$3,$5}' > b73d7Binsertion_flanks_2kb.bed
``
After this I ran the minimap to get paf file with primary and secondary alignment
``
#!/bin/bash
#SBATCH --job-name=array_minimap
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL


#minimap2 -x asm5 -c --cs=long ../../B73D1_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d1.fa > insertion_flanks_vs_species_b73d1.paf
#minimap2 -x asm5 -c --cs=long ../../B73D2_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d2.fa > insertion_flanks_vs_species_b73d2.paf
#minimap2 -x asm5 -c --cs=long ../../B73D4_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d4.fa > insertion_flanks_vs_species_b73d4.paf
#minimap2 -x asm5 -c --cs=long ../../B73D7B_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d7b.fa > insertion_flanks_vs_species_b73d7b.paf
minimap2 -x asm5 -c --cs=long ../../B73D7c1_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d7c.fa > insertion_flanks_vs_species_b73d7c.paf
#minimap2 -x asm5 -c --cs=long ../../B73D10_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d10.fa > insertion_flanks_vs_species_b73d10.paf
#minimap2 -x asm5 -c --cs=long ../../B73all_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73all.fa > insertion_flanks_vs_species_b73all.paf
``
This gives both primary and secondary alignment, lets filter for primary only
``
grep "tp:A:P" insertion_flanks_vs_species_b73all.paf > insertion_flanks_vs_species_b73all_primary.paf
``
Now convert the primary alignment to bed,
If you want a BED of those mapped loci
awk 'BEGIN{FS=OFS="\t"} $0 !~ /^#/ {print $1, $4-1, $5, $3, ".", $7}' ../../edta_after_curation/EDTA_scaffold_lib_B73D10/ragtag.scaffold.fasta.mod.EDTA.TEanno.gff3 > b7310.edta.bed

# Fix insertion BED
``
awk 'BEGIN{FS=OFS="\t"} {sub(/_RagTag(_RagTag)?$/, "", $1); print}' B73d1.insertion_loci_on_own_reference.bed > b73d1_insertion.fixed.bed
``
# Fix EDTA BED
``
awk 'BEGIN{FS=OFS="\t"} {
    sub(/_RagTag.*/, "", $1)
    print
}' b73all.edta.bed > EDTA.fixed_b73all.bed
``
Then overlap with EDTA

bedtools intersect -wa -wb \
  -a B73D7B.insertion_loci_on_own_reference.bed \
  -b EDTA.bed > B73D7B.insertion_loci_vs_EDTA.tsv

If EDTA is GFF3 first convert it:
