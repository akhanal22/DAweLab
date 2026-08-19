# we use following awk and python code to extract TE bed file insertion coordinate and TE bed file
```
awk -f make_Edta_bed.awk ./B73D10_ragtag.scaffold.fasta.mod.EDTA.intact.gff3 > B73D10.intact_TE.all_classes.bed
python extract_insertion.py ./B73all_vs_J702_modified.sorted.paf 100 > B73all_vs_J702.insertions_100bp.bed
```
Identify the active te family. remember to remove the header for insertions_100bp.bed file as bed only recognize integers for second and third column. 
```
bedtools intersect -a B73D4.intact_TE.all_classes.bed -b B73D4_vs_J702.insertions_100bp.bed -wa -wb     > B73D4_active_TE_greater_than_100bp.tsv
```
Identify intact LTR repeat from all te bed file
```
awk 'BEGIN{OFS="\t"}
$5 ~ /^LTR\// && $9=="repeat_region" {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9
}' B73D1.intact_TE.all_classes.bed \
> intact_LTR_repeat_regions.bed
```
now to intersect deletion in B73D1 with respsect to J702 
```
bedtools intersect \
    -a deletions_gt100bp.bed \
    -b intact_LTR_repeat_regions.bed \
    -wa -wb \
    > deletions_in_intact_LTR.tsv
```
