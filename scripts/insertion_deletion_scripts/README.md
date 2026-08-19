This contains awk and python scripts i used for identifyting active tes based on genome to genome alignment
# we use following awk and python code to extract TE bed file insertion coordinate and TE bed file
```
awk -f make_Edta_bed.awk ./B73D10_ragtag.scaffold.fasta.mod.EDTA.intact.gff3 > B73D10.intact_TE.all_classes.bed
python extract_insertion.py ./B73all_vs_J702_modified.sorted.paf 100 > B73all_vs_J702.insertions_100bp.bed
```
Identify the active te family. remember to remove the header for insertions_100bp.bed file as bed only recognize integers for second and third column. 
```
bedtools intersect -a B73D4.intact_TE.all_classes.bed -b B73D4_vs_J702.insertions_100bp.bed -wa -wb     > B73D4_active_TE_greater_than_100bp.tsv

```
Get the deletions bed for deletions > 100bp

```
python extract_deletions.py B73D1_vs_J702_modified.sorted.paf 100     > B73D1_deletions_gt100bp.bed
```
Identify intact LTR repeat with structure from all te bed file
```
python make_LTR_structure_bed.py B73D1_ragtag.scaffold.fasta.mod.EDTA.intact.gff3 > B73D1.LTR_structure.tsv
```
now to classify the deletion 
```
python classify_ltr_deletion.py \
B73D1_deletions_gt100bp.bed \
B73D1.LTR_structure.tsv \
> B73D1_LTR_deletion_classification.tsv

```

