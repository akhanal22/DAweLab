1. The pyhton file dotplot_TE_fixed.py contains the code to get dotplot for soecific region of genome and overap it with the TE annotations. Run it as follows:
python dotplot_TE_auto_xy.py \
  --paf chr4_RagTag_RagTag_1kb.paf \
  --ref_gff J702_raw_intact.renamed.gff3 \
  --qry_gff B73all_raw_intact.renamed.gff3 \
  --ref_chr chr4_RagTag_RagTag \
  --qry_chr chr4_RagTag_RagTag \
  --ref_name J702 \
  --qry_name B73all \
  --start 2200000 \
  --end 2400000 \
  --output b73all_chr4.png
2. snv_calling.md has info about how snvs are called from longcallD
3. 
