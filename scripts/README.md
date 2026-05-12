The pyhton file dotplot_TE_fixed.py contains the code to get dotplot for soecific region of genome and overap it with the TE annotations. Run it as follows:
python dotplot_TE_fixed.py \
  --paf b73all_chr1.paf \
  --ref_gff J702_raw_intact.renamed.gff3 \
  --qry_gff B73all_raw_intact.renamed.gff3 \
  --start 25000000 \
  --end 26000000 \
  --output chr1_25_26Mb_dotplot_TE_fixed.png
