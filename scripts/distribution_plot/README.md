. This folder contains python code for potiing whole genome indels and svs using deepvariant and syri. only validated svs are plotted. all the files needed are also here. You just need to 
load python3 and run this using python3.

  python plot_indel_and_validated_svs.py \
  --deepvariant "upload/private_biallelic_indels.all_sizes.exactly_one_line.chr1-10.vcf.gz" \
  --syri "upload/syri.final.CG_di_tri_then_AT.no_N3.vcf.gz" \
  --insertions "upload/insertion_validation_with_cigar_with500kbfilteredvcfagainstj702.tsv" \
  --deletions "upload/deletion_validation_withcigar_for500kbfilteredvcfagainstj702.tsv" \
  --gff "upload/J702_flip.fa.mod.EDTA.TEanno.gff3" \
  --depth-ratio 0.70 \
  --output-prefix "ma_indel_sv_distributions"

  
