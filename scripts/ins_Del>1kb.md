```
awk 'BEGIN{
OFS="\t";
print "ref_chr","ref_start","ref_end","size_reference_bp","query_chr","query_start","query_end","size_target_bp","variant_id","type"
}

($11=="INS" || $11=="DEL") {

    size_reference = $3 - $2 + 1
    size_target    = $8 - $7 + 1

    # SV size difference
    sv_size = size_target - size_reference
    if(sv_size < 0) sv_size = -sv_size

    # keep only SVs > 1 kb
    if(sv_size >= 1000){

        print $1,$2,$3,size_reference,
              $6,$7,$8,size_target,
              $9,$11
    }
}' B73D1syri.out > B73D1_INS_DEL_gt1kb_coordinates_sizes.tsv
```
