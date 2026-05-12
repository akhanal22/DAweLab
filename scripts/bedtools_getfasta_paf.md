When creating the dotplot for whole genome and whole chromosomes, it doesnot show the small insertions and deletion. So, we need to separate the fasta suence and here we separate it in
1kb bins. To do that, first I create the bed file for each chromosome of each individual.
```
hile read chr
do
    echo "Processing $chr"

    len=$(grep -w "$chr" ../../B73D2_v3_modified_ragtag_re_v2/chr1_10_only.fa.fai | cut -f2)

    echo -e "${chr}\t${len}" > ${chr}.genome

    bedtools makewindows \
        -g ${chr}.genome \
        -w 1000 \
        > B73D2/${chr}_1kb_b73d2.bed

done < ./chr1_10.list
```
Do this for each individual
