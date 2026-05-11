with just whole genome plotsr, everything looks aligned, however there are minor changes contributing to differences in chromosome size which i am investigating. I am focusing on 251-252mb
of chromosome 1 first using following awk command
```
awk -F'\t' '$1=="chr1_RagTag_RagTag" && $6=="chr1_RagTag_RagTag" && (($2<=252000000 && $3>=251000000) || ($7<=252000000 && $8>=251000000))' B73D1syri.out > chr1_251_252Mb.variants.out
```
Now lets plot it
```
plotsr   --sr chr1_251_252Mb.variants.out   --genomes genomes.txt   -o B73D1_chr1_251_252Mb_HDR6.png
```
