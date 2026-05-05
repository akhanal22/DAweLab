We start wit h the vcf generated from sniffles by mapping raw reads agaons J702 reference genome. First, lets remove the variant that is also found in J702.

```
 bcftools isec -p B73D1 -n=1 B73D1_ragtag.vcf J702_ragtag.vcf
```
We do this for all individuals and the output is in B73D1. Within this folder is 0000.vcf  0001.vcf  0002.vcf  0003.vcf  README.txt files and the first one i.e. 0000.vcf is the one that
has vcf for B73D1 which excludes the common variant with J702. To get information about the count of different svs , i use following code:
```
bcftools query -f '%INFO/SVTYPE\n' B73D1/0000.vcf | sort | uniq -c | sort -k1,1nr
```
Now i want different sizes of insertion, for which i use bcftools again
```
bcftools view -i 'SVTYPE="INS"' B73D1/0000.vcf | bcftools query -f '%INFO/SVLEN\n' > ins_sizes.txt
```
