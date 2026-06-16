I do it for each sample individually, you can also create a loop and save the same name in a txt file.
```
SAMPLE=B73D1
OUTDIR=sv_TE_gap_filter_all/$SAMPLE
CHROMSIZE=chromsize.tsv
```
MAke SV delta table: J702-B73D1
```
awk 'BEGIN{OFS="\t"}
FNR==NR{
    sizeB[$4]=$6
    next
}
{
    chrA=$1
    id=$4
    type=$5
    sizeA=$6

    if(id in sizeB){
        delta=sizeA-sizeB[id]
        print chrA,id,type,sizeA,sizeB[id],delta
    }
}' \
$OUTDIR/${SAMPLE}.all_nonSNP_SVs.genomeB.bed \
$OUTDIR/${SAMPLE}.all_nonSNP_SVs.genomeA.bed \
> $OUTDIR/${SAMPLE}.SV_delta_J702_minus_B.tsv
```
Identify TE associated SV
```
cut -f4 $OUTDIR/${SAMPLE}.SVs_noGap_intersect_EDTA.tsv | sort -u > $OUTDIR/${SAMPLE}.TE.ids

awk 'NR==FNR{keep[$1]=1; next}
$2 in keep {bp[$1]+=$6}
END{for(chr in bp) print chr,bp[chr]}' OFS="\t" \
$OUTDIR/${SAMPLE}.TE.ids \
$OUTDIR/${SAMPLE}.SV_delta_J702_minus_B.tsv \
> $OUTDIR/${SAMPLE}.TE_delta_per_chr.tsv
```
IDentify non TE associated SV per chromosome
```
cut -f4 $OUTDIR/${SAMPLE}.SVs_noGap_not_intersect_EDTA.tsv | sort -u > $OUTDIR/${SAMPLE}.nonTE.ids

awk 'NR==FNR{keep[$1]=1; next}
$2 in keep {bp[$1]+=$6}
END{for(chr in bp) print chr,bp[chr]}' OFS="\t" \
$OUTDIR/${SAMPLE}.nonTE.ids \
$OUTDIR/${SAMPLE}.SV_delta_J702_minus_B.tsv \
> $OUTDIR/${SAMPLE}.nonTE_delta_per_chr.tsv
```
Create tsv for not aligned bp of j702 per chromosome. Keep in mind that these are for region in j702 not aligned in b73D1, because we mapped B73D1 against j702. It is not region of B73D1
not aligned against J702.
```
awk '$3 ~ /^NOTAL/ {
    bp[$1]+=$3-$2
}
END{for(chr in bp) print chr,bp[chr]}' OFS="\t" \
$OUTDIR/${SAMPLE}.all_nonSNP_SVs.genomeA.bed \
> $OUTDIR/${SAMPLE}.NOTAL_bp_per_chr.tsv
```
create final summary table. Here chromsize_B73D1.tsv is the table that has chromosome name for both genome and size difference for each chromosome (J702-B73D1). Example of one is also uploaded here.
```
awk 'BEGIN{OFS="\t"}
FNR==1 && ARGIND==1 {next}

ARGIND==1 {
    chr=$1
    gsub(/^Chr/,"chr",chr)
    obs[chr]=$5
}
ARGIND==2 {notal[$1]=$2}
ARGIND==3 {te[$1]=$2}
ARGIND==4 {nonte[$1]=$2}

END{
    print "Chr","Observed_J702_minus_B73D1","NOTAL_bp","TE_SV_delta","nonTE_SV_delta","Explained_delta","Residual"

    for(chr in obs){
        n=notal[chr]+0
        t=te[chr]+0
        nt=nonte[chr]+0
        explained=n+t+nt
        residual=obs[chr]-explained

        print chr,obs[chr],n,t,nt,explained,residual
    }
}' \
chromsize_B73D1.tsv \
$OUTDIR/${SAMPLE}.NOTAL_bp_per_chr.tsv \
$OUTDIR/${SAMPLE}.TE_delta_per_chr.tsv \
$OUTDIR/${SAMPLE}.nonTE_delta_per_chr.tsv \
> $OUTDIR/${SAMPLE}.chromosome_size_accounting.tsv
```
Now I realized, i didnot account for size difference in gap region, so i am adding that
```
awk 'BEGIN{OFS="\t"}
{
    chrA=$8
    id=$4
    type=$5
    sizeB=$6
    sizeA=$11
    delta=sizeA-sizeB

    bp[chrA]+=delta
}
END{
    for(chr in bp)
        print chr,bp[chr]
}' $OUTDIR/${SAMPLE}.SVs_inside_AGP_gaps.genomeB_coordinates.tsv \
> $OUTDIR/${SAMPLE}.AGP_gap_SV_delta_per_chr.tsv

```
Adding this info in final table
```
awk 'BEGIN{OFS="\t"}
FNR==1 && ARGIND==1 {next}

ARGIND==1 {
    chr=$1
    gsub(/^Chr/,"chr",chr)
    obs[chr]=$5
}
ARGIND==2 {notal[$1]=$2}
ARGIND==3 {te[$1]=$2}
ARGIND==4 {nonte[$1]=$2}
ARGIND==5 {gap[$1]=$2}

END{
    print "Chr","Observed_J702_minus_B73D1","NOTAL_bp","AGP_gap_SV_delta","TE_SV_delta","nonTE_SV_delta","Explained_delta","Residual"

    for(chr in obs){
        n=notal[chr]+0
        g=gap[chr]+0
        t=te[chr]+0
        nt=nonte[chr]+0

        explained=n+g+t+nt
        residual=obs[chr]-explained

        print chr,obs[chr],n,g,t,nt,explained,residual
    }
}' \
$CHROMSIZE \
$OUTDIR/${SAMPLE}.NOTAL_bp_per_chr.tsv \
$OUTDIR/${SAMPLE}.TE_delta_per_chr.tsv \
$OUTDIR/${SAMPLE}.nonTE_delta_per_chr.tsv \
$OUTDIR/${SAMPLE}.AGP_gap_SV_delta_per_chr.tsv \
> $OUTDIR/${SAMPLE}.chromosome_size_accounting.with_AGPgapSVs.tsv
```
Okay i think not_al and residual makes it confusing so i will only have residual for all unaligned parts
```
awk 'BEGIN{OFS="\t"}
NR==1{
    print "Chr",
          "Observed_J702_minus_B73D1",
          "AGP_gap_SV_delta",
          "TE_SV_delta",
          "nonTE_SV_delta",
          "Explained_delta",
          "Residual_unexplained"
    next
}
{
    obs=$2
    agp=$4
    te=$5
    nonte=$6

    explained = agp + te + nonte
    residual = obs - explained

    print $1, obs, agp, te, nonte, explained, residual
}' B73D1.chromosome_size_accounting.with_AGPgapSVs.tsv \
> B73D1.chromosome_size_accounting.no_NOTAL.tsv
```
