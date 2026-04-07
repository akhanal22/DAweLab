This is for intersecting vcf insertions with reference to get flanking region 
``
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' B73D7B_ragtag.vcf | awk 'BEGIN{OFS="\t"} $4=="INS" {s=$2-2001; if(s<0)s=0; e=$2+2000; print $1,s,e,$3,$5}' > b73d7Binsertion_flanks_2kb.bed
``
After this I ran the minimap to get paf file with primary and secondary alignment

```
#!/bin/bash
#SBATCH --job-name=array_minimap
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL


#minimap2 -x asm5 -c --cs=long ../../B73D1_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d1.fa > insertion_flanks_vs_species_b73d1.paf
#minimap2 -x asm5 -c --cs=long ../../B73D2_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d2.fa > insertion_flanks_vs_species_b73d2.paf
#minimap2 -x asm5 -c --cs=long ../../B73D4_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d4.fa > insertion_flanks_vs_species_b73d4.paf
#minimap2 -x asm5 -c --cs=long ../../B73D7B_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d7b.fa > insertion_flanks_vs_species_b73d7b.paf
minimap2 -x asm5 -c --cs=long ../../B73D7c1_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d7c.fa > insertion_flanks_vs_species_b73d7c.paf
#minimap2 -x asm5 -c --cs=long ../../B73D10_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73d10.fa > insertion_flanks_vs_species_b73d10.paf
#minimap2 -x asm5 -c --cs=long ../../B73all_v3_modified_ragtag_re/ragtag.scaffold.fasta insertion_flanks_2kb_forb73all.fa > insertion_flanks_vs_species_b73all.paf

```
This gives both primary and secondary alignment, lets filter for primary only
```
grep "tp:A:P" insertion_flanks_vs_species_b73all.paf > insertion_flanks_vs_species_b73all_primary.paf
```
Lets find the exact coordinate of insertions from the paf file
```
#!/usr/bin/env python3

import sys
import re

def parse_cigar(cigar):
    """
    Parse a minimap2 cg:Z: CIGAR string into list of (length, op).
    Example: 2000M124D2001M -> [(2000,'M'), (124,'D'), (2001,'M')]
    """
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

def get_cg_tag(fields):
    """
    Extract cg:Z: tag from optional PAF fields.
    """
    for f in fields[12:]:
        if f.startswith("cg:Z:"):
            return f[5:]
    return None

def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: python paf_target_insertions_to_bed.py input.paf output.bed\n"
        )
        sys.exit(1)

    paf_file = sys.argv[1]
    bed_file = sys.argv[2]

    with open(paf_file) as fin, open(bed_file, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                continue

            query_name = fields[0]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            strand = fields[4]
            target_name = fields[5]
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            mapq = fields[11]

            cigar = get_cg_tag(fields)
            if cigar is None:
                continue

            ops = parse_cigar(cigar)

            # Current positions on query and target
            qpos = query_start
            tpos = target_start

            event_num = 0

            for length, op in ops:
                if op in ("M", "=", "X"):
                    qpos += length
                    tpos += length

                elif op == "I":
                    # insertion in query relative to target
                    # advances query only
                    qpos += length

                elif op == "D":
                    # deletion in query relative to target
                    # == sequence present in target/reference but absent in query
                    # This is the "inserted sequence on target"
                    bed_start = tpos
                    bed_end = tpos + length
                    event_num += 1

                    name = f"{query_name}|target_ins_{event_num}|len={length}|mapq={mapq}"
                    fout.write(
                        f"{target_name}\t{bed_start}\t{bed_end}\t{name}\t{length}\t{strand}\n"
                    )

                    tpos += length

                elif op == "N":
                    # skipped region on target
                    tpos += length

                elif op == "S":
                    # soft clip: query consumes
                    qpos += length

                elif op == "H":
                    # hard clip: consumes neither reference nor query sequence in alignment coordinates
                    pass

                elif op == "P":
                    # padding
                    pass

                else:
                    sys.stderr.write(f"Warning: unhandled CIGAR op {op} in line:\n{line}\n")

if __name__ == "__main__":
    main()
```
now lets run this code
```
python paf_target_insertions_to_bed.py your.paf target_insertions.bed
```
# Fix insertion BED
``
awk 'BEGIN{FS=OFS="\t"} {sub(/_RagTag(_RagTag)?$/, "", $1); print}' taret_insertions.bed > b73d1_insertion.fixed.bed
``
# Fix EDTA BED by renaming chromosomes
``
awk 'BEGIN{FS=OFS="\t"} {
    sub(/_RagTag.*/, "", $1)
    print
}' b73all.edta.bed > EDTA.fixed_b73all.bed
``
Then overlap with EDTA
```
bedtools intersect -wa -wb \
  -a b73d1_insertion.fixed.bed \
  -b EDTA.fixed_b73d1.bed > B73D1.insertion_loci_vs_EDTA.tsv
```
I ran into an issue here, i do get mutiple hits on TEs but it is from teanno.gff3 which has multiple fragmented TEs, so this time I am using raw.gff3 file from this edta_after_curation/EDTA_scaffold_lib_B73D1/ragtag.scaffold.fasta.mod.EDTA.raw/ragtag.scaffold.fasta.mod.EDTA.intact.raw.gff3. This has TEs for which there is atleast almost complete copy instead of very fragmented.
First, i use the following awk code to generated collapsed gff3 which is specifically true for ltr because i want only major element of row and other explain in column
```
BEGIN {
    FS = OFS = "\t"

    print "seqid","element_id","feature_type","classification","strand","start","end", \
          "region_id","region_start","region_end", \
          "lTSD_start","lTSD_end","rTSD_start","rTSD_end", \
          "lLTR_start","lLTR_end","rLTR_start","rLTR_end", \
          "internal_start","internal_end", \
          "TSD","TIR","motif","tsd_seq","ltr_identity","identity","method","Name"
}

$0 !~ /^#/ {

    seqid   = $1
    feature = $3
    start   = $4
    end     = $5
    strand  = $7
    attr    = $9

    id="."; parent="."; name="."; class="."; method="."; motif="."; tsd_seq="."
    tir="."; tsd="."; ltr_identity="."; identity="."

    n = split(attr, a, ";")
    for (i=1; i<=n; i++) {
        if (a[i] ~ /^ID=/)              id = substr(a[i],4)
        else if (a[i] ~ /^Parent=/)     parent = substr(a[i],8)
        else if (a[i] ~ /^Name=/)       name = substr(a[i],6)
        else if (a[i] ~ /^classification=/) class = substr(a[i],16)
        else if (a[i] ~ /^method=/)     method = substr(a[i],8)
        else if (a[i] ~ /^motif=/)      motif = substr(a[i],7)
        else if (a[i] ~ /^tsd=/)        tsd_seq = substr(a[i],5)
        else if (a[i] ~ /^TSD=/)        tsd = substr(a[i],5)
        else if (a[i] ~ /^TIR=/)        tir = substr(a[i],5)
        else if (a[i] ~ /^ltr_identity=/) ltr_identity = substr(a[i],14)
        else if (a[i] ~ /^identity=/)   identity = substr(a[i],10)
    }

    # ------------------------------
    # LTR repeat_region = top-level TE
    # ------------------------------
    if (feature == "repeat_region") {
        key = id

        seqidA[key]   = seqid
        featureA[key] = feature
        classA[key]   = class
        strandA[key]  = strand
        startA[key]   = start
        endA[key]     = end
        region_idA[key]    = id
        region_startA[key] = start
        region_endA[key]   = end
        methodA[key]  = method
        motifA[key]   = motif
        tsd_seqA[key] = tsd_seq
        ltr_idA[key]  = ltr_identity
        nameA[key]    = name
    }

    # ------------------------------
    # LTR subfeatures linked by Parent=repeat_region_X
    # ------------------------------
    else if (feature == "target_site_duplication") {
        key = parent
        if (!(key in seqidA)) {
            seqidA[key] = seqid
            featureA[key] = "repeat_region"
            classA[key] = class
            strandA[key] = strand
            methodA[key] = method
            motifA[key] = motif
            tsd_seqA[key] = tsd_seq
            ltr_idA[key] = ltr_identity
            nameA[key] = name
            region_idA[key] = key
        }

        if (id ~ /^lTSD_/) {
            lTSD_startA[key] = start
            lTSD_endA[key]   = end
        } else if (id ~ /^rTSD_/) {
            rTSD_startA[key] = start
            rTSD_endA[key]   = end
        }
    }

    else if (feature == "long_terminal_repeat") {
        key = parent
        if (!(key in seqidA)) {
            seqidA[key] = seqid
            featureA[key] = "repeat_region"
            classA[key] = class
            strandA[key] = strand
            methodA[key] = method
            motifA[key] = motif
            tsd_seqA[key] = tsd_seq
            ltr_idA[key] = ltr_identity
            nameA[key] = name
            region_idA[key] = key
        }

        if (id ~ /^lLTR_/) {
            lLTR_startA[key] = start
            lLTR_endA[key]   = end
        } else if (id ~ /^rLTR_/) {
            rLTR_startA[key] = start
            rLTR_endA[key]   = end
        }
    }

    else if (feature ~ /_LTR_retrotransposon$/) {
        key = parent
        if (!(key in seqidA)) {
            seqidA[key] = seqid
            featureA[key] = feature
            classA[key] = class
            strandA[key] = strand
            methodA[key] = method
            motifA[key] = motif
            tsd_seqA[key] = tsd_seq
            ltr_idA[key] = ltr_identity
            nameA[key] = name
            region_idA[key] = key
        }
        internal_startA[key] = start
        internal_endA[key]   = end
    }

    # ------------------------------
    # TIR elements as top-level one-row entries
    # ------------------------------
    else if (feature ~ /_TIR_transposon$/) {
        key = id

        seqidA[key]   = seqid
        featureA[key] = feature
        classA[key]   = class
        strandA[key]  = strand
        startA[key]   = start
        endA[key]     = end
        region_idA[key]    = id
        region_startA[key] = start
        region_endA[key]   = end
        tsdA[key]     = tsd
        tirA[key]     = tir
        identityA[key]= identity
        methodA[key]  = method
        nameA[key]    = name
    }

    # ------------------------------
    # Other TE rows: keep one row as top-level
    # ------------------------------
    else if (feature != "target_site_duplication" && feature != "long_terminal_repeat") {
        key = id

        seqidA[key]   = seqid
        featureA[key] = feature
        classA[key]   = class
        strandA[key]  = strand
        startA[key]   = start
        endA[key]     = end
        region_idA[key]    = id
        region_startA[key] = start
        region_endA[key]   = end
        methodA[key]  = method
        motifA[key]   = motif
        tsd_seqA[key] = tsd_seq
        tsdA[key]     = tsd
        tirA[key]     = tir
        ltr_idA[key]  = ltr_identity
        identityA[key]= identity
        nameA[key]    = name
    }
}

END {
    for (key in seqidA) {
        s  = (key in startA ? startA[key] : region_startA[key])
        e  = (key in endA ? endA[key] : region_endA[key])

        print seqidA[key], \
              key, \
              featureA[key], \
              classA[key], \
              strandA[key], \
              s, \
              e, \
              region_idA[key], \
              region_startA[key], \
              region_endA[key], \
              lTSD_startA[key], \
              lTSD_endA[key], \
              rTSD_startA[key], \
              rTSD_endA[key], \
              lLTR_startA[key], \
              lLTR_endA[key], \
              rLTR_startA[key], \
              rLTR_endA[key], \
              internal_startA[key], \
              internal_endA[key], \
              tsdA[key], \
              tirA[key], \
              motifA[key], \
              tsd_seqA[key], \
              ltr_idA[key], \
              identityA[key], \
              methodA[key], \
              nameA[key]
    }
}
```
Then run this code to get the tsv for all eight individual.  
```
awk -f collapse_edta_gff3.awk ../../edta_after_curation/EDTA_scaffold_lib_B73D1/ragtag.scaffold.fasta.mod.EDTA.raw/ragtag.scaffold.fasta.mod.EDTA.intact.raw.gff3 > active_ins/individual1_TE_collapsedb73d1.tsv
```
convert into bed format
```
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1, $6-1, $7, $2 "|" $4 "|" $3, ".", $5}' active_ins/individual1_TE_collapsedb73d1.tsv > active_ins/b73d1_TE_collapsed.bed
```
rename the chromosomes
```
awk 'BEGIN{FS=OFS="\t"} {
    sub(/_RagTag.*/, "", $1)
    print
}'  active_ins/b73d1_TE_collapsed.bed > active_ins/EDTA.fixed_b73d1.bed
```
redo the intersect for this step again, 
```
bedtools intersect -wa -wb   -a b73d1_insertion.fixed.bed   -b active_ins/EDTA.fixed_b73d1.bed > active_ins/B73D1.insertion_loci_vs_EDTA.tsv
```
