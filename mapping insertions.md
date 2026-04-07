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
Now convert the primary alignment to bed,
If you want a BED of those mapped loci
awk 'BEGIN{FS=OFS="\t"} $0 !~ /^#/ {print $1, $4-1, $5, $3, ".", $7}' ../../edta_after_curation/EDTA_scaffold_lib_B73D10/ragtag.scaffold.fasta.mod.EDTA.TEanno.gff3 > b7310.edta.bed

# Fix insertion BED
``
awk 'BEGIN{FS=OFS="\t"} {sub(/_RagTag(_RagTag)?$/, "", $1); print}' taret_insertions.bed > b73d1_insertion.fixed.bed
``
# Fix EDTA BED
``
awk 'BEGIN{FS=OFS="\t"} {
    sub(/_RagTag.*/, "", $1)
    print
}' b73all.edta.bed > EDTA.fixed_b73all.bed
``
Then overlap with EDTA

bedtools intersect -wa -wb \
  -a B73D7B.insertion_loci_on_own_reference.bed \
  -b EDTA.bed > B73D7B.insertion_loci_vs_EDTA.tsv

If EDTA is GFF3 first convert it:
