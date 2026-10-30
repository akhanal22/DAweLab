#!/bin/bash
#SBATCH --job-name=callable
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# Run: sbatch estimate_callable_snps_indels.sh
# Run from the directory containing deepvariant_against_J702.
# Produces a confidence/depth-based reference-position mask, not an
# indel-size-specific sensitivity estimate. No repeat or 500-kb mask is applied.
# SNPs and indels each contribute their single reference anchor position.
# Confident reference calls and qualifying heterozygous calls remain included.
REFERENCE="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"
INPUT_ROOT="deepvariant_against_J702"
OUTPUT_ROOT="callable_snps_indels_re"
SAMPLES=(B73D1 B73D2 B73D4 B73D7c B73D10)
for command in python3 samtools bedtools; do
    command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 1; }
done
[[ -s "$REFERENCE" ]] || { echo "Missing reference: $REFERENCE" >&2; exit 1; }
mkdir -p "$OUTPUT_ROOT"
[[ -s "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
awk 'BEGIN{OFS="\t"} /^chr([1-9]|10)_RagTag\t/ {print $1,$2}' "${REFERENCE}.fai" > "$OUTPUT_ROOT/chr1-10.genome"
[[ $(wc -l < "$OUTPUT_ROOT/chr1-10.genome") -eq 10 ]] || { echo 'Expected ten RagTag chromosomes' >&2; exit 1; }

# Generate the ACGT mask once for all samples.
python3 - "$REFERENCE" "$OUTPUT_ROOT/reference_ACGT.bed" <<'PY'
import re, sys
chroms={f'chr{i}_RagTag' for i in range(1,11)}
name=None; parts=[]
def emit(out):
    if name in chroms:
        for m in re.finditer('[ACGTacgt]+',''.join(parts)):
            out.write(f'{name}\t{m.start()}\t{m.end()}\n')
with open(sys.argv[1]) as src, open(sys.argv[2],'w') as out:
    for line in src:
        if line.startswith('>'):
            emit(out); name=line[1:].split()[0]; parts=[]
        elif name in chroms:
            parts.append(line.strip())
    emit(out)
PY
bedtools sort -g "$OUTPUT_ROOT/chr1-10.genome" -i "$OUTPUT_ROOT/reference_ACGT.bed" > "$OUTPUT_ROOT/reference_ACGT.sorted.bed"
printf 'Sample\tCallable_bases\tBED\n' > "$OUTPUT_ROOT/summary_all_samples.tsv"
for sample in "${SAMPLES[@]}"; do
    gvcf="$INPUT_ROOT/gvcf/${sample}.chr1-10.g.vcf.gz"
    bam="$INPUT_ROOT/prepared_bams/${sample}.sorted.primary_mapq20.RG.bam"
    outdir="$OUTPUT_ROOT/$sample"
    for input in "$gvcf" "$bam"; do
        [[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 1; }
    done
    [[ -s "${bam}.bai" || -s "${bam%.bam}.bai" || -s "${bam}.csi" || -s "${bam%.bam}.csi" ]] || { echo "Missing index: $bam" >&2; exit 1; }
    mkdir -p "$outdir"
    echo "Processing $sample"
    python3 - "$gvcf" "$outdir/gvcf_confident.bed" "$outdir/gvcf_records_summary.tsv" <<'PY'
import collections, gzip, re, sys
source,destination,summary=sys.argv[1:]
chroms={f'chr{i}_RagTag' for i in range(1,11)}
counts=collections.Counter()
def num(x):
    try: return float(x)
    except (TypeError,ValueError): return -1
header=False
with gzip.open(source,'rt') as src, open(destination,'w') as out:
    for line in src:
        if line.startswith('#CHROM'):
            if len(line.rstrip().split('\t'))!=10:
                raise SystemExit('Expected single-sample gVCF')
            header=True; continue
        if line.startswith('#'): continue
        if not header: raise SystemExit('Missing #CHROM header')
        chrom,pos,_,ref,alt,qual,filt,info,fmt,sample=line.rstrip('\n').split('\t')[:10]
        if chrom not in chroms: continue
        pos=int(pos); v=dict(zip(fmt.split(':'),sample.split(':')))
        gt=v.get('GT','.').replace('|','/')
        if num(v.get('GQ'))<30: continue
        if alt in ('<*>','<NON_REF>') and gt=='0/0':
            depth=num(v.get('MIN_DP'))
            if depth<0: depth=num(v.get('DP'))
            if not 40<=depth<=100: continue
            tags=dict(x.split('=',1) for x in info.split(';') if '=' in x)
            end=int(tags.get('END',pos))
            if end<pos: raise SystemExit(f'Invalid END: {chrom}:{pos}')
            out.write(f'{chrom}\t{pos-1}\t{end}\n')
            counts['reference_blocks']+=1
            continue
        # Exactly one sequence ALT; placeholders do not count as real alleles.
        alleles=alt.split(',')
        real=[(i,a) for i,a in enumerate(alleles,1) if a not in ('<*>','<NON_REF>')]
        if len(real)!=1: continue
        alt_index,alternate=real[0]
        if not re.fullmatch('[ACGTacgt]+',ref) or not re.fullmatch('[ACGTacgt]+',alternate): continue
        if ref.upper()==alternate.upper(): continue
        if len(ref)==len(alternate) and len(ref)!=1: continue # exclude MNVs
        if not 40<=num(v.get('DP'))<=100: continue
        if gt=='0/0':
            category='reference_at_candidate'
        else:
            if filt!='PASS': continue
            try:
                genotype=[int(x) for x in gt.split('/')]
                ad=[int(x) for x in v['AD'].split(',')]
                rd,ad_alt=ad[0],ad[alt_index]
            except (ValueError,KeyError,IndexError): continue
            if len(genotype)!=2 or rd<0 or ad_alt<0 or rd+ad_alt==0: continue
            af=ad_alt/(rd+ad_alt)
            if genotype==[alt_index,alt_index] and af>=0.95:
                category='homozygous_alt'
            elif sorted(genotype)==[0,alt_index] and 0.40<=af<=0.60:
                category='heterozygous'
            else: continue
            category+= '_SNP' if len(ref)==len(alternate)==1 else '_indel'
        # Anchor only: do not declare every deleted base callable from one record.
        out.write(f'{chrom}\t{pos-1}\t{pos}\n')
        counts[category]+=1
with open(summary,'w') as out:
    out.write('Category\tRetained_records_before_depth_ACGT_intersection\n')
    for k in sorted(counts): out.write(f'{k}\t{counts[k]}\n')
PY
    # Input BAM must already contain the desired primary MAPQ>=20 reads.
    # These depth settings add no positive base/MAPQ cutoff.
    # Separate chromosome files ensure any samtools failure stops the job.
    : > "$outdir/depth_40_100.bed"
    for chromosome_number in {1..10}; do
        samtools depth -aa -q 0 -Q 0 -r "chr${chromosome_number}_RagTag" "$bam" |
        awk 'BEGIN{OFS="\t"}
        function flush(){if(active){print chromosome,start,stop;active=0}}
        $3>=40 && $3<=100 {
            base=$2-1
            if(active && chromosome==$1 && stop==base){stop=$2}
            else{flush();chromosome=$1;start=base;stop=$2;active=1}
            next
        }
        {flush()}
        END{flush()}' >> "$outdir/depth_40_100.bed"
    done
    for name in gvcf_confident depth_40_100; do
        bedtools sort -g "$OUTPUT_ROOT/chr1-10.genome" -i "$outdir/${name}.bed" > "$outdir/${name}.sorted.bed"
    done
    bedtools intersect -sorted -g "$OUTPUT_ROOT/chr1-10.genome" -a "$outdir/gvcf_confident.sorted.bed" -b "$outdir/depth_40_100.sorted.bed" |
    bedtools intersect -sorted -g "$OUTPUT_ROOT/chr1-10.genome" -a stdin -b "$OUTPUT_ROOT/reference_ACGT.sorted.bed" |
    bedtools sort -g "$OUTPUT_ROOT/chr1-10.genome" -i stdin |
    bedtools merge -i stdin > "$outdir/callable_snps_indels_chr1-10.no_reference_N.bed"
    final_bed="$outdir/callable_snps_indels_chr1-10.no_reference_N.bed"
    callable_bases=$(awk '{n+=$3-$2} END{printf "%.0f",n+0}' "$final_bed")
    printf 'Callable bases (chr1-10): %s\nBED: %s\n' "$callable_bases" "$final_bed" > "$outdir/callable_summary.txt"
    printf '%s\t%s\t%s\n' "$sample" "$callable_bases" "$final_bed" >> "$OUTPUT_ROOT/summary_all_samples.tsv"
done
echo "Finished: $OUTPUT_ROOT/summary_all_samples.tsv"
