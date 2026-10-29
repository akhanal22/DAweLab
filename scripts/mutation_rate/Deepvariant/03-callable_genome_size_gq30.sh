#!/bin/bash
#SBATCH --job-name=callable
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# Usage: from the project directory, bash estimate_callable_from_gvcf.sh
# Callable genome for biallelic SNP analysis on chr1_RagTag..chr10_RagTag.
# Requires Python 3, samtools, and bedtools. 
# Run without arguments from the directory containing deepvariant_against_J702.
# The 4-argument form below is used internally for one sample.

REFERENCE="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"
INPUT_ROOT="deepvariant_against_J702"
OUTPUT_ROOT="callable_biallelic_snps"
SAMPLES=(B73D1 B73D2 B73D4 B73D7c B73D10)

if (( $# == 0 )); then
    [[ -d "$INPUT_ROOT/gvcf" && -d "$INPUT_ROOT/prepared_bams" ]] || {
        echo "Expected $INPUT_ROOT/gvcf and $INPUT_ROOT/prepared_bams beneath: $PWD" >&2
        echo "Change to the directory that contains deepvariant_against_J702 and run again." >&2
        exit 1
    }
    mkdir -p "$OUTPUT_ROOT"
    summary="$OUTPUT_ROOT/summary_all_samples.tsv"
    printf 'Sample\tCallable_bases\tBED\n' > "$summary"
    for sample in "${SAMPLES[@]}"; do
        gvcf_path="$INPUT_ROOT/gvcf/${sample}.chr1-10.g.vcf.gz"
        bam_path="$INPUT_ROOT/prepared_bams/${sample}.sorted.primary_mapq20.RG.bam"
        sample_output="$OUTPUT_ROOT/$sample"
        echo "Processing $sample"
        echo "  gVCF: $gvcf_path"
        echo "  BAM:  $bam_path"
        echo "  REF:  $REFERENCE"
        bash "$0" "$gvcf_path" "$bam_path" "$REFERENCE" "$sample_output"
        callable_bases=$(awk '/^Callable bases/ { print $4 }' "$sample_output/callable_summary.txt")
        printf '%s\t%s\t%s\n' "$sample" "$callable_bases" \
            "$sample_output/callable_biallelic_snps_chr1-10.no_reference_N.bed" >> "$summary"
    done
    echo "All samples finished. Summary: $summary"
    exit 0
fi

if (( $# != 4 )); then
    echo "Usage: bash $0 (no arguments; run where deepvariant_against_J702 exists)" >&2
    exit 2
fi
gvcf=$1
bam=$2
ref=$3
outdir=$4
for input in "$gvcf" "$bam" "$ref"; do
    [[ -s "$input" ]] || { echo "Missing or empty: $input" >&2; exit 1; }
done
[[ -s "${bam}.bai" || -s "${bam%.bam}.bai" || -s "${bam}.csi" ]] || {
    echo "No BAM index found for: $bam" >&2
    echo "On a compute node, run: samtools index -c '$bam'" >&2
    exit 1
}
for command in python3 samtools bedtools; do
    command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 1; }
done
mkdir -p "$outdir"
[[ -s "${ref}.fai" ]] || samtools faidx "$ref"
awk 'BEGIN { OFS="\t" } /^chr([1-9]|10)_RagTag\t/ { print $1, $2 }' \
    "${ref}.fai" > "$outdir/chr1-10.genome"
[[ -s "$outdir/chr1-10.genome" ]] || {
    echo "No chr1_RagTag..chr10_RagTag in ${ref}.fai" >&2
    exit 1
}

# DeepVariant reference bands: 0/0, GQ >=30, MIN_DP 40..100.
# Biallelic SNP calls: one real A/C/G/T ALT (an extra gVCF <*> is ignored),
# PASS, GQ >=30, DP 40..100, and matching allele fraction.
# Heterozygous 0/1: ALT fraction 0.40..0.60; homozygous ALT 1/1: >=0.95.
# Reference 0/0 calls at biallelic SNP candidate sites need GQ >=30, DP 40..100.
# Multiallelic and indel records are excluded; no call needs AD on 0/0 bands.
# The gVCF alone cannot impose a per-base upper-depth limit on reference bands.
python3 - "$gvcf" "$outdir/gvcf_confident.bed" <<'PY'
import gzip
import sys

source, destination = sys.argv[1:]
chromosomes = {f"chr{i}_RagTag" for i in range(1, 11)}

def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return -1

def supports_genotype(gt, sample_values):
    try:
        depths = [int(depth) for depth in sample_values["AD"].split(",")]
    except (KeyError, ValueError):
        return False
    if len(depths) < 2:
        return False
    ref_reads, alt_reads = depths[:2]  # Ignore the symbolic <*> AD column.
    total = ref_reads + alt_reads
    if total == 0 or ref_reads < 0 or alt_reads < 0:
        return False
    if gt in ("1/1",):
        return alt_reads / total >= 0.95
    if gt in ("0/1", "1/0"):
        return 0.40 <= alt_reads / total <= 0.60
    return False

opener = gzip.open if source.endswith(".gz") else open
sample_columns = None
with opener(source, "rt") as stream, open(destination, "w") as output:
    for line in stream:
        if line.startswith("#CHROM"):
            sample_columns = line.rstrip("\n").split("\t")
            if len(sample_columns) != 10:
                raise SystemExit("Expected a single-sample gVCF; found " + str(len(sample_columns) - 9))
            continue
        if line.startswith("#"):
            continue
        if sample_columns is None:
            raise SystemExit("gVCF has no #CHROM header")
        fields = line.rstrip("\n").split("\t")
        chrom, pos_str, _, ref, alt, _, filt, info, fmt, sample = fields[:10]
        if chrom not in chromosomes:
            continue
        pos = int(pos_str)
        values = dict(zip(fmt.split(":"), sample.split(":")))
        gt = values.get("GT", ".").replace("|", "/")
        gq = number(values.get("GQ"))
        if gq < 30:
            continue
        tags = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
        is_ref_band = alt == "<*>" and gt == "0/0"
        if is_ref_band:
            min_dp = number(values.get("MIN_DP"))
            if min_dp < 0:  # Some single-position reference calls store DP only.
                min_dp = number(values.get("DP"))
            if 40 <= min_dp <= 100:
                end = int(tags.get("END", pos))
                if end < pos:
                    raise SystemExit(f"Invalid END at {chrom}:{pos}")
                output.write(f"{chrom}\t{pos - 1}\t{end}\n")
        else:
            alts = [allele for allele in alt.split(",")
                    if allele not in ("<*>", "<NON_REF>")]
            if len(alts) != 1 or len(ref) != 1 or len(alts[0]) != 1:
                continue
            if ref.upper() not in "ACGT" or alts[0].upper() not in "ACGT":
                continue
            if ref.upper() == alts[0].upper():
                continue
            dp = number(values.get("DP"))
            if not 40 <= dp <= 100:
                continue
            if gt == "0/0" or (filt == "PASS" and supports_genotype(gt, values)):
                output.write(f"{chrom}\t{pos - 1}\t{pos}\n")
PY

# Make a BED of all A/C/G/T reference positions. Thus a reference N of *any*
# length is excluded, including NNNN runs, without losing adjacent bases.
python3 - "$ref" "$outdir/reference_ACGT.bed" <<'PY'
import re
import sys

fasta, destination = sys.argv[1:]
chromosomes = {f"chr{i}_RagTag" for i in range(1, 11)}
name = None
parts = []

def emit(output):
    if name not in chromosomes:
        return
    sequence = "".join(parts)
    for match in re.finditer(r"[ACGTacgt]+", sequence):
        output.write(f"{name}\t{match.start()}\t{match.end()}\n")

with open(fasta) as stream, open(destination, "w") as output:
    for line in stream:
        if line.startswith(">"):
            if name is not None:
                emit(output)
            name = line[1:].split()[0]
            parts = []
        elif name in chromosomes:
            parts.append(line.strip())
    if name is not None:
        emit(output)
PY

# Query only chr1_RagTag..chr10_RagTag from the indexed BAM.
# Compress the per-base BAM depth stream into runs. -aa includes uncovered bases.
# Check that samtools depth settings agree with the reads/quality settings used
# to generate the DeepVariant FORMAT/DP; they are not necessarily identical.
for chromosome_number in {1..10}; do
    samtools depth -aa -q 0 -Q 0 -r "chr${chromosome_number}_RagTag" "$bam"
done | awk '
BEGIN { OFS="\t" }
function flush() { if (active) { print chromosome, start, stop; active=0 } }
$1 ~ /^chr([1-9]|10)_RagTag$/ && $3 >= 40 && $3 <= 100 {
    base = $2 - 1
    if (active && chromosome == $1 && stop == base) { stop = $2 }
    else { flush(); chromosome = $1; start = base; stop = $2; active = 1 }
    next
}
{ flush() }
END { flush() }
' > "$outdir/depth_40_100.bed"

for name in gvcf_confident reference_ACGT depth_40_100; do
    bedtools sort -g "$outdir/chr1-10.genome" -i "$outdir/${name}.bed" > "$outdir/${name}.sorted.bed"
done

bedtools intersect -sorted -g "$outdir/chr1-10.genome" \
    -a "$outdir/gvcf_confident.sorted.bed" \
    -b "$outdir/depth_40_100.sorted.bed" |
bedtools intersect -sorted -g "$outdir/chr1-10.genome" \
    -a stdin -b "$outdir/reference_ACGT.sorted.bed" |
bedtools merge -i stdin > "$outdir/callable_biallelic_snps_chr1-10.no_reference_N.bed"

awk -v file="$outdir/callable_biallelic_snps_chr1-10.no_reference_N.bed" '
    { bases += $3 - $2 }
    END { printf "Callable bases (chr1-10): %.0f\nBED: %s\n", bases, file }
' "$outdir/callable_biallelic_snps_chr1-10.no_reference_N.bed" |
tee "$outdir/callable_summary.txt"
echo "Intermediate BEDs: $outdir"
