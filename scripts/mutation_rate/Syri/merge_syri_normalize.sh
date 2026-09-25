#!/usr/bin/env bash
# Normalize SyRI small variants to one reference, then merge samples.
# Usage: bash merge_syri_normalized.sh J702_flip.fa syri_samples.tsv output_dir [chromosome_map.tsv]
set -euo pipefail

if (( $# < 3 || $# > 4 )); then
    echo "Usage: $0 reference.fa syri_samples.tsv output_dir [chromosome_map.tsv]" >&2
    exit 2
fi

ref=$1
manifest=$2
outdir=$3
chrom_map=${4:-}
threads=${SLURM_CPUS_PER_TASK:-8}

for program in bcftools samtools; do
    command -v "$program" >/dev/null || { echo "Missing: $program" >&2; exit 1; }
done
[[ -s "$ref" && -s "$manifest" ]] || { echo "Reference or manifest missing" >&2; exit 1; }
[[ -z "$chrom_map" || -s "$chrom_map" ]] || { echo "Chromosome map missing" >&2; exit 1; }
[[ -s "${ref}.fai" ]] || samtools faidx "$ref"

mkdir -p "$outdir/normalized"
list="$outdir/normalized_files.list"
: > "$list"
declare -A seen=()
count=0

while IFS=$'\t' read -r sample input extra || [[ -n "${sample:-}" ]]; do
    [[ -z "${sample:-}" || "$sample" == \#* ]] && continue
    [[ -z "${extra:-}" && -n "${input:-}" ]] || {
        echo "Each manifest line needs exactly two tab-separated columns: SAMPLE and VCF" >&2
        exit 1
    }
    [[ "$sample" =~ ^[A-Za-z0-9_.-]+$ && -z "${seen[$sample]+x}" ]] || {
        echo "Duplicate/invalid sample name: $sample" >&2
        exit 1
    }
    seen[$sample]=1
    [[ -s "$input" ]] || { echo "Missing SyRI VCF: $input" >&2; exit 1; }

    # SyRI VCFs must each contain exactly one sample. Name it explicitly.
    mapfile -t names < <(bcftools query -l "$input")
    (( ${#names[@]} == 1 )) || { echo "Expected one sample in $input" >&2; exit 1; }
    namefile="$outdir/normalized/$sample.sample_name.txt"
    printf '%s\n' "$sample" > "$namefile"
    renamed="$outdir/normalized/$sample.reheader.vcf.gz"
    bcftools reheader -s "$namefile" -o "$renamed" "$input"

    output="$outdir/normalized/$sample.syri.norm.vcf.gz"
    # Keep originally  SNPs and sequence-level indels; structural
    # types such as INV and NOTAL are not used for exact DV matching.
    # Rely on bcftools variant typing here; a previous REF/ALT regex filter
    # rejected all records with some bcftools builds.
    if [[ -n "$chrom_map" ]]; then
        bcftools annotate --rename-chrs "$chrom_map" -Ou "$renamed" |
            bcftools view -m2 -M2 -v snps,indels -Ou |
            bcftools norm -f "$ref" -c e -Oz -o "$output"
    else
        bcftools view -m2 -M2 -v snps,indels -Ou "$renamed" |
            bcftools norm -f "$ref" -c e -Oz -o "$output"
    fi
    bcftools index -f -c "$output"
    variant_count=$(bcftools index -n "$output")
    (( variant_count > 0 )) || {
        echo "ERROR: $sample normalized to zero variants; inspect chromosome names and input VCF" >&2
        exit 1
    }
    printf '%s\n' "$output" >> "$list"
    count=$((count + 1))
    echo "Normalized $sample from $input: $variant_count variants"
done < "$manifest"

(( count >= 2 )) || { echo "Provide at least two samples to merge" >&2; exit 1; }
merged="$outdir/syri.all_samples.normalized.merged.vcf.gz"
bcftools merge --threads "$threads" -m none -l "$list" -Oz -o "$merged"
bcftools index -f -c "$merged"
echo "Merged VCF: $merged"
bcftools query -l "$merged"
