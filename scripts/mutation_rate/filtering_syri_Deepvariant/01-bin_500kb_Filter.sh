#!/bin/bash
#SBATCH --job-name=compare_private_500kb
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# Run from the directory containing vcf_to_compare, or supply that directory:
#   sbatch deepvariant_raw.sh [vcf_to_compare]
#
# For each of three already merged VCFs, separately:
#   - retain SNP/indel records carried by exactly one sample, regardless of
#     FILTER value or whether the record has multiple ALT alleles;
#   - count private SNPs in each sample in fixed, non-overlapping 500-kb bins;
#   - remove that sample's SNPs AND indels in bins with >5 private SNPs;
#   - merge the individually filtered sample VCFs again.
# Inputs in the filtered_0.95hom_40-60her_dv directory are never used.

INPUT_DIR=${1:-vcf_to_compare}
OUTDIR="${INPUT_DIR}/three_callers_private_500kb"
THREADS=${SLURM_CPUS_PER_TASK:-16}
WINDOW=500000
MAX_SNPS=5

command -v bcftools >/dev/null || { echo "ERROR: bcftools is not on PATH" >&2; exit 1; }
[[ -d "${INPUT_DIR}" ]] || { echo "ERROR: input directory not found: ${INPUT_DIR}" >&2; exit 1; }
mkdir -p "${OUTDIR}"

# The three actual VCFs listed in vcf_to_compare; .csi files are indexes.
CALLERS=(deepvariant anchorwave syri)
INPUTS=(
    "${INPUT_DIR}/five_lines.deepvariant.private_snps_indels.vcf.gz"
    "${INPUT_DIR}/five_lines.private_anchorwave.snps_all_indels.chr1-10.merged.vcf.gz"
    "${INPUT_DIR}/syri.private_snps_indels.exactly_one_sample.vcf.gz"
)

# Counts each *sample* once if its GT contains any ALT allele (1, 2, etc.).
# Handles 0/1, 1/2, 2/2 and partially missing GTs such as ./1 or SyRI .1.
ALT_AWK='
BEGIN { FS=OFS="\t" }
/^#/ { if (!NO_HEADER) print; next }
{
    gt_col=0
    nf=split($9, format_fields, ":")
    for (k=1; k<=nf; k++)
        if (format_fields[k]=="GT") gt_col=k
    if (!gt_col) {
        print "ERROR: VCF row has no FORMAT/GT: " $1 ":" $2 > "/dev/stderr"
        exit 1
    }
    carriers=0
    for (i=10; i<=NF; i++) {
        split($i, sample_fields, ":")
        na=split(sample_fields[gt_col], alleles, /[^0-9]+/)
        for (j=1; j<=na; j++)
            if (alleles[j] ~ /^[0-9]+$/ && (alleles[j]+0)>0) { carriers++; break }
    }
    if (INVERT ? carriers!=REQUIRE : carriers==REQUIRE) print
}'

SUMMARY="${OUTDIR}/three_callers_500kb_filter_summary.tsv"
printf 'CALLER\tSAMPLE\tPRIVATE_SNPS_BEFORE\tPRIVATE_INDELS_BEFORE\tBAD_BINS\tPRIVATE_SNPS_AFTER\tPRIVATE_INDELS_AFTER\n' > "${SUMMARY}"

for idx in "${!CALLERS[@]}"; do
    CALLER=${CALLERS[$idx]}
    INPUT_VCF=${INPUTS[$idx]}
    [[ -s "${INPUT_VCF}" ]] || { echo "ERROR: missing input VCF: ${INPUT_VCF}" >&2; exit 1; }

    CALLER_DIR="${OUTDIR}/${CALLER}"
    PRIVATE_DIR="${CALLER_DIR}/private_per_sample"
    BAD_DIR="${CALLER_DIR}/bad_500kb_bins"
    FILTERED_DIR="${CALLER_DIR}/filtered_per_sample"
    mkdir -p "${PRIVATE_DIR}" "${BAD_DIR}" "${FILTERED_DIR}"

    mapfile -t SAMPLES < <(bcftools query -l "${INPUT_VCF}")
    (( ${#SAMPLES[@]} >= 2 )) || {
        echo "ERROR: expected a multisample VCF: ${INPUT_VCF}" >&2
        exit 1
    }
    echo "Processing ${CALLER}: ${INPUT_VCF} (${#SAMPLES[@]} samples)"

# Recheck privacy for a consistent starting point. Retain the original
# FILTER values and all ALT alleles; no PASS or biallelic restriction.
    # No reference-based normalization is repeated on these merged files.
    PRIVATE_ALL="${CALLER_DIR}/${CALLER}.private.snps_indels.vcf.gz"
    bcftools view -v snps,indels -Ov "${INPUT_VCF}" |
        awk -v REQUIRE=1 "${ALT_AWK}" |
        bcftools view -Oz -o "${PRIVATE_ALL}" -
    bcftools index -f -c "${PRIVATE_ALL}"

    FILTERED_FILES=()
    for SAMPLE in "${SAMPLES[@]}"; do
        SAMPLE_PRIVATE="${PRIVATE_DIR}/${SAMPLE}.private.snps_indels.vcf.gz"
        SNP_POSITIONS="${PRIVATE_DIR}/${SAMPLE}.private_snp_positions.tsv"
        BAD_WITH_COUNT="${BAD_DIR}/${SAMPLE}.bad_500kb_bins.with_snp_count.bed"
        BAD_BED="${BAD_DIR}/${SAMPLE}.bad_500kb_bins.bed"
        SAMPLE_FILTERED="${FILTERED_DIR}/${SAMPLE}.private.500kb.filtered.snps_indels.vcf.gz"

        # Retain only sites at which this sample itself carries ALT.
        bcftools view -s "${SAMPLE}" -Ov "${PRIVATE_ALL}" |
            awk -v REQUIRE=1 "${ALT_AWK}" |
            bcftools view -Oz -o "${SAMPLE_PRIVATE}" -
        bcftools index -f -c "${SAMPLE_PRIVATE}"

        bcftools view -v snps -Ou "${SAMPLE_PRIVATE}" |
            bcftools query -f '%CHROM\t%POS\n' > "${SNP_POSITIONS}"

        # BED is 0-based and half-open. Only a variant's VCF POS determines
        # its bin, consistent with bcftools view -T below.
        awk -v W="${WINDOW}" -v MAX="${MAX_SNPS}" '
            BEGIN { OFS="\t" }
            {
                start=int(($2-1)/W)*W
                end=start+W
                key=$1 OFS start OFS end
                count[key]++
            }
            END {
                for (key in count)
                    if (count[key]>MAX) print key, count[key]
            }
        ' "${SNP_POSITIONS}" | sort -k1,1V -k2,2n > "${BAD_WITH_COUNT}"
        cut -f1-3 "${BAD_WITH_COUNT}" > "${BAD_BED}"

        if [[ -s "${BAD_BED}" ]]; then
            bcftools view -T "^${BAD_BED}" -Oz -o "${SAMPLE_FILTERED}" "${SAMPLE_PRIVATE}"
        else
            bcftools view -Oz -o "${SAMPLE_FILTERED}" "${SAMPLE_PRIVATE}"
        fi
        bcftools index -f -c "${SAMPLE_FILTERED}"
        FILTERED_FILES+=("${SAMPLE_FILTERED}")

        SNP_BEFORE=$(bcftools view -H -v snps "${SAMPLE_PRIVATE}" | wc -l)
        INDEL_BEFORE=$(bcftools view -H -v indels "${SAMPLE_PRIVATE}" | wc -l)
        BAD_COUNT=$(wc -l < "${BAD_BED}")
        SNP_AFTER=$(bcftools view -H -v snps "${SAMPLE_FILTERED}" | wc -l)
        INDEL_AFTER=$(bcftools view -H -v indels "${SAMPLE_FILTERED}" | wc -l)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${CALLER}" "${SAMPLE}" "${SNP_BEFORE}" "${INDEL_BEFORE}" \
            "${BAD_COUNT}" "${SNP_AFTER}" "${INDEL_AFTER}" >> "${SUMMARY}"
    done

    FINAL="${CALLER_DIR}/${CALLER}.private.500kb.filtered.snps_indels.merged.vcf.gz"
    # Keep absent calls missing; an absent SyRI/AnchorWave call is not 0/0.
    bcftools merge --threads "${THREADS}" -m none -Oz -o "${FINAL}" "${FILTERED_FILES[@]}"
    bcftools index -f -c "${FINAL}"

    INVALID=$(bcftools view -H "${FINAL}" |
        awk -v REQUIRE=1 -v INVERT=1 -v NO_HEADER=1 "${ALT_AWK}" | wc -l)
    (( INVALID == 0 )) || {
        echo "ERROR: ${CALLER} final VCF has ${INVALID} non-private records" >&2
        exit 1
    }
    echo "${CALLER} final VCF: ${FINAL}"
    echo "${CALLER} SNPs: $(bcftools view -H -v snps "${FINAL}" | wc -l)"
    echo "${CALLER} indels: $(bcftools view -H -v indels "${FINAL}" | wc -l)"
done

echo "Per-sample summary: ${SUMMARY}"
