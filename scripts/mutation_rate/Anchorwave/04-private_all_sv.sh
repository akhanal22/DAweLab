#!/bin/bash
#SBATCH --job-name=private_anchorwave
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

#######################################################################
# Software
#######################################################################

# Conda is not needed because SAMtools and BCFtools are loaded as modules.
module purge
module load SAMtools/1.23.1-GCC-13.3.0
module load BCFtools/1.23.1-GCC-13.3.0

threads="${SLURM_CPUS_PER_TASK:-8}"

samtools --version | sed -n '1p'
bcftools --version | sed -n '1p'

#######################################################################
# Configuration
#######################################################################

ref="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"

samples=(
    B73D1
    B73D2
    B73D4
    B73D7c1
    B73D10
)

regions="chr1_RagTag,chr2_RagTag,chr3_RagTag,chr4_RagTag,chr5_RagTag,chr6_RagTag,chr7_RagTag,chr8_RagTag,chr9_RagTag,chr10_RagTag"

outdir="private_anchorwave_gvcf_analysis"

normalized_dir="${outdir}/01_gvcf_to_sequence_vcf"
separated_dir="${outdir}/02_SNPs_indels"
private_dir="${outdir}/03_private_by_individual"
summary_dir="${outdir}/04_summary"
merged_dir="${outdir}/05_merged_all_variants"

mkdir -p \
    "$normalized_dir" \
    "$separated_dir" \
    "$private_dir" \
    "$summary_dir" \
    "$merged_dir"

echo "AnchorWave private-variant workflow started"
date

#######################################################################
# STEP 1: Check reference
#######################################################################

echo "STEP 1: Checking reference"

[[ -s "$ref" ]] || {
    echo "ERROR: Missing reference: $ref" >&2
    exit 1
}

[[ -s "${ref}.fai" ]] || samtools faidx "$ref"

for chromosome in \
    chr1_RagTag chr2_RagTag chr3_RagTag chr4_RagTag chr5_RagTag \
    chr6_RagTag chr7_RagTag chr8_RagTag chr9_RagTag chr10_RagTag
do
    if ! awk -v chromosome="$chromosome" '
        $1 == chromosome {
            found = 1
            exit
        }
        END {
            exit(found ? 0 : 1)
        }
    ' "${ref}.fai"
    then
        echo "ERROR: ${chromosome} is absent from ${ref}.fai" >&2
        exit 1
    fi
done

#######################################################################
# STEP 2: Check input VCFs
#######################################################################

echo "STEP 2: Checking AnchorWave VCFs"

input_summary="${summary_dir}/input_VCF_summary.tsv"

printf \
    "Sample\tInput_VCF\tVCF_sample_name\tTotal_records\tGT0_records\tGT_ALT_records\n" \
    > "$input_summary"

for sample in "${samples[@]}"
do
    input="${sample}.chrOnly.g.vcf.gz"

    [[ -s "$input" ]] || {
        echo "ERROR: Missing input: $input" >&2
        exit 1
    }

    bcftools view -h "$input" >/dev/null

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$input"

    vcf_sample_count=$(bcftools query -l "$input" | wc -l)

    if [[ "$vcf_sample_count" -ne 1 ]]
    then
        echo "ERROR: ${input} contains ${vcf_sample_count} samples; expected one" >&2
        exit 1
    fi

    vcf_sample=$(bcftools query -l "$input")

    if [[ "$vcf_sample" != "$sample" ]]
    then
        echo "ERROR: Sample-name mismatch in ${input}" >&2
        echo "Expected: ${sample}" >&2
        echo "Found:    ${vcf_sample}" >&2
        exit 1
    fi

    total_records=$(bcftools view -H "$input" | wc -l)

    gt0_records=$(
        bcftools view \
            -H \
            -i 'GT="ref"' \
            "$input" |
        wc -l
    )

    gt_alt_records=$(
        bcftools view \
            -H \
            -i 'GT="alt"' \
            "$input" |
        wc -l
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$input" \
        "$vcf_sample" \
        "$total_records" \
        "$gt0_records" \
        "$gt_alt_records" \
        >> "$input_summary"

    echo "${sample}: total=${total_records}; GT0=${gt0_records}; ALT=${gt_alt_records}"
done

#######################################################################
# STEP 3: Convert MAF-derived gVCFs to sequence-variant VCFs
#######################################################################

echo "STEP 3: Removing gVCF reference blocks and normalizing sequence variants"

normalization_summary="${summary_dir}/normalization_summary.tsv"

printf \
    "Sample\tInput_records\tGT0_reference_blocks\tGT1_sequence_records\tFinal_SNP_indel_records\n" \
    > "$normalization_summary"

for sample in "${samples[@]}"
do
    input="${sample}.chrOnly.g.vcf.gz"

    output="${normalized_dir}/${sample}.anchorwave.sequence_snps_indels.chr1-10.vcf.gz"

    input_count=$(bcftools view -H "$input" | wc -l)

    input_refblock_count=$(
        bcftools view \
            -H \
            -r "$regions" \
            -i 'GT="ref" && ALT="<NON_REF>"' \
            "$input" |
        wc -l
    )

    input_alt_count=$(
        bcftools view \
            -H \
            -r "$regions" \
            -i 'GT="alt"' \
            "$input" |
        wc -l
    )

    echo "Processing ${sample}"

    # The MAF-to-gVCF conversion represents:
    #   sequence variants as ALT=sequence,<NON_REF> with haploid GT=1
    #   reference blocks as ALT=<NON_REF> with haploid GT=0 and INFO/END
    #
    # 1. GT="alt" removes the GT=0 reference blocks.
    # 2. Remove INFO and all FORMAT fields except GT because ASM coordinates,
    #    artificial AD/DP/PL values, and gVCF END are not variant evidence.
    # 3. Split sequence ALT from the technical <NON_REF> allele and normalize.
    # 4. Retain only the split allele whose genotype remains ALT, explicitly
    #    excluding all symbolic gVCF alleles.
    # 5. Atomize complex substitutions so MNPs become individual SNVs.
    # 6. Retain only biallelic sequence SNVs and indels of any size.
    # 7. Remove exact duplicates and sort because left alignment can move an
    #    indel before a preceding input record.
    bcftools view \
        --threads "$threads" \
        -r "$regions" \
        -i 'GT="alt"' \
        -Ou \
        "$input" |
    bcftools annotate \
        --threads "$threads" \
        -x 'INFO,^FORMAT/GT' \
        -Ou |
    bcftools norm \
        --threads "$threads" \
        -m -any \
        -f "$ref" \
        -Ou |
    bcftools view \
        --threads "$threads" \
        -i '
            GT="alt" &&
            ALT!="<NON_REF>" &&
            ALT!="<*>" &&
            ALT!="*" &&
            ALT!="."
        ' \
        -Ou |
    bcftools norm \
        --threads "$threads" \
        -a \
        --atom-overlaps . \
        -m -any \
        -f "$ref" \
        -Ou |
    bcftools view \
        --threads "$threads" \
        -m2 \
        -M2 \
        -v snps,indels \
        -i 'GT="alt"' \
        -Ou |
    bcftools norm \
        --threads "$threads" \
        -d exact \
        -Ou |
    bcftools sort \
        -Oz \
        -o "$output"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$output"

    normalized_count=$(bcftools view -H "$output" | wc -l)

    invalid_sequence_records=$(
        bcftools view \
            -H \
            -e '
                N_ALT=1 &&
                (TYPE="snp" || TYPE="indel") &&
                GT="alt"
            ' \
            "$output" |
        wc -l
    )

    if [[ "$invalid_sequence_records" -ne 0 ]]
    then
        echo "ERROR: ${sample} cleaned VCF contains ${invalid_sequence_records} non-SNP/indel, symbolic, multiallelic, or non-ALT records" >&2
        exit 1
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$input_count" \
        "$input_refblock_count" \
        "$input_alt_count" \
        "$normalized_count" \
        >> "$normalization_summary"

    echo "${sample}: ${input_refblock_count} reference blocks removed; ${normalized_count} SNP/indel records retained"
done

#######################################################################
# STEP 4: Separate SNPs and indels of all sizes
#######################################################################

echo "STEP 4: Separating SNPs and indels of all sizes"

selection_summary="${summary_dir}/variant_selection_summary.tsv"

printf \
    "Sample\tSNPs\tIndels_all_sizes\tInsertions\tDeletions\n" \
    > "$selection_summary"

for sample in "${samples[@]}"
do
    normalized="${normalized_dir}/${sample}.anchorwave.sequence_snps_indels.chr1-10.vcf.gz"

    snps="${separated_dir}/${sample}.anchorwave.snps.chr1-10.vcf.gz"

    indels="${separated_dir}/${sample}.anchorwave.indels.all_sizes.chr1-10.vcf.gz"

    bcftools view \
        --threads "$threads" \
        -m2 \
        -M2 \
        -v snps \
        -i 'GT="alt"' \
        -Oz \
        -o "$snps" \
        "$normalized"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$snps"

    bcftools view \
        --threads "$threads" \
        -m2 \
        -M2 \
        -v indels \
        -i 'GT="alt" && ABS(ILEN)>=1' \
        -Oz \
        -o "$indels" \
        "$normalized"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$indels"

    snp_count=$(bcftools view -H "$snps" | wc -l)

    indel_count=$(bcftools view -H "$indels" | wc -l)

    insertion_count=$(
        bcftools view \
            -H \
            -i 'ILEN>0' \
            "$indels" |
        wc -l
    )

    deletion_count=$(
        bcftools view \
            -H \
            -i 'ILEN<0' \
            "$indels" |
        wc -l
    )

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$snp_count" \
        "$indel_count" \
        "$insertion_count" \
        "$deletion_count" \
        >> "$selection_summary"
done

#######################################################################
# STEP 5: Prepare exact-allele comparison inputs
#######################################################################

echo "STEP 5: Preparing exact-allele comparison inputs"

snp_inputs=()
indel_inputs=()

for sample in "${samples[@]}"
do
    snp_inputs+=(
        "${separated_dir}/${sample}.anchorwave.snps.chr1-10.vcf.gz"
    )

    indel_inputs+=(
        "${separated_dir}/${sample}.anchorwave.indels.all_sizes.chr1-10.vcf.gz"
    )
done

#######################################################################
# STEP 6: Identify alleles present in exactly one assembly
#######################################################################

echo "STEP 6: Selecting private sequence alleles"

# -n=1:
#   The normalized allele occurs in exactly one of the five VCFs.
#
# -c none:
#   Alleles match only when CHROM, POS, REF and ALT match exactly.
#
# -w:
#   Write the record from the corresponding input VCF.

for array_index in "${!samples[@]}"
do
    sample="${samples[$array_index]}"

    # BCFtools input numbering starts at 1.
    isec_index=$((array_index + 1))

    private_snps="${private_dir}/${sample}.private_anchorwave.snps.chr1-10.vcf.gz"

    private_indels="${private_dir}/${sample}.private_anchorwave.indels.all_sizes.chr1-10.vcf.gz"

    echo "Selecting private SNPs for ${sample}"

    bcftools isec \
        -n=1 \
        -c none \
        -w "$isec_index" \
        "${snp_inputs[@]}" \
        -Oz \
        -o "$private_snps"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$private_snps"

    echo "Selecting private indels for ${sample}"

    bcftools isec \
        -n=1 \
        -c none \
        -w "$isec_index" \
        "${indel_inputs[@]}" \
        -Oz \
        -o "$private_indels"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$private_indels"
done

#######################################################################
# STEP 7: Count private sequence alleles
#######################################################################

echo "STEP 7: Counting private sequence alleles"

count_summary="${summary_dir}/private_anchorwave_variant_counts.tsv"

printf \
    "Sample\tPrivate_SNPs\tPrivate_indels_all_sizes\tPrivate_insertions\tPrivate_deletions\n" \
    > "$count_summary"

for sample in "${samples[@]}"
do
    private_snps="${private_dir}/${sample}.private_anchorwave.snps.chr1-10.vcf.gz"

    private_indels="${private_dir}/${sample}.private_anchorwave.indels.all_sizes.chr1-10.vcf.gz"

    snp_count=$(bcftools view -H "$private_snps" | wc -l)

    indel_count=$(bcftools view -H "$private_indels" | wc -l)

    insertion_count=$(
        bcftools view \
            -H \
            -i 'ILEN>0' \
            "$private_indels" |
        wc -l
    )

    deletion_count=$(
        bcftools view \
            -H \
            -i 'ILEN<0' \
            "$private_indels" |
        wc -l
    )

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$snp_count" \
        "$indel_count" \
        "$insertion_count" \
        "$deletion_count" \
        >> "$count_summary"
done

#######################################################################
# STEP 8: Create readable TSV files
#######################################################################

echo "STEP 8: Creating TSV files"

for sample in "${samples[@]}"
do
    private_snps="${private_dir}/${sample}.private_anchorwave.snps.chr1-10.vcf.gz"

    private_indels="${private_dir}/${sample}.private_anchorwave.indels.all_sizes.chr1-10.vcf.gz"

    snp_table="${private_dir}/${sample}.private_anchorwave.snps.chr1-10.tsv"

    indel_table="${private_dir}/${sample}.private_anchorwave.indels.all_sizes.chr1-10.tsv"

    printf \
        "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tGT\n" \
        > "$snp_table"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' \
        "$private_snps" \
        >> "$snp_table"

    printf \
        "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tGT\tSIGNED_ILEN\tABS_SIZE\tTYPE\n" \
        > "$indel_table"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' \
        "$private_indels" |
    awk '
        BEGIN {
            FS = OFS = "\t"
        }
        {
            # For a normalized biallelic sequence indel:
            # signed indel length = ALT length - REF length.
            ilen = length($4) - length($3)
            size = ilen

            if (size < 0) {
                size = -size
            }

            if (ilen > 0) {
                type = "insertion"
            } else if (ilen < 0) {
                type = "deletion"
            } else {
                type = "complex"
            }

            print $0, ilen, size, type
        }
    ' >> "$indel_table"
done

#######################################################################
# STEP 9: Merge private SNPs and all-size indels
#######################################################################

echo "STEP 9: Merging all private variant types"

merged_inputs=()
merged_summary="${summary_dir}/merged_all_variant_counts.tsv"

printf \
    "Sample\tMerged_SNPs\tMerged_indels_all_sizes\tMerged_total_records\n" \
    > "$merged_summary"

for sample in "${samples[@]}"
do
    private_snps="${private_dir}/${sample}.private_anchorwave.snps.chr1-10.vcf.gz"
    private_indels="${private_dir}/${sample}.private_anchorwave.indels.all_sizes.chr1-10.vcf.gz"
    merged_sample="${merged_dir}/${sample}.private_anchorwave.snps_all_indels.chr1-10.vcf.gz"

    # SNP and indel VCFs overlap in genomic coordinate ranges, so concatenate
    # with -a and immediately restore coordinate order with bcftools sort.
    bcftools concat \
        --threads "$threads" \
        -a \
        -Ou \
        "$private_snps" \
        "$private_indels" |
    bcftools sort \
        -Oz \
        -o "$merged_sample"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$merged_sample"

    merged_inputs+=("$merged_sample")

    merged_snp_count=$(bcftools view -v snps -H "$merged_sample" | wc -l)
    merged_indel_count=$(bcftools view -v indels -H "$merged_sample" | wc -l)
    merged_total_count=$(bcftools view -H "$merged_sample" | wc -l)

    printf "%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$merged_snp_count" \
        "$merged_indel_count" \
        "$merged_total_count" \
        >> "$merged_summary"
done

all_samples_merged="${merged_dir}/five_lines.private_anchorwave.snps_all_indels.chr1-10.merged.vcf.gz"

# Keep different alleles at the same position as separate biallelic records.
bcftools merge \
    --threads "$threads" \
    -m none \
    -Oz \
    -o "$all_samples_merged" \
    "${merged_inputs[@]}"

bcftools index \
    --threads "$threads" \
    -f \
    -c \
    "$all_samples_merged"

#######################################################################
# STEP 10: Validate final files
#######################################################################

echo "STEP 10: Validating private and merged variants"

for sample in "${samples[@]}"
do
    private_snps="${private_dir}/${sample}.private_anchorwave.snps.chr1-10.vcf.gz"

    private_indels="${private_dir}/${sample}.private_anchorwave.indels.all_sizes.chr1-10.vcf.gz"

    invalid_snps=$(
        bcftools view \
            -H \
            -e '
                N_ALT=1 &&
                TYPE="snp" &&
                GT="alt"
            ' \
            "$private_snps" |
        wc -l
    )

    invalid_indels=$(
        bcftools view \
            -H \
            -e '
                N_ALT=1 &&
                TYPE="indel" &&
                GT="alt" &&
                ABS(ILEN)>=1
            ' \
            "$private_indels" |
        wc -l
    )

    if [[ "$invalid_snps" -ne 0 ]]
    then
        echo "ERROR: ${sample} has ${invalid_snps} invalid SNP records" >&2
        exit 1
    fi

    if [[ "$invalid_indels" -ne 0 ]]
    then
        echo "ERROR: ${sample} has ${invalid_indels} invalid indel records" >&2
        exit 1
    fi

    expected_total=$(
        awk -v sample="$sample" '
            BEGIN { FS = "\t" }
            $1 == sample { print $2 + $3 }
        ' "$count_summary"
    )

    merged_sample="${merged_dir}/${sample}.private_anchorwave.snps_all_indels.chr1-10.vcf.gz"
    observed_total=$(bcftools view -H "$merged_sample" | wc -l)

    if [[ "$observed_total" -ne "$expected_total" ]]
    then
        echo "ERROR: ${sample} merged record count is ${observed_total}; expected ${expected_total}" >&2
        exit 1
    fi
done

merged_sample_count=$(bcftools query -l "$all_samples_merged" | wc -l)

if [[ "$merged_sample_count" -ne "${#samples[@]}" ]]
then
    echo "ERROR: Final merged VCF contains ${merged_sample_count} samples; expected ${#samples[@]}" >&2
    exit 1
fi

echo
echo "Input VCF summary:"
column -t "$input_summary" || cat "$input_summary"

echo
echo "Normalization summary:"
column -t "$normalization_summary" || cat "$normalization_summary"

echo
echo "Variant-selection summary:"
column -t "$selection_summary" || cat "$selection_summary"

echo
echo "Private AnchorWave variant counts:"
column -t "$count_summary" || cat "$count_summary"

echo
echo "Merged SNP and all-indel counts:"
column -t "$merged_summary" || cat "$merged_summary"

echo
echo "Final five-sample merged VCF:"
echo "$all_samples_merged"

echo
echo "Workflow completed successfully"
date
