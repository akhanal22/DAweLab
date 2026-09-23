#!/bin/bash
#SBATCH --job-name=private_deepvariant
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

source ~/ENTER/etc/profile.d/conda.sh
conda activate base

module load SAMtools/1.23.1-GCC-13.3.0
module load BCFtools/1.23.1-GCC-13.3.0

threads="${SLURM_CPUS_PER_TASK:-8}"

#######################################################################
# Configuration
#######################################################################

ref="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"

all_samples=(
    B73D1
    B73D2
    B73D4
    B73D7c
    B73D10
    J702
)

ma_samples=(
    B73D1
    B73D2
    B73D4
    B73D7c
    B73D10
)

regions="chr1_RagTag,chr2_RagTag,chr3_RagTag,chr4_RagTag,chr5_RagTag,chr6_RagTag,chr7_RagTag,chr8_RagTag,chr9_RagTag,chr10_RagTag"

outdir="private_deepvariant_analysis"

normalized_dir="${outdir}/01_normalized"
separated_dir="${outdir}/02_biallelic_variants"
filtered_dir="${outdir}/03_filtered_MA_lines"
j702_dir="${outdir}/04_unfiltered_J702"
noj702_dir="${outdir}/05_noJ702"
merged_dir="${outdir}/06_merged"
private_dir="${outdir}/07_private"
individual_dir="${outdir}/08_private_by_individual"
summary_dir="${outdir}/09_summary"

mkdir -p \
    "$normalized_dir" \
    "$separated_dir" \
    "$filtered_dir" \
    "$j702_dir" \
    "$noj702_dir" \
    "$merged_dir" \
    "$private_dir" \
    "$individual_dir" \
    "$summary_dir"

echo "Workflow started"
date

#######################################################################
# Check reference
#######################################################################

[[ -s "$ref" ]] || {
    echo "ERROR: Missing reference: $ref" >&2
    exit 1
}

[[ -s "${ref}.fai" ]] || samtools faidx "$ref"

for chromosome in \
    chr1_RagTag chr2_RagTag chr3_RagTag chr4_RagTag chr5_RagTag \
    chr6_RagTag chr7_RagTag chr8_RagTag chr9_RagTag chr10_RagTag
do
    if ! grep -Fxq "$chromosome" <(cut -f1 "${ref}.fai")
    then
        echo "ERROR: ${chromosome} is absent from ${ref}.fai" >&2
        exit 1
    fi
done

#######################################################################
# STEP 1: Check VCF files and sample names
#######################################################################

echo "STEP 1: Checking DeepVariant VCF files"

for sample in "${all_samples[@]}"
do
    input="${sample}.deepvariant.vcf.gz"

    [[ -s "$input" ]] || {
        echo "ERROR: Missing input: $input" >&2
        exit 1
    }

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$input"

    vcf_sample=$(bcftools query -l "$input")

    if [[ "$vcf_sample" != "$sample" ]]
    then
        echo "ERROR: Sample-name mismatch" >&2
        echo "File:          $input" >&2
        echo "Expected:      $sample" >&2
        echo "Found in VCF:  $vcf_sample" >&2
        exit 1
    fi

    echo "${input}: sample name verified as ${sample}"
done

#######################################################################
# STEP 2: Retain originally biallelic records, then normalize
#######################################################################

echo "STEP 2: Removing original multiallelic records before normalization"

for sample in "${all_samples[@]}"
do
    input="${sample}.deepvariant.vcf.gz"
    output="${normalized_dir}/${sample}.deepvariant.norm.chr1-10.vcf.gz"

    echo "Selecting originally biallelic records and normalizing ${sample}"

    # IMPORTANT:
    # Select records containing exactly one REF and one ALT allele from the
    # ORIGINAL DeepVariant VCF before normalization. This removes original
    # multiallelic genotypes such as 1/2 for both SNPs and indels.
    #
    # Do not use "bcftools norm -m -any" here. Splitting first would convert
    # an original 1/2 genotype into separate 1/0 and 0/1 records, which would
    # then look biallelic and could be counted as independent mutations.
    bcftools view \
        --threads "$threads" \
        -r "$regions" \
        -m2 \
        -M2 \
        -Ou \
        "$input" |
    bcftools norm \
        --threads "$threads" \
        -f "$ref" \
        -d exact \
        -Oz \
        -o "$output"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$output"

    echo -n "${sample} originally biallelic normalized records: "
    bcftools view -H "$output" | wc -l
done

#######################################################################
# STEP 3: Select biallelic SNPs and biallelic 1–50 bp indels
#######################################################################

echo "STEP 3: Separating SNPs and 1–50 bp indels"

selection_summary="${summary_dir}/variant_selection_summary.tsv"

printf \
    "Sample\tBiallelic_SNPs\tBiallelic_indels_1_50bp\n" \
    > "$selection_summary"

for sample in "${all_samples[@]}"
do
    normalized="${normalized_dir}/${sample}.deepvariant.norm.chr1-10.vcf.gz"

    snps="${separated_dir}/${sample}.deepvariant.norm.biallelic_snps.chr1-10.vcf.gz"

    indels="${separated_dir}/${sample}.deepvariant.norm.biallelic_indels.size1-50.chr1-10.vcf.gz"

    # Select biallelic SNPs
    bcftools view \
        --threads "$threads" \
        -m2 \
        -M2 \
        -v snps \
        -Oz \
        -o "$snps" \
        "$normalized"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$snps"

    # Select biallelic indels of 1–50 bp
    bcftools view \
        --threads "$threads" \
        -m2 \
        -M2 \
        -v indels \
        -i 'ABS(ILEN)>=1 && ABS(ILEN)<=50' \
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

    printf "%s\t%s\t%s\n" \
        "$sample" \
        "$snp_count" \
        "$indel_count" \
        >> "$selection_summary"
done

#######################################################################
# STEP 4: Prepare normalized, unfiltered J702 comparison callsets
#######################################################################

echo "STEP 4: Preparing unfiltered J702 comparison files"

j702_snps="${j702_dir}/J702.deepvariant.norm.unfiltered.biallelic_snps.chr1-10.vcf.gz"

j702_indels="${j702_dir}/J702.deepvariant.norm.unfiltered.biallelic_indels.size1-50.chr1-10.vcf.gz"

bcftools view \
    --threads "$threads" \
    -Oz \
    -o "$j702_snps" \
    "${separated_dir}/J702.deepvariant.norm.biallelic_snps.chr1-10.vcf.gz"

bcftools index \
    --threads "$threads" \
    -f \
    -c \
    "$j702_snps"

bcftools view \
    --threads "$threads" \
    -Oz \
    -o "$j702_indels" \
    "${separated_dir}/J702.deepvariant.norm.biallelic_indels.size1-50.chr1-10.vcf.gz"

bcftools index \
    --threads "$threads" \
    -f \
    -c \
    "$j702_indels"

echo -n "Unfiltered normalized J702 SNPs: "
bcftools view -H "$j702_snps" | wc -l

echo -n "Unfiltered normalized J702 1–50 bp indels: "
bcftools view -H "$j702_indels" | wc -l

#######################################################################
# STEP 5: Filter MA-line variants
#######################################################################

echo "STEP 5: Applying rigorous filters to MA lines"

filter_summary="${summary_dir}/rigorous_filtering_summary.tsv"

printf \
    "Sample\tVariant_type\tBefore_filtering\tAfter_filtering\n" \
    > "$filter_summary"

for sample in "${ma_samples[@]}"
do
    for variant_type in snps indels
    do
        if [[ "$variant_type" == "snps" ]]
        then
            input="${separated_dir}/${sample}.deepvariant.norm.biallelic_snps.chr1-10.vcf.gz"

            output="${filtered_dir}/${sample}.deepvariant.norm.biallelic_snps.DP40-100.GQ30.rigorous.chr1-10.vcf.gz"
        else
            input="${separated_dir}/${sample}.deepvariant.norm.biallelic_indels.size1-50.chr1-10.vcf.gz"

            output="${filtered_dir}/${sample}.deepvariant.norm.biallelic_indels.size1-50.DP40-100.GQ30.rigorous.chr1-10.vcf.gz"
        fi

        echo "Filtering ${sample} ${variant_type}"

        bcftools view \
            --threads "$threads" \
            -i '
                FILTER="PASS" &&
                FORMAT/DP>=40 &&
                FORMAT/DP<=100 &&
                FORMAT/GQ>=30 &&
                FORMAT/AD[0:0]!="." &&
                FORMAT/AD[0:1]!="." &&
                (
                    FORMAT/AD[0:0] +
                    FORMAT/AD[0:1]
                ) > 0 &&
                (
                    (
                        GT="AA" &&
                        FORMAT/AD[0:1] /
                        (
                            FORMAT/AD[0:0] +
                            FORMAT/AD[0:1]
                        ) >= 0.90
                    )
                    ||
                    (
                        GT="het" &&
                        FORMAT/AD[0:1] /
                        (
                            FORMAT/AD[0:0] +
                            FORMAT/AD[0:1]
                        ) >= 0.30 &&
                        FORMAT/AD[0:1] /
                        (
                            FORMAT/AD[0:0] +
                            FORMAT/AD[0:1]
                        ) <= 0.70
                    )
                )
            ' \
            -Oz \
            -o "$output" \
            "$input"

        bcftools index \
            --threads "$threads" \
            -f \
            -c \
            "$output"

        before_count=$(bcftools view -H "$input" | wc -l)
        after_count=$(bcftools view -H "$output" | wc -l)

        printf "%s\t%s\t%s\t%s\n" \
            "$sample" \
            "$variant_type" \
            "$before_count" \
            "$after_count" \
            >> "$filter_summary"

        echo "${sample} ${variant_type}: ${after_count} retained"
    done
done

#######################################################################
# STEP 6: Remove exact variants represented in unfiltered J702
#######################################################################

echo "STEP 6: Removing J702-shared variants"

j702_removal_summary="${summary_dir}/J702_removal_summary.tsv"

printf \
    "Sample\tVariant_type\tBefore_J702_removal\tAfter_J702_removal\tRemoved\n" \
    > "$j702_removal_summary"

for sample in "${ma_samples[@]}"
do
    for variant_type in snps indels
    do
        if [[ "$variant_type" == "snps" ]]
        then
            input="${filtered_dir}/${sample}.deepvariant.norm.biallelic_snps.DP40-100.GQ30.rigorous.chr1-10.vcf.gz"

            j702_comparison="$j702_snps"

            output="${noj702_dir}/${sample}.deepvariant.norm.biallelic_snps.DP40-100.GQ30.rigorous.noJ702.chr1-10.vcf.gz"
        else
            input="${filtered_dir}/${sample}.deepvariant.norm.biallelic_indels.size1-50.DP40-100.GQ30.rigorous.chr1-10.vcf.gz"

            j702_comparison="$j702_indels"

            output="${noj702_dir}/${sample}.deepvariant.norm.biallelic_indels.size1-50.DP40-100.GQ30.rigorous.noJ702.chr1-10.vcf.gz"
        fi

        echo "Removing J702-shared ${variant_type} from ${sample}"

        bcftools isec \
            -C \
            -c none \
            -w1 \
            "$input" \
            "$j702_comparison" \
            -Oz \
            -o "$output"

        bcftools index \
            --threads "$threads" \
            -f \
            -c \
            "$output"

        before_count=$(bcftools view -H "$input" | wc -l)
        after_count=$(bcftools view -H "$output" | wc -l)
        removed_count=$((before_count - after_count))

        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$sample" \
            "$variant_type" \
            "$before_count" \
            "$after_count" \
            "$removed_count" \
            >> "$j702_removal_summary"
    done
done

#######################################################################
# STEP 7: Merge the five MA lines
#######################################################################

echo "STEP 7: Merging MA-line VCFs"

for variant_type in snps indels
do
    list_file="${merged_dir}/${variant_type}.noJ702.list"
    merged="${merged_dir}/five_MA_lines.${variant_type}.noJ702.merged.vcf.gz"

    : > "$list_file"

    for sample in "${ma_samples[@]}"
    do
        if [[ "$variant_type" == "snps" ]]
        then
            input="${noj702_dir}/${sample}.deepvariant.norm.biallelic_snps.DP40-100.GQ30.rigorous.noJ702.chr1-10.vcf.gz"
        else
            input="${noj702_dir}/${sample}.deepvariant.norm.biallelic_indels.size1-50.DP40-100.GQ30.rigorous.noJ702.chr1-10.vcf.gz"
        fi

        [[ -s "$input" ]] || {
            echo "ERROR: Missing merge input: $input" >&2
            exit 1
        }

        printf "%s\n" "$input" >> "$list_file"
    done

    # Merge identical alleles across samples but do not create new
    # multiallelic records from different ALT alleles at the same position.
    bcftools merge \
        --threads "$threads" \
        -m none \
        -l "$list_file" \
        -Oz \
        -o "$merged"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$merged"
done

#######################################################################
# STEP 8: Select private variants
#######################################################################

echo "STEP 8: Selecting variants found in exactly one MA line"

private_snps="${private_dir}/private_biallelic_snps.exactly_one_line.chr1-10.vcf.gz"

private_indels="${private_dir}/private_biallelic_indels.size1-50.exactly_one_line.chr1-10.vcf.gz"

bcftools view \
    --threads "$threads" \
    -m2 \
    -M2 \
    -v snps \
    -i 'COUNT(GT="alt")=1' \
    -Oz \
    -o "$private_snps" \
    "${merged_dir}/five_MA_lines.snps.noJ702.merged.vcf.gz"

bcftools index \
    --threads "$threads" \
    -f \
    -c \
    "$private_snps"

bcftools view \
    --threads "$threads" \
    -m2 \
    -M2 \
    -v indels \
    -i '
        ABS(ILEN)>=1 &&
        ABS(ILEN)<=50 &&
        COUNT(GT="alt")=1
    ' \
    -Oz \
    -o "$private_indels" \
    "${merged_dir}/five_MA_lines.indels.noJ702.merged.vcf.gz"

bcftools index \
    --threads "$threads" \
    -f \
    -c \
    "$private_indels"

#######################################################################
# STEP 9: Produce individual VCFs and counts
#######################################################################

echo "STEP 9: Producing individual private-variant files"

count_summary="${summary_dir}/private_variant_counts.tsv"

printf \
    "Sample\tPrivate_SNPs\tSNP_hom_ALT\tSNP_het\tPrivate_indels_1_50bp\tIndel_hom_ALT\tIndel_het\n" \
    > "$count_summary"

for sample in "${ma_samples[@]}"
do
    sample_snps="${individual_dir}/${sample}.private_biallelic_snps.chr1-10.vcf.gz"

    sample_indels="${individual_dir}/${sample}.private_biallelic_indels.size1-50.chr1-10.vcf.gz"

    bcftools view \
        -s "$sample" \
        -Ou \
        "$private_snps" |
    bcftools view \
        -i 'GT="alt"' \
        -Oz \
        -o "$sample_snps"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$sample_snps"

    bcftools view \
        -s "$sample" \
        -Ou \
        "$private_indels" |
    bcftools view \
        -i 'GT="alt"' \
        -Oz \
        -o "$sample_indels"

    bcftools index \
        --threads "$threads" \
        -f \
        -c \
        "$sample_indels"

    snp_total=$(bcftools view -H "$sample_snps" | wc -l)
    snp_hom=$(bcftools view -H -i 'GT="AA"' "$sample_snps" | wc -l)
    snp_het=$(bcftools view -H -i 'GT="het"' "$sample_snps" | wc -l)

    indel_total=$(bcftools view -H "$sample_indels" | wc -l)
    indel_hom=$(bcftools view -H -i 'GT="AA"' "$sample_indels" | wc -l)
    indel_het=$(bcftools view -H -i 'GT="het"' "$sample_indels" | wc -l)

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$snp_total" \
        "$snp_hom" \
        "$snp_het" \
        "$indel_total" \
        "$indel_hom" \
        "$indel_het" \
        >> "$count_summary"
done

#######################################################################
# STEP 10: Create readable TSV variant lists
#######################################################################

echo "STEP 10: Creating TSV lists"

for sample in "${ma_samples[@]}"
do
    sample_snps="${individual_dir}/${sample}.private_biallelic_snps.chr1-10.vcf.gz"

    sample_indels="${individual_dir}/${sample}.private_biallelic_indels.size1-50.chr1-10.vcf.gz"

    snp_table="${individual_dir}/${sample}.private_biallelic_snps.chr1-10.tsv"

    indel_table="${individual_dir}/${sample}.private_biallelic_indels.size1-50.chr1-10.tsv"

    printf \
        "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tGT\tDP\tGQ\tAD\n" \
        > "$snp_table"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ\t%AD]\n' \
        "$sample_snps" \
        >> "$snp_table"

    printf \
        "CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tGT\tDP\tGQ\tAD\tSIGNED_ILEN\tABS_SIZE\tTYPE\n" \
        > "$indel_table"

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%GQ\t%AD]\t%ILEN\n' \
        "$sample_indels" |
    awk '
        BEGIN {
            OFS = "\t"
        }
        {
            ilen = $11
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

            print $0, size, type
        }
    ' >> "$indel_table"
done

#######################################################################
# STEP 11: Validate final variant types and sizes
#######################################################################

echo "STEP 11: Validating final results"

invalid_snps=$(
    bcftools view \
        -H \
        -e 'N_ALT=1 && TYPE="snp"' \
        "$private_snps" |
    wc -l
)

invalid_indels=$(
    bcftools view \
        -H \
        -e '
            N_ALT=1 &&
            TYPE="indel" &&
            ABS(ILEN)>=1 &&
            ABS(ILEN)<=50
        ' \
        "$private_indels" |
    wc -l
)

if [[ "$invalid_snps" -ne 0 ]]
then
    echo "ERROR: Found ${invalid_snps} invalid SNP records" >&2
    exit 1
fi

if [[ "$invalid_indels" -ne 0 ]]
then
    echo "ERROR: Found ${invalid_indels} invalid indel records" >&2
    exit 1
fi

echo
echo "Private-variant counts:"
column -t "$count_summary" || cat "$count_summary"

echo
echo "Workflow completed successfully"
date
