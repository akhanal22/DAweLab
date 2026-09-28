#!/bin/bash
#SBATCH --job-name=deepvariant_J702
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --array=1-8
#SBATCH --output=deepvariant_%A_%a.out
#SBATCH --error=deepvariant_%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL


module load SAMtools/1.23.1-GCC-13.3.0
# ============================================================
# USER SETTINGS
# ============================================================

# Reference genome
REFERENCE="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"

# Directory containing the original BAM files
BAM_DIR="."

# Sample list, one sample name per line
SAMPLE_LIST="allname.lst"

# DeepVariant output directory
OUTPUT_DIR="deepvariant_against_J702"

# DeepVariant Apptainer/Singularity image
DEEPVARIANT_IMAGE="/scratch/ak18787/software/deepvariant/deepvariant_1.10.0.sif"

# Number of CPUs
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# DeepVariant model
# Use PACBIO for PacBio HiFi reads
MODEL_TYPE="PACBIO"


# ============================================================
# LOAD SOFTWARE
# ============================================================


# Try to load Apptainer if it is available as a module.
# If Apptainer was installed locally, this section can be ignored.
if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_CMD="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_CMD="singularity"
else
    echo "ERROR: Neither apptainer nor singularity was found."
    exit 1
fi


# ============================================================
# READ SAMPLE NAME
# ============================================================

if [[ ! -s "$SAMPLE_LIST" ]]; then
    echo "ERROR: Sample list not found or empty: $SAMPLE_LIST"
    exit 1
fi

sample=$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" |
    tr -d '\r' |
    awk '{$1=$1};1'
)

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


# ============================================================
# DEFINE FILES
# ============================================================

original_bam="${BAM_DIR}/${sample}.sorted.primary_mapq20.bam"

prepared_dir="${OUTPUT_DIR}/prepared_bams"
vcf_dir="${OUTPUT_DIR}/vcf"
gvcf_dir="${OUTPUT_DIR}/gvcf"
log_dir="${OUTPUT_DIR}/logs"
tmp_dir="${OUTPUT_DIR}/tmp/${sample}"

rg_bam="${prepared_dir}/${sample}.sorted.primary_mapq20.RG.bam"

output_vcf="${vcf_dir}/${sample}.deepvariant.vcf.gz"
output_gvcf="${gvcf_dir}/${sample}.deepvariant.g.vcf.gz"

mkdir -p \
    "$prepared_dir" \
    "$vcf_dir" \
    "$gvcf_dir" \
    "$log_dir" \
    "$tmp_dir"


# ============================================================
# JOB INFORMATION
# ============================================================

echo "=========================================="
echo "Sample:       $sample"
echo "Host:         $(hostname)"
echo "Task:         ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:         $THREADS"
echo "Input BAM:    $original_bam"
echo "RG BAM:       $rg_bam"
echo "Reference:    $REFERENCE"
echo "Model:        $MODEL_TYPE"
echo "Container:    $DEEPVARIANT_IMAGE"
echo "Start time:   $(date)"
echo "=========================================="


# ============================================================
# CHECK INPUT FILES
# ============================================================

if [[ ! -s "$REFERENCE" ]]; then
    echo "ERROR: Reference genome not found or empty:"
    echo "$REFERENCE"
    exit 1
fi

if [[ ! -s "${REFERENCE}.fai" ]]; then
    echo "Reference index not found. Creating:"
    echo "${REFERENCE}.fai"

    samtools faidx "$REFERENCE"
fi

if [[ ! -s "$original_bam" ]]; then
    echo "ERROR: Original BAM not found or empty:"
    echo "$original_bam"
    exit 1
fi

if [[ ! -s "$DEEPVARIANT_IMAGE" ]]; then
    echo "ERROR: DeepVariant image not found:"
    echo "$DEEPVARIANT_IMAGE"
    exit 1
fi


# ============================================================
# CHECK ORIGINAL BAM
# ============================================================

echo
echo "Checking original BAM integrity..."

if ! samtools quickcheck -v "$original_bam"; then
    echo "ERROR: Original BAM failed samtools quickcheck:"
    echo "$original_bam"
    exit 1
fi

echo "Original BAM passed integrity check."


# ============================================================
# CREATE BAM WITH READ-GROUP TAGS
# ============================================================

create_rg_bam=false

if [[ ! -s "$rg_bam" ]]; then
    create_rg_bam=true
else
    echo
    echo "Existing RG-tagged BAM found. Checking it..."

    if ! samtools quickcheck -v "$rg_bam"; then
        echo "Existing RG BAM is invalid and will be recreated."
        create_rg_bam=true
    else
        # Check whether mapped alignments contain the correct RG tag
        rg_count=$(
            samtools view -F 4 "$rg_bam" |
            head -n 1000 |
            awk -v expected="RG:Z:${sample}" '
            {
                for (i=12; i<=NF; i++) {
                    if ($i == expected) {
                        count++
                        break
                    }
                }
            }
            END {
                print count+0
            }'
        )

        mapped_checked=$(
            samtools view -F 4 "$rg_bam" |
            head -n 1000 |
            wc -l
        )

        if [[ "$mapped_checked" -eq 0 ]]; then
            echo "ERROR: No mapped reads found in the existing RG BAM."
            exit 1
        fi

        if [[ "$rg_count" -ne "$mapped_checked" ]]; then
            echo "Existing BAM does not have the correct RG tag on every checked read."
            echo "Correct tags: $rg_count"
            echo "Reads checked: $mapped_checked"
            echo "The BAM will be recreated."
            create_rg_bam=true
        else
            echo "Existing RG BAM passed validation."
        fi
    fi
fi


if [[ "$create_rg_bam" == true ]]; then

    echo
    echo "Removing old RG BAM and index..."

    rm -f \
        "$rg_bam" \
        "${rg_bam}.bai" \
        "${rg_bam%.bam}.bai" \
        "${rg_bam}.tmp.bam" \
        "${rg_bam}.tmp.bam.bai"

    echo
    echo "Creating BAM with read-group information..."
    echo "Adding RG:Z:${sample} to every alignment..."

    samtools addreplacerg \
        -@ "$THREADS" \
        -m overwrite_all \
        -r "ID:${sample}" \
        -r "SM:${sample}" \
        -r "LB:${sample}" \
        -r "PL:PACBIO" \
        -r "PU:${sample}" \
        -o "${rg_bam}.tmp.bam" \
        "$original_bam"

    echo
    echo "Sorting RG-tagged BAM..."

    samtools sort \
        -@ "$THREADS" \
        -m 4G \
        -T "${tmp_dir}/${sample}.sort" \
        -o "$rg_bam" \
        "${rg_bam}.tmp.bam"

    rm -f "${rg_bam}.tmp.bam"

    echo
    echo "RG-tagged BAM created:"
    echo "$rg_bam"
fi


# ============================================================
# INDEX RG-TAGGED BAM
# ============================================================

echo
echo "Indexing RG-tagged BAM..."

rm -f \
    "${rg_bam}.bai" \
    "${rg_bam%.bam}.bai"

samtools index \
    -@ "$THREADS" \
    "$rg_bam"


# ============================================================
# VALIDATE RG-TAGGED BAM
# ============================================================

echo
echo "Checking RG-tagged BAM integrity..."

if ! samtools quickcheck -v "$rg_bam"; then
    echo "ERROR: RG-tagged BAM failed samtools quickcheck."
    exit 1
fi

if ! samtools idxstats "$rg_bam" >/dev/null; then
    echo "ERROR: BAM index could not be read."
    exit 1
fi


echo
echo "Read-group header:"

rg_header=$(
    samtools view -H "$rg_bam" |
    grep '^@RG' || true
)

if [[ -z "$rg_header" ]]; then
    echo "ERROR: No @RG line found in BAM header."
    exit 1
fi

echo "$rg_header"


echo
echo "Checking RG tags among the first 1,000 mapped alignments..."

mapped_checked=$(
    samtools view -F 4 "$rg_bam" |
    head -n 1000 |
    wc -l
)

if [[ "$mapped_checked" -eq 0 ]]; then
    echo "ERROR: No mapped alignments were found."
    exit 1
fi

correct_rg_count=$(
    samtools view -F 4 "$rg_bam" |
    head -n 1000 |
    awk -v expected="RG:Z:${sample}" '
    {
        found=0

        for (i=12; i<=NF; i++) {
            if ($i == expected) {
                found=1
                break
            }
        }

        if (found == 1) {
            count++
        }
    }
    END {
        print count+0
    }'
)

echo "Mapped alignments checked: $mapped_checked"
echo "Alignments with RG:Z:${sample}: $correct_rg_count"

if [[ "$correct_rg_count" -ne "$mapped_checked" ]]; then
    echo "ERROR: Some alignments do not contain RG:Z:${sample}."
    exit 1
fi

echo "All checked alignments contain the correct RG tag."


echo
echo "Example RG tags:"

samtools view -F 4 "$rg_bam" |
head -n 5 |
awk '{
    for (i=12; i<=NF; i++) {
        if ($i ~ /^RG:Z:/) {
            print $1, $i
            break
        }
    }
}'


# ============================================================
# PREPARE CONTAINER BIND PATHS
# ============================================================

reference_abs=$(readlink -f "$REFERENCE")
bam_abs=$(readlink -f "$rg_bam")
output_vcf_abs=$(readlink -m "$output_vcf")
output_gvcf_abs=$(readlink -m "$output_gvcf")
tmp_abs=$(readlink -m "$tmp_dir")

reference_dir=$(dirname "$reference_abs")
bam_dir=$(dirname "$bam_abs")
vcf_dir_abs=$(dirname "$output_vcf_abs")
gvcf_dir_abs=$(dirname "$output_gvcf_abs")

mkdir -p \
    "$vcf_dir_abs" \
    "$gvcf_dir_abs" \
    "$tmp_abs"


# ============================================================
# RUN DEEPVARIANT
# ============================================================

echo
echo "=========================================="
echo "Running DeepVariant"
echo "=========================================="

"$CONTAINER_CMD" exec \
    --bind "${reference_dir}:${reference_dir}" \
    --bind "${bam_dir}:${bam_dir}" \
    --bind "${vcf_dir_abs}:${vcf_dir_abs}" \
    --bind "${gvcf_dir_abs}:${gvcf_dir_abs}" \
    --bind "${tmp_abs}:${tmp_abs}" \
    "$DEEPVARIANT_IMAGE" \
    /opt/deepvariant/bin/run_deepvariant \
        --model_type="$MODEL_TYPE" \
        --ref="$reference_abs" \
        --reads="$bam_abs" \
        --output_vcf="$output_vcf_abs" \
        --output_gvcf="$output_gvcf_abs" \
        --num_shards="$THREADS" \
        --intermediate_results_dir="$tmp_abs"


# ============================================================
# CHECK DEEPVARIANT OUTPUT
# ============================================================

echo
echo "Checking DeepVariant output..."

if [[ ! -s "$output_vcf" ]]; then
    echo "ERROR: DeepVariant VCF was not created:"
    echo "$output_vcf"
    exit 1
fi

if [[ ! -s "$output_gvcf" ]]; then
    echo "ERROR: DeepVariant gVCF was not created:"
    echo "$output_gvcf"
    exit 1
fi

echo
echo "DeepVariant finished successfully."
echo "VCF:"
echo "$output_vcf"

echo
echo "gVCF:"
echo "$output_gvcf"

echo
echo "Prepared RG BAM:"
echo "$rg_bam"

echo
echo "End time: $(date)"
echo "=========================================="
