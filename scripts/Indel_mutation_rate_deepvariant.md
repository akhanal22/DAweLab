This is the code used for calling snvs using long read pacbio hifi data, I already had the bam file that was made using minimap2 by mapping raw reads against the 
J702 modified reference genomes. So I used deepvariant (Version: 0.0.8) to run following script. But before that i selected for primary reads and that passed 
mapping quality of 20.
Step 1.
```
#!/bin/bash
#SBATCH --job-name=array_minimap
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# 1. Pull the name directly from your list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" all.lst)

# ---- 1) activate conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

# 2. I want to map the raw hifi reads against the maie v5 reference genome for each lines, converted it to bam files and sorted it.
# -t 8: Uses the 8 CPUs you requested
# -o: Creates a unique folder for this specific file to avoid conflicts

#minimap2 -t 8 -ax map-hifi ../J702_v3_modified_ragtag_re/J702_flip.fa ../fastX/${file}.hifireads.fastq.gz \
#| samtools view -@ 8 -b \
#| samtools sort -@ 8 -o ./${file}.sorted.bam


samtools view \
    -@ 8 \
    -b \
    -q 20 \
    -F 2308 \
    ${file}.sorted.bam \
> ${file}.sorted.primary_mapq20.bam


```
Step 2. Now lets run deepvariant
```
#!/bin/bash
#SBATCH --job-name=deepvariant
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ==================================================
# 1. Activate environment
# ==================================================

module load SAMtools/1.23.1-GCC-13.3.0

threads="${SLURM_CPUS_PER_TASK}"

# ==================================================
# 2. Read sample name from all.lst
# ==================================================

sample=$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" all.lst |
    tr -d '\r' |
    xargs
)

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=========================================="
echo "Sample: $sample"
echo "Host:   $(hostname)"
echo "Task:   ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:   $threads"
echo "=========================================="

# ==================================================
# 3. Define paths
# ==================================================

# Reference genome
reference="../J702_v3_modified_ragtag_re/J702_flip.fa"

# Original BAM without read groups
input_bam="${sample}.sorted.primary_mapq20.bam"

# DeepVariant image
# Change this if your SIF is stored somewhere else
deepvariant_image="/scratch/ak18787/software/deepvariant/deepvariant_1.10.0.sif"

# Output directories
output_dir="deepvariant_against_J702"
prepared_bam_dir="${output_dir}/prepared_bams"
log_dir="${output_dir}/logs/${sample}"

mkdir -p "$output_dir"
mkdir -p "$prepared_bam_dir"
mkdir -p "$log_dir"

# BAM containing read groups
rg_bam="${prepared_bam_dir}/${sample}.sorted.primary_mapq20.RG.bam"
rg_bai="${rg_bam}.bai"

# DeepVariant outputs
output_vcf="${output_dir}/${sample}.deepvariant.vcf.gz"
output_gvcf="${output_dir}/${sample}.deepvariant.g.vcf.gz"

# ==================================================
# 4. Check input files
# ==================================================

if [[ ! -s "$reference" ]]; then
    echo "ERROR: Reference genome not found:"
    echo "$reference"
    exit 1
fi

if [[ ! -s "${reference}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found:"
    echo "${reference}.fai"
    exit 1
fi

if [[ ! -s "$input_bam" ]]; then
    echo "ERROR: Input BAM not found:"
    echo "$input_bam"
    exit 1
fi

if [[ ! -s "$deepvariant_image" ]]; then
    echo "ERROR: DeepVariant image not found:"
    echo "$deepvariant_image"
    exit 1
fi

# ==================================================
# 5. Check original BAM
# ==================================================

echo "Checking original BAM integrity..."

if ! samtools quickcheck -v "$input_bam"; then
    echo "ERROR: Input BAM failed samtools quickcheck:"
    echo "$input_bam"
    exit 1
fi

# Check whether original BAM has an index
if [[ ! -s "${input_bam}.bai" ]]; then
    echo "Input BAM index not found. Creating it..."

    samtools index \
        -@ "$threads" \
        "$input_bam"
fi

if ! samtools idxstats "$input_bam" >/dev/null 2>&1; then
    echo "ERROR: Input BAM index cannot be read:"
    echo "${input_bam}.bai"
    exit 1
fi

# ==================================================
# 6. Add read groups
# ==================================================

# Reuse the prepared BAM if it already exists and is valid.
if [[ -s "$rg_bam" ]] &&
   samtools quickcheck "$rg_bam" &&
   samtools view -H "$rg_bam" | grep -q '^@RG'; then

    echo "Valid RG-tagged BAM already exists."
    echo "Reusing: $rg_bam"

else

    echo "Creating BAM with read-group information..."

    # Remove incomplete files left by an earlier failed run
    rm -f "$rg_bam" "$rg_bai"

    samtools addreplacerg \
        -@ "$threads" \
        -r "ID:${sample}" \
        -r "SM:${sample}" \
        -r "LB:${sample}" \
        -r "PL:PACBIO" \
        -r "PU:${sample}" \
        -o "$rg_bam" \
        "$input_bam"

    echo "RG-tagged BAM created:"
    echo "$rg_bam"

fi

# ==================================================
# 7. Index and verify RG BAM
# ==================================================

if [[ ! -s "$rg_bai" ]]; then
    echo "Indexing RG-tagged BAM..."

    samtools index \
        -@ "$threads" \
        "$rg_bam"
fi

echo "Checking RG-tagged BAM integrity..."

if ! samtools quickcheck -v "$rg_bam"; then
    echo "ERROR: RG-tagged BAM failed samtools quickcheck:"
    echo "$rg_bam"
    exit 1
fi

if ! samtools idxstats "$rg_bam" >/dev/null 2>&1; then
    echo "ERROR: RG-tagged BAM index cannot be read:"
    echo "$rg_bai"
    exit 1
fi

echo "Read-group header:"

if ! samtools view -H "$rg_bam" | grep '^@RG'; then
    echo "ERROR: @RG header is missing from prepared BAM."
    exit 1
fi

# Confirm that at least one alignment contains an RG tag
if ! samtools view "$rg_bam" |
    awk '
        NR <= 1000 {
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^RG:Z:/) {
                    found = 1
                    exit
                }
            }
        }
        END {
            exit(found ? 0 : 1)
        }
    '
then
    echo "ERROR: No RG tag found among the first 1,000 alignments."
    exit 1
fi

echo "RG tags are present on the alignments."

# ==================================================
# 8. Convert paths to absolute paths
# ==================================================

reference_abs=$(readlink -f "$reference")
rg_bam_abs=$(readlink -f "$rg_bam")
output_dir_abs=$(readlink -f "$output_dir")
image_abs=$(readlink -f "$deepvariant_image")

reference_dir=$(dirname "$reference_abs")
reference_name=$(basename "$reference_abs")

bam_dir=$(dirname "$rg_bam_abs")
bam_name=$(basename "$rg_bam_abs")

echo "Reference: $reference_abs"
echo "BAM:       $rg_bam_abs"
echo "Image:     $image_abs"
echo "Output:    $output_dir_abs"

# ==================================================
# 9. Find the container program
# ==================================================

if command -v apptainer >/dev/null 2>&1; then
    container_command="apptainer"

elif command -v singularity >/dev/null 2>&1; then
    container_command="singularity"

else
    echo "ERROR: Neither Apptainer nor Singularity is available."
    echo "This script must run on a Sapelo2 compute node."
    exit 1
fi

echo "Container command: $container_command"

# ==================================================
# 10. Run DeepVariant
# ==================================================

echo "Starting DeepVariant for $sample..."

"$container_command" exec \
    --bind "${reference_dir}:/reference" \
    --bind "${bam_dir}:/reads" \
    --bind "${output_dir_abs}:/output" \
    "$image_abs" \
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref="/reference/${reference_name}" \
        --reads="/reads/${bam_name}" \
        --output_vcf="/output/${sample}.deepvariant.vcf.gz" \
        --output_gvcf="/output/${sample}.deepvariant.g.vcf.gz" \
        --num_shards="$threads" \
        --logging_dir="/output/logs/${sample}" \
        --vcf_stats_report=true

# ==================================================
# 11. Verify DeepVariant output
# ==================================================

if [[ ! -s "$output_vcf" ]]; then
    echo "ERROR: DeepVariant VCF was not generated:"
    echo "$output_vcf"
    exit 1
fi

if [[ ! -s "$output_gvcf" ]]; then
    echo "ERROR: DeepVariant gVCF was not generated:"
    echo "$output_gvcf"
    exit 1
fi

# ==================================================
# 12. Index VCF and gVCF
# ==================================================

echo "Indexing DeepVariant outputs..."

bcftools index \
    --threads "$threads" \
    --tbi \
    --force \
    "$output_vcf"

bcftools index \
    --threads "$threads" \
    --tbi \
    --force \
    "$output_gvcf"

# ==================================================
# 13. Count output variants
# ==================================================

total_records=$(
    bcftools view -H "$output_vcf" |
    wc -l
)

snp_records=$(
    bcftools view -H -v snps "$output_vcf" |
    wc -l
)

indel_records=$(
    bcftools view -H -v indels "$output_vcf" |
    wc -l
)

echo "=========================================="
echo "DeepVariant completed successfully"
echo "Sample:            $sample"
echo "Prepared BAM:      $rg_bam"
echo "VCF:               $output_vcf"
echo "gVCF:              $output_gvcf"
echo "Total VCF records: $total_records"
echo "SNP records:       $snp_records"
echo "Indel records:     $indel_records"
echo "=========================================="
```
step 3. Select fr indels only
```
for sample in B73D1 B73D2 B73D4 B73D7B B73D7c B73D10
do
    echo "Processing $sample ..."

    bcftools view \
        -v indels \
        -f PASS \
        ${sample}.deepvariant.vcf.gz \
        -Oz \
        -o indel_pipeline/${sample}.deepvariant.pass.indels.vcf.gz

    bcftools index -t indel_pipeline/${sample}.deepvariant.pass.indels.vcf.gz

    echo "Finished $sample"
done
```
step 4. now lets normalize each vcf according to the reference genome 

```
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

source ~/ENTER/etc/profile.d/conda.sh
conda activate base

# reference used for mapping/calling
ref="/scratch/ak18787/B73_project/J702_v3_modified_ragtag_re/J702_flip.fa"

# sample name from list
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst | tr -d '\r' | xargs)


gzvcf="${file}.deepvariant.pass.indels.vcf.gz"
normvcf="${file}.deepvariant.pass.indels.norm.vcf.gz"

echo "Sample: $file"
#echo "Input VCF: $invcf"
echo "Reference: $ref"
date

# checks
[[ -s "$ref" ]] || { echo "ERROR: missing reference $ref"; exit 1; }
[[ -s "${ref}.fai" ]] || samtools faidx "$ref"

# compress if needed
#if [[ ! -s "$gzvcf" ]]; then
 #   [[ -s "$invcf" ]] || { echo "ERROR: missing input VCF $invcf"; exit 1; }
  #  bgzip -c "$invcf" > "$gzvcf"
#fi

# index input
if [[ ! -s "${gzvcf}.csi" ]]; then
    bcftools index -c "$gzvcf"
fi

# normalize
bcftools norm \
  -f "$ref" \
  -Oz \
  -o "$normvcf" \
  "$gzvcf"

# index normalized output
bcftools index -c "$normvcf"

echo "Finished normalization for $file"
date
```
Step 5. filtering for snps od depth <40 , GQ <30, QUAL <50 and homozygous reads to have coverage > 90% and heterozygous to have read support of atleast 30-70%. 
```
samples="B73D1 B73D2 B73D4 B73D7c B73D10 J702"

for sample in $samples
do
    echo "Processing ${sample}..."

    bcftools filter \
        -e 'QUAL<50 || FILTER!="PASS" || FORMAT/DP<40 || FORMAT/GQ<30' \
        ${sample}.deepvariant.pass.indels.norm.vcf.gz \
        -Ou | \
    bcftools view \
        -i '(
            GT="1/1" &&
            FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.90
        ) || (
            GT="0/1" &&
            FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.30 &&
            FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) <= 0.70
        )' \
        -Oz \
        -o ${sample}.deepvariant.indels.norm.rigorous.vcf.gz

    bcftools index -t ${sample}.deepvariant.indels.norm.rigorous.vcf.gz

    echo -n "${sample} retained indels: "
    bcftools view -H ${sample}.deepvariant.indels.norm.rigorous.vcf.gz | wc -l

done
```
Step 6. Remove the indels that is shared with J702
```
while read file; do
    invcf="${file}.deepvariant.indels.norm.rigorous.vcf.gz"
    outvcf="${file}.deepvariant.indels.rigorous.norm.noJ702.vcf.gz"

    echo "Removing J702-shared SNPs from $file"
    echo "Input:  $invcf"
    echo "Output: $outvcf"

    bcftools isec \
      -n=1 -w1 \
      "$invcf" \
      J702.ldeepvariant.indels.norm.rigorous.vcf.gz \
      -Oz -o "$outvcf"

    bcftools index -c "$outvcf"

done < allname.lst
```
step 7. Merge all and select the unique indels only
```
module load BCFtools/1.21-GCC-13.3.0
bcftools merge   -l vcf.list   -Oz   -o five_lines.merged.vcf.gz # vcflist is the list of vcf files after rigorous filtering
bcftools index -c seven_lines.merged.vcf.gz
bcftools view   -m2 -M2 -v indels   -i 'COUNT(GT="alt")=1'   -Oz   -o private_indels.only_one_line.vcf.gz   seven_lines.merged.vcf.gz
bcftools index -c private_indels.only_one_line.vcf.gz
 awk '
BEGIN{OFS="\t"}
/^#CHROM/{
  for(i=10;i<=NF;i++) sample[i]=$i
  next
}
!/^#/{
  for(i=10;i<=NF;i++){
    gt=$i
    split(gt,a,":")
    if(a[1]=="0/1" || a[1]=="1/0" || a[1]=="1/1" || a[1]=="0|1" || a[1]=="1|0" || a[1]=="1|1"){
      count[sample[i]]++
    }
  }
}
END{
  for(s in count) print s, count[s]
}' private_indels.only_one_line.vcf > private_indels_counts.tsv
