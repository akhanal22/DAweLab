This is the code used for calling snvs using long read pacbio hifi data, I already had the bam file that was made using minimap2 by mapping raw reads against the J702 modified reference genomes.
So I used longcall (Version: 0.0.8) to run following script 
```
#!/bin/bash
#SBATCH --job-name=longcall
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --array=1-8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source ~/ENTER/etc/profile.d/conda.sh
conda activate base

file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" linename.lst | tr -d '\r' | xargs)

ref="../J702_v3_modified_ragtag_re/ragtag.scaffold.fasta"
bam="${file}.sorted.bam"
outdir="vcf_against_j702_snvs_longcall"

mkdir -p "$outdir"

echo "Sample: $file"
echo "Reference: $ref"
echo "BAM: $bam"
echo "Threads: $SLURM_CPUS_PER_TASK"
date

[[ -s "$ref" ]] || { echo "ERROR: missing reference $ref"; exit 1; }
[[ -s "$bam" ]] || { echo "ERROR: missing BAM $bam"; exit 1; }

if [[ ! -s "${bam}.bai" ]]; then
    samtools index "$bam"
fi

longcallD call \
  -t ${SLURM_CPUS_PER_TASK} \
  "$ref" \
  "$bam" \
  > "${outdir}/${file}.longcallD.vcf"
  ```
  After this i only select for snps
  
  ```
  bcftools view -v snps B73D1.longcallD.vcf -Oz -o B73D1.longcallD.snps.vcf
```
  Now, the rigourous filtering only for confident snps
  ```
   bcftools filter   -e 'QUAL<50 || FILTER!="PASS" || FORMAT/DP<40 || FORMAT/GQ<30'   B73D1.longcallD.snps.vcf -Ou | bcftools view   -i '(
        GT="1|1" && FORMAT/AD[0:1]/FORMAT/DP >= 0.90
      ) || (
        (GT="0|1" || GT="1|0") &&
        FORMAT/AD[0:1]/FORMAT/DP >= 0.30 &&
        FORMAT/AD[0:1]/FORMAT/DP <= 0.70
      )'   -Oz -o B73D1.longcallD.snps.rigorous.vcf

  ```

Same filtering was applied for all other individuals
Now to be sure that the variant are accurate and not generated due to sequencing error, i filtered out the snvs that were present in J702 using following code:
```
 While read file; do bcftools isec   -n=1 -w1   ${file}.longcallD.snps.rigorous.vcf.gz   J702.longcallD.snps.vcf.gz   -Oz -o ${file}.noJ702.snps.vcf.gz ; done < linename.lst
``
