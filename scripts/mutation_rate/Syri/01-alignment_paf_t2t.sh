#!/bin/bash
#SBATCH --job-name=paf_genome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=400G
#SBATCH --array=1-7
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu
#SBATCH --mail-type=END,FAIL

# --------------------
# Load environment
# --------------------
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools


# --------------------
# 1. Get sample name
# --------------------
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" allname.lst | tr -d '\r' | xargs)

if [[ -z "$file" ]]; then
  echo "ERROR: Empty sample name for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

echo "[$(date)] Processing sample: $file"

# --------------------
# 2. Input files
# --------------------
ref="J702_v3_modified_ragtag_re/flipped/J702_chr1-10.fa"
qry="${file}_v3_modified_ragtag_re_after_telomere_rearrangement/${file}_chr1-10.fa"


#outdir="paf_mapped_against_J702_modified_re_chrrename"
#outfile="${outdir}/B73D2_vs_J702_modified_v2.paf"

outdir="paf_genome_mapped_against_J702_modified_re_after_t2t"
outfile="${outdir}/${file}_vs_J702_modified.sorted.paf"

# --------------------
# 3. Sanity checks
# --------------------
[[ -s "$ref" ]] || { echo "ERROR: Missing reference $ref"; exit 1; }
[[ -s "$qry" ]] || { echo "ERROR: Missing query $qry"; exit 1; }

mkdir -p "$outdir"
# --------------------
# 4. Create samtools FASTA index (.fai)
# --------------------
#if [[ ! -s "${qry}.fai" ]]; then
 #   echo "[$(date)] Creating FASTA index (.fai)"
#    samtools faidx "$qry"
#fi

# --------------------
# 5. Create minimap2 index (.mmi)
# --------------------
#ref_index="${ref}.mmi"

#if [[ ! -s "$ref_index" ]]; then
 #   echo "[$(date)] Creating minimap2 index (.mmi)"
  #  minimap2 -d "$ref_index" "$ref"
#fi

# --------------------
# 6. Run minimap2 alignment
# --------------------
echo "[$(date)] Running minimap2..."

minimap2 \
  -cx asm5 \
  --eqx \
  --secondary=no \
  -t "${SLURM_CPUS_PER_TASK:-4}" \
  "$ref" "$qry" > "$outfile"

echo "[$(date)] Finished: $outfile"

