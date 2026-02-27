#!/bin/bash
#SBATCH --job-name=b73d2_binplot_png
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --array=1-1
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --mail-user=ak18787@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL            # Mail events (BEGIN, END, FAIL, ALL)

# ---- activate conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate base

# ---- params ----
BIN=100        # 50 kb bins
WINDOW=100    # 50 kb sliding window
STEP=100       # 50 kb step
MAXY=""         # optional, e.g. MAXY=120
DELETE_BINS=1   # 1 = delete bin files after plotting; 0 = keep bin files

# ---- sample for this array task ----
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" b73d2.lst)
outdir="./${file}"
contigs_txt="${outdir}/chr_contigs.txt"

if [[ ! -f "${contigs_txt}" ]]; then
  echo "ERROR: contigs file not found: ${contigs_txt}"
  exit 1
fi

if [[ ! -x "./plot_from_bins.py" ]]; then
  echo "ERROR: ./plot_from_bins.py not found or not executable"
  echo "Run: chmod +x plot_from_bins.py"
  exit 1
fi

echo "Sample: ${file}"
echo "BIN=${BIN} WINDOW=${WINDOW} STEP=${STEP} (PNG only)"

while IFS= read -r contig; do
  [[ -z "${contig}" ]] && continue
  safe_contig=$(echo "${contig}" | sed 's/[^A-Za-z0-9._-]/_/g')

  primary_depth="${outdir}/depth_primary/split_primary/${safe_contig}.primary.depth"
  secondary_depth="${outdir}/depth_secondary/split_secondary/${safe_contig}.secondary.depth"

  if [[ ! -s "${primary_depth}" ]]; then
    echo "WARN: missing/empty primary depth: ${primary_depth} (skip)"
    continue
  fi
  if [[ ! -s "${secondary_depth}" ]]; then
    echo "WARN: missing/empty secondary depth: ${secondary_depth} (skip)"
    continue
  fi

  # Temporary binned files (created to plot fast)
  primary_bins="${outdir}/${file}_${safe_contig}.bin${BIN}.primary.tsv"
  secondary_bins="${outdir}/${file}_${safe_contig}.bin${BIN}.secondary.tsv"

  # Final PNG (saved in MAIN folder)
  plot_png="${outdir}/${file}_${safe_contig}.w${WINDOW}_s${STEP}.depth.png"

  echo "Contig: ${contig}"

  # 1) Bin PRIMARY (only if not already done)
  if [[ ! -s "${primary_bins}" ]]; then
    awk -v BIN="${BIN}" '{
      b=int(($2-1)/BIN);
      sum[b]+=$3; cnt[b]++
    }
    END{
      for (b in sum){
        mid=(b*BIN + BIN/2)/1000;   # midpoint kbp
        print mid, sum[b]/cnt[b]
      }
    }' OFS="\t" "${primary_depth}" | sort -k1,1n > "${primary_bins}"
  fi

  # 2) Bin SECONDARY (only if not already done)
  if [[ ! -s "${secondary_bins}" ]]; then
    awk -v BIN="${BIN}" '{
      b=int(($2-1)/BIN);
      sum[b]+=$3; cnt[b]++
    }
    END{
      for (b in sum){
        mid=(b*BIN + BIN/2)/1000;   # midpoint kbp
        print mid, sum[b]/cnt[b]
      }
    }' OFS="\t" "${secondary_depth}" | sort -k1,1n > "${secondary_bins}"
  fi

  # 3) Plot from bins (NO TSV output)
  if [[ -n "${MAXY}" ]]; then
    ./plot_from_bins.py \
      --primary_bins "${primary_bins}" \
      --secondary_bins "${secondary_bins}" \
      --bin_bp "${BIN}" \
      --window_bp "${WINDOW}" \
      --step_bp "${STEP}" \
      --out_png "${plot_png}" \
      --title "${safe_contig}" \
      --maxy "${MAXY}"
  else
    ./plot_from_bins.py \
      --primary_bins "${primary_bins}" \
      --secondary_bins "${secondary_bins}" \
      --bin_bp "${BIN}" \
      --window_bp "${WINDOW}" \
      --step_bp "${STEP}" \
      --out_png "${plot_png}" \
      --title "${safe_contig}"
  fi

  # 4) Optional cleanup: keep only PNGs
  if [[ "${DELETE_BINS}" -eq 1 ]]; then
    rm -f "${primary_bins}" "${secondary_bins}"
  fi

done < "${contigs_txt}"

echo "Done. PNGs saved in: ${outdir}/"
