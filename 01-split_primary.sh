#!/bin/bash
#SBATCH --job-name=splits           # Job name
#SBATCH --partition=batch              # Partition name (batch, highmem_p, hugemem_p, or gpu_p)
#SBATCH --ntasks=1                      # 1 task (process) for below commands
#SBATCH --cpus-per-task=8               # CPU core count per task, by default 1 CPU core per task
#SBATCH --time=10:00:00                  # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --mem=200G                        # Memory per node (4GB); by default using M as unit
#SBATCH --output=%x.%j.out            # Standard output log
#SBATCH --error=%x.%j.err             # Standard error log
#SBATCH --mail-user=ak18787@uga.edu    # Where to send mail
#SBATCH --mail-type=END,FAIL            # Mail events (BEGIN, END, FAIL, ALL)


# ---- activate conda ----
source ~/ENTER/etc/profile.d/conda.sh
conda activate samtools

in="B73D1.primary.depth"   # change to your file
outdir="split_primary"
mkdir -p "$outdir"

awk -v outdir="$outdir" '
{
  # $1 = contig name
  fn = outdir "/" $1 ".primary.depth"
  print > fn
}
' "$in"
