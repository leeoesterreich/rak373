#!/bin/bash
#SBATCH --job-name=ichorCNA_snakemake
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/ichorCNA/1000kb/ichorCNAsnake_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/ichorCNA/1000kb/ichorCNA_%j.err

set -euo pipefail

# Shared runtime setup (modules + R libs + TMPDIR + thread caps)
source /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Scripts/ichorCNA_env_setup.sh

# Override these when needed:
#   RUN_DIR=/path/to/run_dir
#   SNAKEMAKE_MAX_JOBS=1
RUN_DIR="${RUN_DIR:-/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/ichorCNA/1000kb}"

cd "$RUN_DIR"
echo "Snakemake version: $(snakemake --version)"
echo "RUN_DIR=$RUN_DIR"
echo "TMPDIR=$TMPDIR"
echo "R_LIBS_USER=$R_LIBS_USER"

snakemake -s "$RUN_DIR/ichorCNA.snakefile" \
  --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --jobs "${SNAKEMAKE_MAX_JOBS:-1}" \
  --rerun-incomplete \
  --latency-wait 60

