#!/bin/bash
#
#SBATCH --job-name=ichorCNAsnakemake
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # Request that ncpus be allocated per process.
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/ichorCNA/ichorCNAsnake_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/ichorCNA/ichorCNA_%j.err
set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/ichorCNA_SOP}"
RUN_DIR="${RUN_DIR:-$SUBMIT_DIR}"
source "${SUBMIT_DIR}/Scripts/ichorCNA_env_setup.sh"

cd "${RUN_DIR}"
echo "Snakemake version: $(snakemake --version)"

snakemake -s "${RUN_DIR}/ichorCNA.snakefile" \
  --cores "${SLURM_CPUS_PER_TASK:-16}" \
  --jobs "${SNAKEMAKE_MAX_JOBS:-1}" \
  --rerun-incomplete \
  --latency-wait 60