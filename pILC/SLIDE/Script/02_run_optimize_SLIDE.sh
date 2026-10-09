#!/bin/bash
#SBATCH --job-name=pILC_SLIDE_retry
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/SLIDE/Logs_pleomorphic/pILC_SLIDE_retry_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/SLIDE/Logs_pleomorphic/pILC_SLIDE_retry_%j.err

set -euo pipefail

module load r/4.6.0

SLIDE_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/SLIDE"

mkdir -p "${SLIDE_DIR}/Logs_pleomorphic"
mkdir -p "${SLIDE_DIR}/Results_pleomorphic_stable_retry"

cd "${SLIDE_DIR}"

which Rscript
Rscript --version

Rscript --vanilla "${SLIDE_DIR}/02_run_optimize_SLIDE.R"
