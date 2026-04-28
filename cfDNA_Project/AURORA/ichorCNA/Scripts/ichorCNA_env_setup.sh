#!/usr/bin/env bash
set -euo pipefail

# Shared runtime bootstrap for ichorCNA Snakemake jobs on CRC.
# Keeps module stack + env vars deterministic across workflows.

module purge
module load r/4.5.0
module load hmmcopy_utils/dev
module load snakemake/8.25.2

# User R library path used by ichorCNA and dependencies.
export R_LIBS_USER="${R_LIBS_USER:-/ihome/alee/dal413/R/x86_64-pc-linux-gnu-library/4.5}"

# Prefer writable local temp dir across login/compute nodes.
export TMPDIR="${TMPDIR:-/tmp/ichorCNA_slurm-${SLURM_JOB_ID:-$$}}"
mkdir -p "$TMPDIR"

# Keep BLAS/OpenMP deterministic and avoid oversubscription.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

