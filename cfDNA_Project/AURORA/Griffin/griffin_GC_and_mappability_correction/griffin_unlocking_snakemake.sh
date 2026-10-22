#!/usr/bin/env bash

set -Eeuo pipefail

PROJECT_DIR=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_nucleosome_profiling
SNAKEMAKE=/ihome/alee/rak373/.conda/envs/griffin_updated/bin/snakemake

cd "${PROJECT_DIR}"

"${SNAKEMAKE}" \
  --snakefile "${PROJECT_DIR}/griffin_nucleosome_profiling.snakefile" \
  --configfile "${PROJECT_DIR}/config/samples.GC.yaml" \
  --configfile "${PROJECT_DIR}/config/sites.yaml" \
  --configfile "${PROJECT_DIR}/config/config.yaml" \
  --unlock