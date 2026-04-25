#!/bin/bash
#SBATCH --job-name=facets
#SBATCH --output=slurm_log/facets_%A_%a.out
#SBATCH --error=slurm_log/facets_%A_%a.err
#SBATCH --array=0-15
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Load R if needed 
module load gcc/12.2.0 r/4.4.0

# Define input list
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/snp-pileup"
FILES=($(ls ${INPUT_DIR}/*.snp_pileup.gz))
COUNTS_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
run_facets="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Script/run-facets-wrapper.R"

# Sample ID from filename, stripping extension
BASENAME=$(basename "$COUNTS_FILE" .snp_pileup.gz)
SAMPLE_ID="${BASENAME//./_}"  # sanitize for output filenames

OUTDIR=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Result/$SAMPLE_ID
mkdir -pv $OUTDIR

export R_LIBS=/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4

# Run wrapper
Rscript $run_facets \
    --counts-file "$COUNTS_FILE" \
    --sample-id "$SAMPLE_ID" \
    --purity-cval 150 \
    --cval 150 \
    --everything \
    --facets-lib-path /ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4 \
    --directory $OUTDIR

