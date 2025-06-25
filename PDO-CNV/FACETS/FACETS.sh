# This script is used to run FACETS by using run-facets-wrapper.R
#!/bin/bash
#SBATCH --job-name=facets
#SBATCH --output=slurm_log/facets_%A_%a.out
#SBATCH --error=slurm_log/facets_%A_%a.err
#SBATCH --array=0-5
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4

cd /path

# Load R if needed (depends on your environment)
module load gcc/12.2.0 r/4.4.0

# Define input list
INPUT_DIR="Folder"
FILES=($(ls ${INPUT_DIR}/*.snp_pileup.gz))
COUNTS_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

# Sample ID from filename, stripping extension
BASENAME=$(basename "$COUNTS_FILE" .snp_pileup.gz)
SAMPLE_ID="${BASENAME//./_}"  # sanitize for output filenames

OUTDIR=./facets_out/$SAMPLE_ID
mkdir -pv $OUTDIR

export R_LIBS=/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.3

# Run wrapper
Rscript run-facets-wrapper.R \
    --counts-file "$COUNTS_FILE" \
    --sample-id "$SAMPLE_ID" \
    --purity-cval 1000 \
    --cval 500 \
    --everything \
    --facets-lib-path /ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.3 \
    --directory $OUTDIR
