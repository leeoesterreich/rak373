#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --output=slurm_log/trimgalore/trimgalore_%A_%a.out
#SBATCH --error=slurm_log/trimgalore/trimgalore_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1

# Load  module
module purge
module load trimgalore/0.6.10

# Paths
SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/2_TrimGalore/Sample_sheet.csv"
OUTPUT_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/2_TrimGalore"
mkdir -p "$OUTPUT_BASE" "slurm_log/trimgalore"

# Read line from CSV (skip header using tail +2)
LINE=$(tail -n +2 "$SHEET" | sed -n "${SLURM_ARRAY_TASK_ID}p" | tr -d '\r')
IFS=',' read -r SAMPLE_NAME FQ1 FQ2 <<< "$LINE"

echo "[INFO] Processing sample: $SAMPLE_NAME"
echo "  R1: $FQ1"
echo "  R2: $FQ2"

# Check if input files exist
if [[ ! -f "$FQ1" ]]; then
    echo "[ERROR] R1 file missing: $FQ1"
    exit 1
fi
if [[ ! -f "$FQ2" ]]; then
    echo "[ERROR] R2 file missing: $FQ2"
    exit 1
fi

# Create sample-specific output dir
OUT_DIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
mkdir -p "$OUT_DIR"

# Run Trim Galore
trim_galore --paired --cores 4 --gzip -o "$OUT_DIR" "$FQ1" "$FQ2"

# Report completion
if [[ $? -eq 0 ]]; then
    echo "[SUCCESS] Trim Galore completed for $SAMPLE_NAME"
else
    echo "[FAILURE] Trim Galore failed for $SAMPLE_NAME"
    exit 1
fi
