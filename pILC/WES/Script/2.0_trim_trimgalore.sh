#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --output=slurm_log/trimgalore/trimgalore_%A_%a.out
#SBATCH --error=slurm_log/trimgalore/trimgalore_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-47

set -euo pipefail
shopt -s nullglob

# Load module
module purge
module load trimgalore/0.6.10

# Paths
BASE_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20250611_HTGC_WES_pILC_Jagmohan"
OUTPUT_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/2_TrimGalore"

mkdir -p "$OUTPUT_BASE"
mkdir -p slurm_log/trimgalore

# Get sample folders
mapfile -t SAMPLE_DIRS < <(find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

SAMPLE_DIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "${SAMPLE_DIR}" || ! -d "${SAMPLE_DIR}" ]]; then
    echo "[ERROR] No valid sample directory found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Find lane-split FASTQs
mapfile -t R1_FILES < <(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" | sort)
mapfile -t R2_FILES < <(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*_R2_001.fastq.gz" | sort)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No R1 FASTQ files found in $SAMPLE_DIR"
    exit 1
fi

if [[ ${#R2_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No R2 FASTQ files found in $SAMPLE_DIR"
    exit 1
fi

if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
    echo "[ERROR] Mismatch in number of R1 and R2 files in $SAMPLE_DIR"
    echo "R1 count: ${#R1_FILES[@]}"
    echo "R2 count: ${#R2_FILES[@]}"
    exit 1
fi

# Derive sample name from FASTQ prefix
SAMPLE_NAME=$(basename "${R1_FILES[0]}" | sed -E 's/_OES.*$//')

echo "[INFO] Processing sample: $SAMPLE_NAME"
echo "[INFO] Sample folder: $SAMPLE_DIR"

echo "[INFO] R1 files:"
printf '  %s\n' "${R1_FILES[@]}"

echo "[INFO] R2 files:"
printf '  %s\n' "${R2_FILES[@]}"

# Output and temp working dir
OUT_DIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
mkdir -p "$OUT_DIR"

WORK_DIR="${SLURM_TMPDIR:-${OUT_DIR}/tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$WORK_DIR"

MERGED_R1="${WORK_DIR}/${SAMPLE_NAME}_R1.merged.fastq.gz"
MERGED_R2="${WORK_DIR}/${SAMPLE_NAME}_R2.merged.fastq.gz"

echo "[INFO] Merging lane FASTQs..."
cat "${R1_FILES[@]}" > "$MERGED_R1"
cat "${R2_FILES[@]}" > "$MERGED_R2"

# Run Trim Galore, use Illumina instead of the default auto-detection of adapter sequence.
trim_galore \
    --paired \
    --illumina \
    --gzip \
    --cores 8 \
    -o "$OUT_DIR" \
    "$MERGED_R1" "$MERGED_R2"

# Cleanup merged temporary files
rm -f "$MERGED_R1" "$MERGED_R2"

echo "[SUCCESS] Trim Galore completed for $SAMPLE_NAME"