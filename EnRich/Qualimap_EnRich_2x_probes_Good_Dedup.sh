#!/bin/bash
#SBATCH --job-name=qualimap_qc
#SBATCH --output=slurm_logs/qualimap/qualimap_qc_%A_%a.out
#SBATCH --error=slurm_logs/qualimap/qualimap_qc_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=0-7   # 8 samples (0-based)

set -euo pipefail

# Load required modules
module load qualimap/2.3
module load samtools/1.22.1

# Disable GUI requirements
unset DISPLAY
export _JAVA_OPTIONS="-Djava.awt.headless=true"

# -------------------------
# Input paths
# -------------------------
DEDUP_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Dedup_010526"
BED_FILE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Bed_file/Lee_2x_probes_Good.bed.bed"
QUALIMAP_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Qualimap/20260105/qualimap_2x_probes_Good_Dedup"

mkdir -p "$QUALIMAP_ROOT" slurm_logs/qualimap

# -------------------------
# Build BAM list 
# -------------------------
mapfile -t BAM_LIST < <(find "$DEDUP_ROOT" -maxdepth 2 -type f -name "*.dedup.bam" | sort)

N=${#BAM_LIST[@]}
echo "[INFO] Found $N dedup BAM files under: $DEDUP_ROOT"

if [[ "$N" -eq 0 ]]; then
  echo "[ERROR] No *.dedup.bam found in $DEDUP_ROOT" >&2
  exit 1
fi

# Validate array bounds
TASK_ID="${SLURM_ARRAY_TASK_ID}"
if [[ "$TASK_ID" -lt 0 || "$TASK_ID" -ge "$N" ]]; then
  echo "[SKIP] Task $TASK_ID out of range (0..$((N-1))). Found $N samples." >&2
  exit 0
fi

BAM_PATH="${BAM_LIST[$TASK_ID]}"
SAMPLE="$(basename "$(dirname "$BAM_PATH")")"

# If you prefer sample name from BAM filename:
# SAMPLE="$(basename "$BAM_PATH" .dedup.bam)"

OUTPUT_DIR="${QUALIMAP_ROOT}/${SAMPLE}"
mkdir -p "$OUTPUT_DIR"

# -------------------------
# Run Qualimap
# -------------------------
qualimap bamqc \
  -bam "$BAM_PATH" \
  -outdir "$OUTPUT_DIR" \
  -nt "${SLURM_CPUS_PER_TASK}" \
  --feature-file "$BED_FILE" \
  --java-mem-size=128G

echo "[SUCCESS] Qualimap QC completed for ${SAMPLE}"
