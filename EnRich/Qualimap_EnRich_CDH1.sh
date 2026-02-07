#!/bin/bash
#SBATCH --job-name=qualimap_qc
#SBATCH --output=slurm_logs/qualimap/qualimap_qc_%A_%a.out
#SBATCH --error=slurm_logs/qualimap/qualimap_qc_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=2-9  # Adjust according to sample count (line 1 is header)

# Load required modules
module load qualimap/2.3
module load samtools/1.22.1

# FIX: Disable GUI requirements
unset DISPLAY
export _JAVA_OPTIONS="-Djava.awt.headless=true"

# Input paths
SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet.csv"
BED_FILE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Bed_file/03924250_Lee_CDH1_Targets_JT.bed.bed"
QUALIMAP_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Qualimap/20260105/qualimap_CDH1"

# Read current line from sample sheet
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SHEET")
SAMPLE=$(echo "$LINE" | cut -d',' -f1)
BAM_PATH=$(echo "$LINE" | cut -d',' -f2)

# Define centralized output folder
OUTPUT_DIR="${QUALIMAP_ROOT}/${SAMPLE}"
mkdir -p "$OUTPUT_DIR"

# Logging
echo "Running Qualimap QC for sample: ${SAMPLE}"
echo "BAM path: ${BAM_PATH}"
echo "Output directory: ${OUTPUT_DIR}"

# Optional: Check if BAM exists
if [ ! -f "$BAM_PATH" ]; then
    echo "[ERROR] BAM file not found: $BAM_PATH" >&2
    exit 1
fi

# Run Qualimap
qualimap bamqc \
    -bam "$BAM_PATH" \
    -outdir "$OUTPUT_DIR" \
    -nt ${SLURM_CPUS_PER_TASK} \
    --feature-file "$BED_FILE" \
    --java-mem-size=128G

# Report status
if [ $? -eq 0 ]; then
    echo "[SUCCESS] Qualimap QC completed for ${SAMPLE}"
else
    echo "[ERROR] Qualimap QC failed for ${SAMPLE}" >&2
fi
