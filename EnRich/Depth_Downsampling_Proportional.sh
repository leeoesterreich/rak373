#!/bin/bash
# Fraction-wise Downsampling
#SBATCH --job-name=RAW_Downsample_0.5
#SBATCH --output=slurm_log/Downsample/RAW_DS05_%A_%a.out
#SBATCH --error=slurm_log/Downsample/RAW_DS05_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-7

set -euo pipefail

module purge
module load samtools/1.22.1

# =========================
# Paths
# =========================
SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet_010526.csv"
OUTDIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Downsampling/0.5X/010526_RAW_DS05"
mkdir -p "$OUTDIR" slurm_log/Downsample

# =========================
# Fixed downsampling params
# =========================
SEED=42
FRAC=0.5
SARG="${SEED}.5"   # samtools format: seed.fraction

# =========================
# Parse sample sheet
# =========================
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$SHEET" | tr -d '\r')
[[ -z "$LINE" ]] && exit 0

SAMPLE=$(echo "$LINE" | cut -d',' -f1 | xargs)
INPUT_BAM=$(echo "$LINE" | cut -d',' -f2 | xargs)

if [[ ! -f "$INPUT_BAM" ]]; then
  echo "[ERROR] BAM not found: $INPUT_BAM"
  exit 1
fi

SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
mkdir -p "$SAMPLE_DIR"

# =========================
# Downsampling
# =========================
OUT_BAM="${SAMPLE_DIR}/${SAMPLE}.ds0.5.sorted.bam"

echo "[INFO] Sample: $SAMPLE"
echo "[INFO] Downsampling fraction: $FRAC"
echo "[INFO] Seed argument: $SARG"

samtools view \
  -@ "$SLURM_CPUS_PER_TASK" \
  -b -s "$SARG" \
  "$INPUT_BAM" > "$OUT_BAM"

samtools index "$OUT_BAM"

echo "[DONE] 0.5x downsampling completed for $SAMPLE"
