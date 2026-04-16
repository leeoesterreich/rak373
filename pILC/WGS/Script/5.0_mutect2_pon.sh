#!/bin/bash

# This script is used to run mutect2 per-normal to generate VCF file

#SBATCH --job-name=mutect2_pon
#SBATCH --output=slurm_log/mutect2/pon/mutect2_pon_%A_%a.out
#SBATCH --error=slurm_log/mutect2/pon/mutect2_pon_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-9  # Make sure this matches your sample count

module purge
module load singularity/4.3.2

# ==== CONFIG ====
NORMAL_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/5_PON/Normal_Sample_sheet.csv"
VCF_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/5_PON/vcfs"
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$VCF_DIR" "slurm_log/mutect2/pon"

# ==== READ CSV LINE ====
# Ensure file uses Unix line endings (in-memory fix)
LINE=$(tail -n +2 "$NORMAL_SHEET" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" | tr -d '\r')

# Debug: Check if line was read
echo "[DEBUG] Raw CSV line: '$LINE'"

# Parse CSV line: sample_name,BQSR_path
IFS=',' read -r SAMPLE BQSR_PATH <<< "$LINE"

# Debug: Check parsed values
echo "[DEBUG] Parsed SAMPLE: '$SAMPLE'"
echo "[DEBUG] Parsed BQSR_PATH: '$BQSR_PATH'"

# Validate
if [[ -z "$SAMPLE" || -z "$BQSR_PATH" ]]; then
  echo "[ERROR] Failed to parse line or sample is missing. Exiting."
  exit 1
fi

if [[ ! -f "$BQSR_PATH" ]]; then
  echo "[ERROR] BAM not found: $BQSR_PATH"
  exit 1
fi

# ==== OUTPUT PATHS ====
OUT_DIR="${VCF_DIR}/${SAMPLE}"
OUT_VCF="${OUT_DIR}/${SAMPLE}.vcf.gz"

mkdir -p "$OUT_DIR"

if [[ -f "$OUT_VCF" ]]; then
  echo "[SKIP] VCF already exists for $SAMPLE: $OUT_VCF"
  exit 0
fi

# ==== RUN MUTECT2 ====
echo "[INFO] Running Mutect2 for $SAMPLE"

singularity exec $BIND $SIF \
  gatk Mutect2 \
    -R "$REF" \
    -I "$BQSR_PATH" \
    -max-mnp-distance 0 \
    -O "$OUT_VCF"

# ==== CHECK RESULT ====
if [[ -f "$OUT_VCF" ]]; then
  echo "[DONE] VCF generated for $SAMPLE: $OUT_VCF"
  exit 0
else
  echo "[ERROR] Mutect2 failed for $SAMPLE"
  exit 1
fi
