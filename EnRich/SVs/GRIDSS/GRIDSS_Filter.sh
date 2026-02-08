#!/bin/bash

# This script is used to filter GRIDSS output vcf file. 
# High-confidence filter:
#   - FILTER == PASS
#   - SR >= 1
#   - RP + SR >= 8
#   - BQ >= 30

#SBATCH --job-name=GRIDSS_filter
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --array=1-8
#SBATCH --output=slurm_logs/GRIDSS/GRIDSS_filter_%A_%a.out
#SBATCH --error=slurm_logs/GRIDSS/GRIDSS_filter_%A_%a.err

set -euo pipefail

# ============================================================
# Modules
# ============================================================
module purge
module load bcftools/1.22
module load htslib/1.22

# ============================================================
# Paths
# ============================================================
base_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/GRIDSS_012425"

mkdir -p slurm_logs/GRIDSS

# ============================================================
# Resolve sample directory from array index
# ============================================================
mapfile -t sample_dirs < <(find "$base_dir" -mindepth 1 -maxdepth 1 -type d | sort)

idx=$((SLURM_ARRAY_TASK_ID - 1))
if [[ "$idx" -lt 0 || "$idx" -ge "${#sample_dirs[@]}" ]]; then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID out of range (n=${#sample_dirs[@]})"
  exit 1
fi

sample_dir="${sample_dirs[$idx]}"
sample="$(basename "$sample_dir")"

echo "============================================================"
echo "[INFO] SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "[INFO] Processing sample:   $sample"
echo "[INFO] Sample directory:    $sample_dir"
echo "============================================================"

in_vcf="${sample_dir}/${sample}.gridss.vcf.gz"
out_vcf="${sample_dir}/${sample}.gridss.filtered.vcf.gz"

# ============================================================
# Checks
# ============================================================
if [[ ! -s "$in_vcf" ]]; then
  echo "[ERROR] Input VCF not found or empty: $in_vcf"
  exit 1
fi

if [[ -s "$out_vcf" ]]; then
  echo "[SKIP] Filtered VCF already exists: $out_vcf"
  exit 0
fi

if [[ ! -f "${in_vcf}.tbi" ]]; then
  echo "[INFO] Indexing input VCF: $in_vcf"
  bcftools index -t "$in_vcf" || true
fi

# ============================================================
# Filtering
# NOTE:
# - GRIDSS usually stores SR/RP/BQ in FORMAT (per-sample), not INFO.
# - We filter on the sample FORMAT fields.
# - Also require a non-missing genotype call (GT not "./.") to avoid empty sample fields.
# ============================================================
echo "[INFO] Filtering GRIDSS VCF with rules:"
echo "       FILTER=PASS, FORMAT/SR>=1, (FORMAT/RP+FORMAT/SR)>=8, FORMAT/BQ>=30"

bcftools view \
  -i 'FILTER="PASS" && FORMAT/GT!="./." && FORMAT/SR>=1 && (FORMAT/RP + FORMAT/SR) >= 8 && FORMAT/BQ >= 30' \
  -Oz -o "$out_vcf" "$in_vcf"

bcftools index -t "$out_vcf"

echo "[SUCCESS] Filtered VCF created:"
echo "          $out_vcf"
