#!/bin/bash

# This script is used to filter Lumpy output vcf file.
# High-confidence filter:
            #  - SR>=3
            #  - total support read >=8

#SBATCH --job-name=Lumpy_filter
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH -t 01:00:00
#SBATCH --array=1-48
#SBATCH --output=slurm_log/Lumpy/Lumpy_filter_%A_%a.out
#SBATCH --error=slurm_log/Lumpy/Lumpy_filter_%A_%a.err

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
base_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20260213/Batch2_021326/Lumpy_021326"

# SLURM logs
mkdir -p slurm_log/Lumpy

# ============================================================
# Resolve sample directory from array index
# ============================================================
mapfile -t sample_dirs < <(find "$base_dir" -mindepth 1 -maxdepth 1 -type d -name "TP*" | sort)

idx=$((SLURM_ARRAY_TASK_ID - 1))
if [ "$idx" -lt 0 ] || [ "$idx" -ge "${#sample_dirs[@]}" ]; then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID out of range (n=${#sample_dirs[@]})"
  exit 1
fi

sample_dir="${sample_dirs[$idx]}"
sample=$(basename "$sample_dir")

echo "============================================================"
echo "[INFO] SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "[INFO] Processing sample:   $sample"
echo "[INFO] Sample directory:    $sample_dir"
echo "============================================================"

# ============================================================
# Input / Output
# ============================================================
# Prefer .vcf.gz if it exists; otherwise use .vcf
if [ -s "${sample_dir}/${sample}.lumpy.vcf.gz" ]; then
  in_vcf="${sample_dir}/${sample}.lumpy.vcf.gz"
elif [ -s "${sample_dir}/${sample}.lumpy.vcf" ]; then
  in_vcf="${sample_dir}/${sample}.lumpy.vcf"
else
  echo "[ERROR] No input VCF found:"
  echo "        ${sample_dir}/${sample}.lumpy.vcf(.gz)"
  exit 1
fi

out_vcf="${sample_dir}/${sample}.lumpy_filtered.vcf"

# ============================================================
# Skip if already done
# ============================================================
if [ -s "$out_vcf" ]; then
  echo "[SKIP] Filtered VCF already exists: $out_vcf"
  exit 0
fi

# ============================================================
# Filtering
#   NOTE: Lumpy VCF often has FILTER="." not "PASS", so we do NOT require PASS here.
# ============================================================
echo "[INFO] Filtering Lumpy VCF with rules:"
echo "       INFO/SU>=8 && INFO/SR>=3"

bcftools view \
  -i 'INFO/SU>=8 && INFO/SR>=3' \
  -Ov -o "$out_vcf" "$in_vcf"

echo "[SUCCESS] Filtered Lumpy VCF created:"
echo "          $out_vcf"
