#!/bin/bash

# This script is used to filter Delly output vcf file.
# High-confidence filter:
            #  - FILTER must be PASS
            #  - at least 1 split-read (SR>=3)
            #  - PE + SR >= 8 total support
            #  - mapping quality for paired-ends and split-reads >= 30

#SBATCH --job-name=Delly_filter_merged
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --array=1-8
#SBATCH --output=slurm_log/Delly/Delly_filterMerged_%A_%a.out
#SBATCH --error=slurm_log/Delly/Delly_filterMerged_%A_%a.err

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
base_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Delly_010526"

# SLURM logs (relative to enrich_root)
mkdir -p slurm_logs/Delly

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

in_vcf="${sample_dir}/${sample}_delly_merged.vcf.gz"
out_vcf="${sample_dir}/${sample}_delly_filtered_merged.vcf.gz"

# ============================================================
# Checks
# ============================================================
if [ ! -s "$in_vcf" ]; then
  echo "[ERROR] Input merged VCF not found or empty: $in_vcf"
  exit 1
fi

if [ -s "$out_vcf" ]; then
  echo "[SKIP] Filtered merged VCF already exists: $out_vcf"
  exit 0
fi

# Make sure input is indexed (bcftools can still read without it, but good practice)
if [ ! -f "${in_vcf}.tbi" ]; then
  echo "[INFO] Indexing input merged VCF: $in_vcf"
  bcftools index -t "$in_vcf" || true
fi

# ============================================================
# Filtering (same criteria you used earlier)
# ============================================================
echo "[INFO] Filtering merged VCF with high-confidence rules:"
echo "       FILTER=PASS, SR>=3, PE+SR>=8, MAPQ>=30, SRMAPQ>=30"

bcftools view \
  -i 'FILTER="PASS" && INFO/SR>=3 && (INFO/PE + INFO/SR) >= 8 && INFO/MAPQ>=30 && INFO/SRMAPQ>=30' \
  -Oz -o "$out_vcf" "$in_vcf"

bcftools index -t "$out_vcf"

echo "[SUCCESS] Filtered merged VCF created:"
echo "          $out_vcf"
