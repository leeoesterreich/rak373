#!/bin/bash
#SBATCH --job-name=Delly_merge
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --array=1-8
#SBATCH --output=slurm_log/Delly/Delly_merge_%A_%a.out
#SBATCH --error=slurm_log/Delly/Delly_merge_%A_%a.err

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
enrich_root="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"

# Ensure relative slurm_logs path works
cd "$enrich_root"

# Make sure SLURM log folder exists
mkdir -p slurm_logs/Delly

# ============================================================
# Resolve sample ID
# ============================================================
mapfile -t sample_dirs < <(find "$base_dir" -mindepth 1 -maxdepth 1 -type d | sort)

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

out_vcf="${sample_dir}/${sample}_delly_merged.vcf.gz"

if [ -f "$out_vcf" ]; then
  echo "[SKIP] Merged VCF already exists: $out_vcf"
  exit 0
fi

# ============================================================
# bgzip + tabix individual SVTYPE VCFs
# ============================================================
declare -a inputs=()

for sv in DEL DUP INV INS BND; do
  vcf="${sample_dir}/${sample}_${sv}.vcf"

  if [ -s "$vcf" ]; then
    echo "[INFO] Compressing + indexing: $vcf"
    bgzip -f "$vcf"
    tabix -f -p vcf "${vcf}.gz"
    inputs+=("${vcf}.gz")
  elif [ -s "${vcf}.gz" ]; then
    echo "[INFO] Found gzipped VCF: ${vcf}.gz (re-indexing)"
    tabix -f -p vcf "${vcf}.gz"
    inputs+=("${vcf}.gz")
  else
    echo "[WARN] Missing or empty: $vcf (or ${vcf}.gz)"
  fi
done

if [ "${#inputs[@]}" -eq 0 ]; then
  echo "[ERROR] No SVTYPE VCFs found for $sample"
  exit 1
fi

# ============================================================
# Merge SVTYPE VCFs into one Delly VCF
# ============================================================
echo "[INFO] Merging ${#inputs[@]} SVTYPE VCFs for $sample"
printf " - %s\n" "${inputs[@]}"

bcftools concat -a "${inputs[@]}" -Oz -o "$out_vcf"
bcftools index -t "$out_vcf"

echo "[SUCCESS] Created merged Delly VCF:"
echo "          $out_vcf"
