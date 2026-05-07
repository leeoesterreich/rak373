#!/bin/bash
#SBATCH --job-name=vcf2maf_tiers
#SBATCH --output=slurm_log/vcf2maf/vcf2maf_%A_%a.out
#SBATCH --error=slurm_log/vcf2maf/vcf2maf_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-2

set -euo pipefail

module load singularity/4.3.2

# ============================================================
# Purpose:
#   Convert tier-specific VCFs into MAF files
#
#   Tier_1_common          = common variants
#   Tier_2_With_PON_only   = With_PON-only variants
#   Tier_3_Tumor_only_only = Tumor-only-only variants
#
# Each Slurm array task processes one tier folder.
# Each tier folder contains multiple sample VCFs.
# ============================================================

# ============================================================
# Main directories
# ============================================================

COMPARE_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON"

OUT_ROOT="${COMPARE_DIR}/mafs"

mkdir -p "${OUT_ROOT}"

# ============================================================
# Tier folders
# ============================================================

TIERS=(
  "Tier_1_common"
  "Tier_2_With_PON_only"
  "Tier_3_Tumor_only_only"
)

TIER_PATTERNS=(
  "*.tier1.common.vcf.gz"
  "*.tier2.with_pon_only.vcf.gz"
  "*.tier3.tumor_only_only.vcf.gz"
)

TIER_SUFFIXES=(
  ".tier1.common.vcf.gz"
  ".tier2.with_pon_only.vcf.gz"
  ".tier3.tumor_only_only.vcf.gz"
)

MAF_SUFFIXES=(
  "tier1.common"
  "tier2.with_pon_only"
  "tier3.tumor_only_only"
)

# ============================================================
# vcf2maf / VEP / Singularity configuration
# ============================================================

VEP_DIR="/ix1/alee/LO_LAB/Personal/Brent_Schlegel/tools/VEP"
SIF="/ix1/alee/LO_LAB/Personal/Brent_Schlegel/envs/pcgr_2.2.3.sif"
REF_FASTA="/ix1/alee/LO_LAB/Personal/Brent_Schlegel/tools/VEP/homo_sapiens/113_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"

# ============================================================
# Select tier based on Slurm array task
# ============================================================

TIER_NAME="${TIERS[$SLURM_ARRAY_TASK_ID]}"
TIER_PATTERN="${TIER_PATTERNS[$SLURM_ARRAY_TASK_ID]}"
TIER_SUFFIX="${TIER_SUFFIXES[$SLURM_ARRAY_TASK_ID]}"
MAF_SUFFIX="${MAF_SUFFIXES[$SLURM_ARRAY_TASK_ID]}"

VCF_DIR="${COMPARE_DIR}/${TIER_NAME}"
OUT_DIR="${OUT_ROOT}/${TIER_NAME}"

mkdir -p "${OUT_DIR}"

SUMMARY_FILE="${OUT_DIR}/${TIER_NAME}_vcf2maf_summary.csv"

echo "============================================================"
echo "[$(date)] Starting VCF to MAF conversion"
echo "Tier:        ${TIER_NAME}"
echo "VCF_DIR:     ${VCF_DIR}"
echo "Pattern:     ${TIER_PATTERN}"
echo "OUT_DIR:     ${OUT_DIR}"
echo "SIF:         ${SIF}"
echo "VEP_DIR:     ${VEP_DIR}"
echo "REF_FASTA:   ${REF_FASTA}"
echo "============================================================"

# ============================================================
# Checks
# ============================================================

if [[ ! -d "${VCF_DIR}" ]]; then
  echo "[ERROR] VCF directory not found: ${VCF_DIR}"
  exit 1
fi

if [[ ! -f "${SIF}" ]]; then
  echo "[ERROR] Singularity image not found: ${SIF}"
  exit 1
fi

if [[ ! -d "${VEP_DIR}" ]]; then
  echo "[ERROR] VEP directory not found: ${VEP_DIR}"
  exit 1
fi

if [[ ! -f "${REF_FASTA}" ]]; then
  echo "[ERROR] Reference FASTA not found: ${REF_FASTA}"
  exit 1
fi

# ============================================================
# Find VCF files
# ============================================================

cd "${VCF_DIR}"

shopt -s nullglob
VCF_FILES=( ${TIER_PATTERN} )
shopt -u nullglob

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No VCF files found in ${VCF_DIR} with pattern ${TIER_PATTERN}"
  exit 1
fi

echo "[INFO] Number of VCF files found: ${#VCF_FILES[@]}"

# Summary header
echo "Sample,Tier,Input_VCF,Output_MAF,Status" > "${SUMMARY_FILE}"

# ============================================================
# Convert each VCF to MAF
# ============================================================

for VCF_FILE in "${VCF_FILES[@]}"; do

  SAMPLE="${VCF_FILE%${TIER_SUFFIX}}"
  tumor_base="${VCF_FILE%.vcf.gz}"

  OUT_MAF="${OUT_DIR}/${SAMPLE}.${MAF_SUFFIX}.maf"

  echo
  echo "------------------------------------------------------------"
  echo "[$(date)] Processing sample: ${SAMPLE}"
  echo "Input VCF:  ${VCF_DIR}/${VCF_FILE}"
  echo "Output MAF: ${OUT_MAF}"
  echo "------------------------------------------------------------"

  if [[ -f "${OUT_MAF}" ]]; then
    echo "[INFO] Removing existing MAF: ${OUT_MAF}"
    rm -f "${OUT_MAF}"
  fi

  set +e

  singularity exec \
    --writable-tmpfs \
    -B "${VCF_DIR}:/vcf_dir" \
    -B "${OUT_DIR}:/out_dir" \
    -B "${VEP_DIR}:/vep" \
    "${SIF}" \
    bash -c "
      set -euo pipefail

      gunzip -c /vcf_dir/${VCF_FILE} > /tmp/${tumor_base}.vcf

      vcf2maf.pl \
        --input-vcf /tmp/${tumor_base}.vcf \
        --output-maf /out_dir/${SAMPLE}.${MAF_SUFFIX}.maf \
        --tumor-id ${SAMPLE} \
        --vep-path /opt/mambaforge/envs/pcgr/bin/vep \
        --vep-data /vep/.vep \
        --ncbi-build GRCh38 \
        --ref-fasta ${REF_FASTA} \
        --inhibit-vep
    "

  STATUS_CODE=$?
  set -e

  if [[ ${STATUS_CODE} -eq 0 && -f "${OUT_MAF}" ]]; then
    echo "[DONE] Created MAF: ${OUT_MAF}"
    echo "${SAMPLE},${TIER_NAME},${VCF_DIR}/${VCF_FILE},${OUT_MAF},OK" >> "${SUMMARY_FILE}"
  else
    echo "[ERROR] Failed to create MAF for sample: ${SAMPLE}"
    echo "${SAMPLE},${TIER_NAME},${VCF_DIR}/${VCF_FILE},${OUT_MAF},ERROR" >> "${SUMMARY_FILE}"
  fi

done

echo
echo "============================================================"
echo "[$(date)] Completed tier: ${TIER_NAME}"
echo "Summary file: ${SUMMARY_FILE}"
echo "Output folder: ${OUT_DIR}"
echo "============================================================"