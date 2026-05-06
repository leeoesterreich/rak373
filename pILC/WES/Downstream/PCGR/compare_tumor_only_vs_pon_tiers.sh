#!/bin/bash
#SBATCH --job-name=pcgr_tier_compare
#SBATCH --output=slurm_log/pcgr_tier_compare_%A_%a.out
#SBATCH --error=slurm_log/pcgr_tier_compare_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-31

set -euo pipefail

# ============================================================
# Script: compare_tumor_only_vs_pon_tiers.sh
#
# Purpose:
#   Compare PCGR PASS VCFs from two Mutect2/PCGR workflows:
#
#   1. Tumor-only Mutect2 mode
#   2. Mutect2 with Panel of Normals (With_PON)
#
# Tier definitions:
#   Tier 1 = variants common in Tumor-only and With_PON
#            Tumor-only n With_PON
#
#   Tier 2 = variants present only in With_PON
#            With_PON - Tumor-only
#
#   Tier 3 = variants present only in Tumor-only
#            Tumor-only - With_PON
#
# Matching is based on exact CHROM, POS, REF, ALT using:
#   bcftools isec -c none
#
# ============================================================

# ============================================================
# Load modules
# ============================================================
module load bcftools/1.22

# ============================================================
# Directories
# ============================================================

BASE_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR"

TUMOR_ONLY_DIR="${BASE_DIR}/Tumor_only/pcgr_outputs"
WITH_PON_DIR="${BASE_DIR}/With_PON/pcgr_outputs"

OUT_DIR="${BASE_DIR}/Compare_TumorOnly_vs_WithPON"

SAMPLE_LIST="${OUT_DIR}/sample_list.txt"

TIER1_DIR="${OUT_DIR}/Tier_1_common"
TIER2_DIR="${OUT_DIR}/Tier_2_With_PON_only"
TIER3_DIR="${OUT_DIR}/Tier_3_Tumor_only_only"

SUMMARY_DIR="${OUT_DIR}/summaries"
LOG_DIR="${OUT_DIR}/logs"
TEMP_BASE="${OUT_DIR}/temp"

mkdir -p "${OUT_DIR}"
mkdir -p "${TIER1_DIR}" "${TIER2_DIR}" "${TIER3_DIR}"
mkdir -p "${SUMMARY_DIR}" "${LOG_DIR}" "${TEMP_BASE}"


# ============================================================
# Temporary directory for this task
# ============================================================

JOB_ID="${SLURM_JOB_ID:-manual}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-manual}"

TEMP_DIR="${TEMP_BASE}/tmp_${JOB_ID}_${TASK_ID}"
mkdir -p "${TEMP_DIR}"

cleanup() {
    rm -rf "${TEMP_DIR}"
}
trap cleanup EXIT


# ============================================================
# Helper functions
# ============================================================

log_info() {
    echo -e "[INFO] $(date '+%F %T') $1"
}

log_warn() {
    echo -e "[WARNING] $(date '+%F %T') $1"
}

log_error() {
    echo -e "[ERROR] $(date '+%F %T') $1"
}

ensure_index() {
    local vcf="$1"

    if [[ ! -f "${vcf}" ]]; then
        log_error "VCF not found: ${vcf}"
        return 1
    fi

    if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
        log_info "Index not found. Creating tabix index for: ${vcf}"
        bcftools index -t "${vcf}"
    else
        log_info "Index exists for: ${vcf}"
    fi
}

count_variants() {
    local vcf="$1"

    if [[ -f "${vcf}" ]]; then
        bcftools view -H "${vcf}" | wc -l | awk '{print $1}'
    else
        echo "NA"
    fi
}

run_isec() {
    local tier_name="$1"
    local n_rule="$2"
    local write_from="$3"
    local output_vcf="$4"
    local tumor_only_vcf="$5"
    local with_pon_vcf="$6"

    log_info "Creating ${tier_name}"
    log_info "Output: ${output_vcf}"

    # -c none requires exact CHROM, POS, REF, ALT match
    # -n=2 means present in both files
    # -n=1 means present in exactly one file
    # -w 1 writes records from file 1: Tumor-only
    # -w 2 writes records from file 2: With_PON

    bcftools isec \
        -c none \
        "${n_rule}" \
        -w "${write_from}" \
        -Oz \
        -o "${output_vcf}" \
        "${tumor_only_vcf}" \
        "${with_pon_vcf}"

    bcftools index -t -f "${output_vcf}"

    local count
    count=$(count_variants "${output_vcf}")

    log_info "${tier_name} variant count: ${count}"
}


# ============================================================
# Check sample list
# ============================================================

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    log_error "Sample list not found: ${SAMPLE_LIST}"
    log_error "Please create the sample list before submitting the job."
    log_error "Example:"
    log_error "find ${TUMOR_ONLY_DIR} -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort > ${SAMPLE_LIST}"
    exit 1
fi


# ============================================================
# Get sample for this Slurm array task
# ============================================================

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    log_error "SLURM_ARRAY_TASK_ID is not set. This script is intended for Slurm array execution."
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

if [[ -z "${SAMPLE}" ]]; then
    log_error "No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

log_info "============================================"
log_info "Processing sample: ${SAMPLE}"
log_info "SLURM_JOB_ID: ${SLURM_JOB_ID}"
log_info "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
log_info "============================================"


# ============================================================
# Define input VCFs
# ============================================================

TUMOR_ONLY_VCF="${TUMOR_ONLY_DIR}/${SAMPLE}/${SAMPLE}.pcgr.grch38.pass.vcf.gz"
WITH_PON_VCF="${WITH_PON_DIR}/${SAMPLE}/${SAMPLE}.pcgr.grch38.pass.vcf.gz"

log_info "Tumor-only VCF: ${TUMOR_ONLY_VCF}"
log_info "With_PON VCF:   ${WITH_PON_VCF}"


# ============================================================
# Define output VCFs
# ============================================================

TIER1_VCF="${TIER1_DIR}/${SAMPLE}.tier1.common.vcf.gz"
TIER2_VCF="${TIER2_DIR}/${SAMPLE}.tier2.with_pon_only.vcf.gz"
TIER3_VCF="${TIER3_DIR}/${SAMPLE}.tier3.tumor_only_only.vcf.gz"

SUMMARY_FILE="${SUMMARY_DIR}/${SAMPLE}.tier_summary.csv"


# ============================================================
# Check input files
# ============================================================

STATUS="PASS"

if [[ ! -f "${TUMOR_ONLY_VCF}" ]]; then
    log_error "Missing Tumor-only PASS VCF: ${TUMOR_ONLY_VCF}"
    STATUS="MISSING_TUMOR_ONLY"
fi

if [[ ! -f "${WITH_PON_VCF}" ]]; then
    log_error "Missing With_PON PASS VCF: ${WITH_PON_VCF}"

    if [[ "${STATUS}" == "PASS" ]]; then
        STATUS="MISSING_WITH_PON"
    else
        STATUS="MISSING_BOTH"
    fi
fi

if [[ "${STATUS}" != "PASS" ]]; then
    {
        echo "Sample,Status,Tumor_only_PASS,With_PON_PASS,Tier1_common,Tier2_With_PON_only,Tier3_Tumor_only_only,Tumor_only_VCF,With_PON_VCF"
        echo "${SAMPLE},${STATUS},NA,NA,NA,NA,NA,${TUMOR_ONLY_VCF},${WITH_PON_VCF}"
    } > "${SUMMARY_FILE}"

    log_error "Input file check failed for sample ${SAMPLE}. Summary written: ${SUMMARY_FILE}"
    exit 1
fi


# ============================================================
# Ensure VCFs are indexed
# ============================================================

ensure_index "${TUMOR_ONLY_VCF}"
ensure_index "${WITH_PON_VCF}"


# ============================================================
# Count original PASS variants
# ============================================================

TUMOR_ONLY_COUNT=$(count_variants "${TUMOR_ONLY_VCF}")
WITH_PON_COUNT=$(count_variants "${WITH_PON_VCF}")

log_info "Tumor-only PASS variant count: ${TUMOR_ONLY_COUNT}"
log_info "With_PON PASS variant count:   ${WITH_PON_COUNT}"


# ============================================================
# Create Tier 1, Tier 2, Tier 3 VCFs
# ============================================================

# ------------------------------------------------------------
# Tier 1:
# Common variants in both workflows
#
# Use -n=2 because the variant must be present in both VCFs.
# Use -w 2 to write the With_PON version of the record.
# This keeps the PON-informed PCGR annotation where the variant is shared.
# ------------------------------------------------------------

run_isec \
    "Tier 1 common: Tumor-only intersect With_PON" \
    "-n=2" \
    "2" \
    "${TIER1_VCF}" \
    "${TUMOR_ONLY_VCF}" \
    "${WITH_PON_VCF}"


# ------------------------------------------------------------
# Tier 2:
# With_PON-only variants
#
# Use -n=1 because the variant must be present in exactly one VCF.
# Use -w 2 to write variants from file 2, which is With_PON.
# This gives: With_PON - Tumor-only
# ------------------------------------------------------------

run_isec \
    "Tier 2 With_PON-only: With_PON minus Tumor-only" \
    "-n=1" \
    "2" \
    "${TIER2_VCF}" \
    "${TUMOR_ONLY_VCF}" \
    "${WITH_PON_VCF}"


# ------------------------------------------------------------
# Tier 3:
# Tumor-only-only variants
#
# Use -n=1 because the variant must be present in exactly one VCF.
# Use -w 1 to write variants from file 1, which is Tumor-only.
# This gives: Tumor-only - With_PON
# ------------------------------------------------------------

run_isec \
    "Tier 3 Tumor-only-only: Tumor-only minus With_PON" \
    "-n=1" \
    "1" \
    "${TIER3_VCF}" \
    "${TUMOR_ONLY_VCF}" \
    "${WITH_PON_VCF}"


# ============================================================
# Count tiered variants
# ============================================================

TIER1_COUNT=$(count_variants "${TIER1_VCF}")
TIER2_COUNT=$(count_variants "${TIER2_VCF}")
TIER3_COUNT=$(count_variants "${TIER3_VCF}")


# ============================================================
# Write per-sample summary
# ============================================================

{
    echo "Sample,Status,Tumor_only_PASS,With_PON_PASS,Tier1_common,Tier2_With_PON_only,Tier3_Tumor_only_only,Tumor_only_VCF,With_PON_VCF,Tier1_VCF,Tier2_VCF,Tier3_VCF"
    echo "${SAMPLE},${STATUS},${TUMOR_ONLY_COUNT},${WITH_PON_COUNT},${TIER1_COUNT},${TIER2_COUNT},${TIER3_COUNT},${TUMOR_ONLY_VCF},${WITH_PON_VCF},${TIER1_VCF},${TIER2_VCF},${TIER3_VCF}"
} > "${SUMMARY_FILE}"

log_info "Summary written: ${SUMMARY_FILE}"


# ============================================================
# Final check
# ============================================================

log_info "============================================"
log_info "Completed sample: ${SAMPLE}"
log_info "Tumor-only PASS variants: ${TUMOR_ONLY_COUNT}"
log_info "With_PON PASS variants:   ${WITH_PON_COUNT}"
log_info "Tier 1 common:            ${TIER1_COUNT}"
log_info "Tier 2 With_PON only:     ${TIER2_COUNT}"
log_info "Tier 3 Tumor-only only:   ${TIER3_COUNT}"
log_info "============================================"