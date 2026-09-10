#!/usr/bin/env bash
#SBATCH --job-name=TH2_pairs
#SBATCH --array=1-10%5
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2/logs/TH2_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2/logs/TH2_%A_%a.err

set -euo pipefail

###############################################################################
# Paths
###############################################################################

BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2"

CSV="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/4_BQSR/Paired_samples.csv"

TH2="${BASE}/software/telomerehunter2_env/bin/telomerehunter2"

OUTPUT_ROOT="${BASE}/results"

mkdir -p "${OUTPUT_ROOT}"

###############################################################################
# Raise open-file limit
###############################################################################

echo "Open-file limits before adjustment:"
echo "Soft: $(ulimit -Sn)"
echo "Hard: $(ulimit -Hn)"

HARD_LIMIT="$(ulimit -Hn)"

if [[ "${HARD_LIMIT}" == "unlimited" ]]; then
    TARGET_LIMIT=65535
elif (( HARD_LIMIT >= 65535 )); then
    TARGET_LIMIT=65535
else
    TARGET_LIMIT="${HARD_LIMIT}"
fi

ulimit -Sn "${TARGET_LIMIT}"

echo "Open-file limits after adjustment:"
echo "Soft: $(ulimit -Sn)"
echo "Hard: $(ulimit -Hn)"

if (( TARGET_LIMIT < 8192 )); then
    echo "ERROR: Open-file limit is too low: ${TARGET_LIMIT}"
    exit 1
fi

###############################################################################
# Check software
###############################################################################

if [[ ! -x "${TH2}" ]]; then
    echo "ERROR: TelomereHunter2 executable was not found:"
    echo "${TH2}"
    exit 1
fi

if [[ ! -s "${CSV}" ]]; then
    echo "ERROR: Paired-sample CSV was not found or is empty:"
    echo "${CSV}"
    exit 1
fi

###############################################################################
# Extract the CSV row corresponding to this array task
#
# Array task 1 = CSV line 2
# Array task 2 = CSV line 3
# ...
# Array task 10 = CSV line 11
###############################################################################

CSV_LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))

ROW="$(
    awk -v line_number="${CSV_LINE_NUMBER}" \
        'NR == line_number {gsub(/\r$/, ""); print}' \
        "${CSV}"
)"

if [[ -z "${ROW}" ]]; then
    echo "ERROR: No CSV entry found for array task ${SLURM_ARRAY_TASK_ID}"
    echo "Expected CSV line: ${CSV_LINE_NUMBER}"
    exit 1
fi

IFS=',' read -r \
    PAIR_ID \
    TUMOR_ID \
    NORMAL_ID \
    TUMOR_BAM \
    NORMAL_BAM <<< "${ROW}"

###############################################################################
# Validate parsed fields
###############################################################################

for value in \
    "${PAIR_ID}" \
    "${TUMOR_ID}" \
    "${NORMAL_ID}" \
    "${TUMOR_BAM}" \
    "${NORMAL_BAM}"
do
    if [[ -z "${value}" ]]; then
        echo "ERROR: One or more CSV fields are empty."
        echo "CSV row: ${ROW}"
        exit 1
    fi
done

echo "============================================================"
echo "Slurm job ID:       ${SLURM_JOB_ID}"
echo "Array task ID:      ${SLURM_ARRAY_TASK_ID}"
echo "Pair ID:            ${PAIR_ID}"
echo "Tumor ID:           ${TUMOR_ID}"
echo "Normal ID:          ${NORMAL_ID}"
echo "Tumor BAM:          ${TUMOR_BAM}"
echo "Normal BAM:         ${NORMAL_BAM}"
echo "CPUs allocated:     ${SLURM_CPUS_PER_TASK}"
echo "Reference:          hg38"
echo "============================================================"

###############################################################################
# Confirm filenames agree with sample IDs
###############################################################################

EXPECTED_TUMOR_NAME="${TUMOR_ID}.bqsr.bam"
EXPECTED_NORMAL_NAME="${NORMAL_ID}.bqsr.bam"

if [[ "$(basename "${TUMOR_BAM}")" != "${EXPECTED_TUMOR_NAME}" ]]; then
    echo "ERROR: Tumor BAM filename does not match tumor sample ID."
    echo "Expected: ${EXPECTED_TUMOR_NAME}"
    echo "Observed: $(basename "${TUMOR_BAM}")"
    exit 1
fi

if [[ "$(basename "${NORMAL_BAM}")" != "${EXPECTED_NORMAL_NAME}" ]]; then
    echo "ERROR: Normal BAM filename does not match normal sample ID."
    echo "Expected: ${EXPECTED_NORMAL_NAME}"
    echo "Observed: $(basename "${NORMAL_BAM}")"
    exit 1
fi

###############################################################################
# Confirm BAM files and indexes
###############################################################################

check_bam()
{
    local sample_type="$1"
    local bam="$2"

    if [[ ! -f "${bam}" ]]; then
        echo "ERROR: ${sample_type} BAM does not exist:"
        echo "${bam}"
        exit 1
    fi

    if [[ ! -s "${bam}" ]]; then
        echo "ERROR: ${sample_type} BAM is empty:"
        echo "${bam}"
        exit 1
    fi

    if [[ ! -r "${bam}" ]]; then
        echo "ERROR: ${sample_type} BAM is not readable:"
        echo "${bam}"
        exit 1
    fi

    if [[ -s "${bam}.bai" ]]; then
        index="${bam}.bai"
    elif [[ -s "${bam%.bam}.bai" ]]; then
        index="${bam%.bam}.bai"
    elif [[ -s "${bam}.csi" ]]; then
        index="${bam}.csi"
    elif [[ -s "${bam%.bam}.csi" ]]; then
        index="${bam%.bam}.csi"
    else
        echo "ERROR: Index not found for ${sample_type} BAM:"
        echo "${bam}"
        exit 1
    fi

    echo "${sample_type} BAM passed validation."
    echo "${sample_type} index: ${index}"
}

check_bam "Tumor" "${TUMOR_BAM}"
check_bam "Normal" "${NORMAL_BAM}"

###############################################################################
# Handle existing output
###############################################################################

PAIR_OUTPUT="${OUTPUT_ROOT}/${PAIR_ID}"

# This is the correct TelomereHunter2 combined-summary filename.
SUMMARY="${PAIR_OUTPUT}/${PAIR_ID}_summary.tsv"

if [[ -s "${SUMMARY}" ]]; then
    echo "Completed result already exists."
    echo "Skipping this pair to avoid overwriting it."
    echo "Summary: ${SUMMARY}"
    exit 0
fi

# Preserve an incomplete output directory from an earlier failed run.
if [[ -d "${PAIR_OUTPUT}" ]]; then
    TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
    BACKUP_OUTPUT="${PAIR_OUTPUT}.incomplete_${SLURM_JOB_ID}_${TIMESTAMP}"

    echo "Incomplete output directory detected."
    echo "Moving it to:"
    echo "${BACKUP_OUTPUT}"

    mv "${PAIR_OUTPUT}" "${BACKUP_OUTPUT}"
fi

###############################################################################
# Run TelomereHunter2
###############################################################################

START_TIME="$(date '+%Y-%m-%d %H:%M:%S')"

echo "TelomereHunter2 started: ${START_TIME}"

"${TH2}" \
    --pid "${PAIR_ID}" \
    --inputBamTumor "${TUMOR_BAM}" \
    --inputBamControl "${NORMAL_BAM}" \
    --outPath "${OUTPUT_ROOT}" \
    --bandingFile hg38 \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --parallel \
    --plotFileFormat html

###############################################################################
# Confirm successful completion
###############################################################################

if [[ ! -s "${SUMMARY}" ]]; then
    echo "ERROR: TelomereHunter2 finished without generating the summary:"
    echo "${SUMMARY}"
    exit 1
fi

END_TIME="$(date '+%Y-%m-%d %H:%M:%S')"

echo "============================================================"
echo "TelomereHunter2 completed successfully."
echo "Pair:       ${PAIR_ID}"
echo "Started:    ${START_TIME}"
echo "Finished:   ${END_TIME}"
echo "Summary:    ${SUMMARY}"
echo "HTML report: ${PAIR_OUTPUT}/all_reports.html"
echo "============================================================"

echo "Summary contents:"
column -t -s $'\t' "${SUMMARY}" || cat "${SUMMARY}"