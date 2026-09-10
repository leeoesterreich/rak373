#!/usr/bin/env bash
#SBATCH --job-name=TH2_pilot
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2/logs/TH2_pilot_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2/logs/TH2_pilot_%j.err

set -euo pipefail

# ============================================================
# Configuration
# ============================================================

BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2"

TH2="${BASE}/software/telomerehunter2_env/bin/telomerehunter2"

PAIR_ID="JC287-TP08-525_TP08-446"

TUMOR_BAM="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/4_BQSR/JC287-TP08-525/JC287-TP08-525.bqsr.bam"

NORMAL_BAM="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/4_BQSR/TP08-446/TP08-446.bqsr.bam"

OUTPUT_ROOT="${BASE}/results"
PAIR_OUTPUT="${OUTPUT_ROOT}/${PAIR_ID}"
SUMMARY="${PAIR_OUTPUT}/summary.tsv"

mkdir -p \
    "${OUTPUT_ROOT}" \
    "${BASE}/logs"

# ============================================================
# Increase open-file limit
# ============================================================

echo "Open-file limits before adjustment:"
echo "  Soft limit: $(ulimit -Sn)"
echo "  Hard limit: $(ulimit -Hn)"

HARD_LIMIT=$(ulimit -Hn)

if [[ "${HARD_LIMIT}" == "unlimited" ]]; then
    TARGET_LIMIT=65535
elif [[ "${HARD_LIMIT}" -ge 65535 ]]; then
    TARGET_LIMIT=65535
else
    TARGET_LIMIT="${HARD_LIMIT}"
fi

if [[ "${TARGET_LIMIT}" -lt 4096 ]]; then
    echo "ERROR: The maximum open-file limit is only ${TARGET_LIMIT}."
    echo "TelomereHunter2 may need more than 3,366 file descriptors."
    exit 1
fi

ulimit -Sn "${TARGET_LIMIT}"

echo "Open-file limits after adjustment:"
echo "  Soft limit: $(ulimit -Sn)"
echo "  Hard limit: $(ulimit -Hn)"

# ============================================================
# Check the TelomereHunter2 installation
# ============================================================

if [[ ! -x "${TH2}" ]]; then
    echo "ERROR: TelomereHunter2 executable not found:"
    echo "${TH2}"
    exit 1
fi

echo
echo "TelomereHunter2 executable:"
echo "${TH2}"

"${TH2}" --help >/dev/null

# ============================================================
# Check the input BAM files and indexes
# ============================================================

for bam in "${TUMOR_BAM}" "${NORMAL_BAM}"
do
    if [[ ! -s "${bam}" ]]; then
        echo "ERROR: BAM is missing or empty:"
        echo "${bam}"
        exit 1
    fi

    if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
        echo "ERROR: BAM index is missing:"
        echo "${bam}"
        exit 1
    fi
done

# ============================================================
# Handle previous output
# ============================================================

if [[ -s "${SUMMARY}" ]]; then
    echo
    echo "A completed result already exists:"
    echo "${SUMMARY}"
    echo "Nothing will be rerun."
    exit 0
fi

if [[ -d "${PAIR_OUTPUT}" ]]; then
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    FAILED_OUTPUT="${PAIR_OUTPUT}.failed_${TIMESTAMP}"

    echo
    echo "Incomplete output from a previous run was detected."
    echo "Moving:"
    echo "  ${PAIR_OUTPUT}"
    echo "to:"
    echo "  ${FAILED_OUTPUT}"

    mv "${PAIR_OUTPUT}" "${FAILED_OUTPUT}"
fi

# ============================================================
# Run paired tumor-normal TelomereHunter2 analysis
# ============================================================

echo
echo "============================================================"
echo "Starting TelomereHunter2 pilot"
echo "Slurm job: ${SLURM_JOB_ID:-unknown}"
echo "Pair:      ${PAIR_ID}"
echo "Genome:    hg38"
echo "Tumor:     ${TUMOR_BAM}"
echo "Normal:    ${NORMAL_BAM}"
echo "Output:    ${PAIR_OUTPUT}"
echo "CPUs:      ${SLURM_CPUS_PER_TASK}"
echo "============================================================"
echo

"${TH2}" \
    --pid "${PAIR_ID}" \
    --inputBamTumor "${TUMOR_BAM}" \
    --inputBamControl "${NORMAL_BAM}" \
    --outPath "${OUTPUT_ROOT}" \
    --bandingFile hg38 \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --parallel \
    --plotFileFormat html

# ============================================================
# Verify completion
# ============================================================

if [[ ! -s "${SUMMARY}" ]]; then
    echo
    echo "ERROR: TelomereHunter2 finished without creating:"
    echo "${SUMMARY}"
    exit 1
fi

echo
echo "============================================================"
echo "TELOMEREHUNTER2 PILOT COMPLETED SUCCESSFULLY"
echo "Summary: ${SUMMARY}"
echo "============================================================"