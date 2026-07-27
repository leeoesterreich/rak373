#!/usr/bin/env bash
#SBATCH --job-name=OncodriveFML_pILC
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveFML/logs/OncodriveFML_pILC_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveFML/logs/OncodriveFML_pILC_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -Eeuo pipefail
trap 'echo "ERROR: Script failed at line ${LINENO}" >&2' ERR

# ============================================================
# Standalone OncodriveFML 2.5.0
# pILC WES Tier-1 common cohort
#
# Genome build: GRCh38
# Score file: CADD v1.7 GRCh38 PHRED
# Analysis: coding SNVs only
# ============================================================

PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveFML"

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveFML"
ENV_DIR="${TOOL_DIR}/env_oncodrivefml_2.5.0"
BGDATA_DIR="${TOOL_DIR}/bgdata"

INPUT="${PROJECT_DIR}/input/pILC_WES_oncodrivefml_input.tsv.gz"
ELEMENTS="${PROJECT_DIR}/reference/cds_hg38.tsv.gz"
CONFIG="${PROJECT_DIR}/config/oncodrivefml_hg38_cadd_v1.7.conf"

CADD_FILE="${TOOL_DIR}/resources/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
CADD_INDEX="${CADD_FILE}.tbi"

RESULTS_DIR="${PROJECT_DIR}/results"
LOG_DIR="${PROJECT_DIR}/logs"

RUN_ID="${SLURM_JOB_ID:-manual_$(date +%Y%m%d_%H%M%S)}"
RUN_OUTPUT="${RESULTS_DIR}/pILC_WES_coding_SNV_CADD_v1.7_${RUN_ID}"
RUN_LOG="${LOG_DIR}/pILC_WES_coding_SNV_CADD_v1.7_${RUN_ID}.log"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

# ------------------------------------------------------------
# Activate standalone Conda environment
# ------------------------------------------------------------

if ! type module >/dev/null 2>&1; then
    if [[ -r /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
    elif [[ -r /usr/share/lmod/lmod/init/bash ]]; then
        source /usr/share/lmod/lmod/init/bash
    else
        echo "ERROR: Could not initialize the cluster module system."
        exit 1
    fi
fi

module load anaconda3/2023.09-0-python_3.11.5

CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$ENV_DIR"

export BGDATA_LOCAL="$BGDATA_DIR"
export PYTHONNOUSERSITE=1
export PATH="${ENV_DIR}/bin:${PATH}"

hash -r

ONCODRIVEFML="${ENV_DIR}/bin/oncodrivefml"
TABIX="${ENV_DIR}/bin/tabix"

echo "============================================================"
echo "Standalone OncodriveFML run"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-manual}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-4}"
echo "Environment: ${CONDA_PREFIX:-not active}"
echo "============================================================"

# ------------------------------------------------------------
# Executable checks
# ------------------------------------------------------------

if [[ ! -x "$ONCODRIVEFML" ]]; then
    echo "ERROR: OncodriveFML executable is missing:"
    echo "  $ONCODRIVEFML"
    exit 1
fi

if [[ ! -x "$TABIX" ]]; then
    echo "ERROR: tabix executable is missing:"
    echo "  $TABIX"
    exit 1
fi

echo "OncodriveFML:"
"$ONCODRIVEFML" --version

echo "Tabix:"
"$TABIX" --version | head -n 1 || true

# ------------------------------------------------------------
# Required-file checks
# ------------------------------------------------------------

echo
echo "Checking required files..."

REQUIRED_FILES=(
    "$INPUT"
    "$ELEMENTS"
    "$CONFIG"
    "$CADD_FILE"
    "$CADD_INDEX"
)

for FILE in "${REQUIRED_FILES[@]}"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Required file is missing or empty:"
        echo "  $FILE"
        exit 1
    fi

    echo "FOUND: $FILE"
done

# ------------------------------------------------------------
# Configuration checks
# ------------------------------------------------------------

if ! grep -Eiq "^[[:space:]]*build[[:space:]]*=[[:space:]]*['\"]hg38['\"]" "$CONFIG"; then
    echo "ERROR: Configuration does not specify build = 'hg38'."
    exit 1
fi

if grep -Eiq "hg19" "$CONFIG"; then
    echo "ERROR: HG19 text was detected in the GRCh38 configuration."
    exit 1
fi

if ! grep -Fq "$CADD_FILE" "$CONFIG"; then
    echo "ERROR: Configuration does not point to the expected CADD file:"
    echo "  $CADD_FILE"
    exit 1
fi

if ! grep -Eiq "^[[:space:]]*format[[:space:]]*=[[:space:]]*['\"]tabix['\"]" "$CONFIG"; then
    echo "ERROR: Configuration does not specify format = 'tabix'."
    exit 1
fi

echo "PASS: GRCh38 configuration is consistent."

# ------------------------------------------------------------
# Input-format checks
# ------------------------------------------------------------

INPUT_HEADER="$(zgrep -m1 '^' "$INPUT")"
EXPECTED_INPUT_HEADER=$'CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE'

if [[ "$INPUT_HEADER" != "$EXPECTED_INPUT_HEADER" ]]; then
    echo "ERROR: Unexpected mutation-input header."
    echo "Observed: $INPUT_HEADER"
    echo "Expected: $EXPECTED_INPUT_HEADER"
    exit 1
fi

ELEMENT_HEADER="$(zgrep -m1 '^' "$ELEMENTS")"

case "$ELEMENT_HEADER" in
    $'CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL')
        ;;
    $'CHROMOSOME\tSTART\tEND\tELEMENT'*)
        ;;
    *)
        echo "ERROR: Unexpected genomic-elements header:"
        echo "  $ELEMENT_HEADER"
        exit 1
        ;;
esac

echo "PASS: Mutation and coding-elements headers are valid."

# ------------------------------------------------------------
# CADD Tabix check
# ------------------------------------------------------------

echo
echo "Checking CADD Tabix index..."

CONTIG_FILE="${LOG_DIR}/CADD_v1.7_GRCh38_contigs_${RUN_ID}.txt"
"$TABIX" -l "$CADD_FILE" > "$CONTIG_FILE"

if ! grep -qx "1" "$CONTIG_FILE"; then
    echo "ERROR: Chromosome 1 is absent from the CADD Tabix index."
    exit 1
fi

CADD_TEST="$("$TABIX" "$CADD_FILE" 1:1044014-1044014)"

if [[ -z "$CADD_TEST" ]]; then
    echo "ERROR: The CADD test query returned no records."
    exit 1
fi

echo "PASS: CADD Tabix query succeeded."
printf '%s\n' "$CADD_TEST" | sed -n '1p'

# ------------------------------------------------------------
# GRCh38 reference check
# ------------------------------------------------------------

echo
echo "Checking GRCh38 reference through bgreference..."

python - <<'PY'
from bgreference import refseq

base = refseq("hg38", "1", 10001, 1)

if not base or base.upper() not in {"A", "C", "G", "T", "N"}:
    raise RuntimeError(f"Unexpected GRCh38 reference result: {base!r}")

print(f"PASS: hg38 chr1:10001 reference base = {base}")
PY

# ------------------------------------------------------------
# Record run inputs
# ------------------------------------------------------------

mkdir -p "$RUN_OUTPUT"

cp "$CONFIG" "${RUN_OUTPUT}/used_oncodrivefml_hg38_cadd_v1.7.conf"

{
    echo "Run date: $(date)"
    echo "Host: $(hostname)"
    echo "SLURM job ID: ${SLURM_JOB_ID:-manual}"
    echo "OncodriveFML: $("$ONCODRIVEFML" --version 2>&1)"
    echo "Environment: $ENV_DIR"
    echo "Input: $INPUT"
    echo "Elements: $ELEMENTS"
    echo "Configuration: $CONFIG"
    echo "CADD: $CADD_FILE"
    echo "Output: $RUN_OUTPUT"
    echo "Cores: ${SLURM_CPUS_PER_TASK:-4}"
    echo "Seed: 123"
    echo "Sequencing: wes"
    echo "Signature correction: wx"
    echo "Indels: excluded"
} > "${RUN_OUTPUT}/run_manifest.txt"

# ------------------------------------------------------------
# Run OncodriveFML
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Running OncodriveFML..."
echo "Output: $RUN_OUTPUT"
echo "Log: $RUN_LOG"
echo "============================================================"

"$ONCODRIVEFML" \
    --input "$INPUT" \
    --elements "$ELEMENTS" \
    --sequencing wes \
    --signature-correction wx \
    --no-indels \
    --configuration "$CONFIG" \
    --output "$RUN_OUTPUT" \
    --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --seed 123 \
    --debug \
    --force \
    2>&1 | tee "$RUN_LOG"

# ------------------------------------------------------------
# Completion check
# ------------------------------------------------------------

echo
echo "Checking output files..."

mapfile -t OUTPUT_FILES < <(
    find "$RUN_OUTPUT" -maxdepth 2 -type f -size +0c | sort
)

if (( ${#OUTPUT_FILES[@]} == 0 )); then
    echo "ERROR: No nonempty output files were produced."
    exit 1
fi

printf '  %s\n' "${OUTPUT_FILES[@]}"

ln -sfn "$RUN_OUTPUT" "${RESULTS_DIR}/latest"

echo

# ------------------------------------------------------------
# Filter significant genes: Q_VALUE < 0.05
# ------------------------------------------------------------

FML_RESULTS="$(
    find "$RUN_OUTPUT" \
        -maxdepth 2 \
        -type f \
        -name '*-oncodrivefml.tsv' \
        -size +0c \
        | head -n 1
)"

if [[ -z "$FML_RESULTS" ]]; then
    echo "ERROR: OncodriveFML result TSV was not found."
    exit 1
fi

SIGNIFICANT_GENES="${RUN_OUTPUT}/significant_genes_Q_VALUE_lt_0.05.tsv"

awk -F'\t' '
    BEGIN { OFS = "\t" }

    NR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "Q_VALUE") {
                q_col = i
                break
            }
        }

        if (!q_col) {
            print "ERROR: Q_VALUE column not found." > "/dev/stderr"
            exit 1
        }

        print
        next
    }

    $q_col != "" &&
    $q_col != "NA" &&
    ($q_col + 0) < 0.05 {
        print
    }
' "$FML_RESULTS" > "$SIGNIFICANT_GENES"

SIGNIFICANT_GENE_COUNT="$(
    awk 'END {print NR > 0 ? NR - 1 : 0}' "$SIGNIFICANT_GENES"
)"

echo "Significant genes (Q_VALUE < 0.05): $SIGNIFICANT_GENE_COUNT"
echo "Filtered result: $SIGNIFICANT_GENES"

echo "============================================================"
echo "OncodriveFML completed successfully."
echo "Finished: $(date)"
echo "Results: $RUN_OUTPUT"
echo "Log: $RUN_LOG"
echo "============================================================"