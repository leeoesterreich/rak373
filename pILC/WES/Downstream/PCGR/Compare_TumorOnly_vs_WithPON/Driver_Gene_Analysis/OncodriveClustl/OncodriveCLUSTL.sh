#!/usr/bin/env bash
#SBATCH --job-name=ODCLUSTL_pILC
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveCLUSTL/logs/ODCLUSTL_pILC_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveCLUSTL/logs/ODCLUSTL_pILC_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -Eeuo pipefail
trap 'echo "ERROR: Script failed at line ${LINENO}" >&2' ERR

# ============================================================
# Standalone OncodriveCLUSTL 1.2.0
# pILC WES Tier-1 common cohort
#
# Input: coding SNVs only
# Genome: GRCh38/hg38
# Output: native OncodriveCLUSTL results only
# ============================================================

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveCLUSTL"

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveCLUSTL"

INPUT="${PROJECT_DIR}/input/pILC_WES_oncodriveclustl_SNV.tsv.gz"

# One canonical CDS-region file shared by FML and CLUSTL.
REGIONS="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveFML/resources/GRCh38/regions/cds.regions.gz"

RESULTS_DIR="${PROJECT_DIR}/results"
RUN_OUTPUT="${RESULTS_DIR}/pILC_WES_coding_SNV_CLUSTL"

ACTIVATE="${TOOL_DIR}/activate_oncodriveclustl.sh"
EXPECTED_VERSION="1.2.0"

mkdir -p "$RESULTS_DIR"

# ------------------------------------------------------------
# Activate OncodriveCLUSTL
# ------------------------------------------------------------

if [[ ! -r "$ACTIVATE" ]]; then
    echo "ERROR: Activation helper is missing:"
    echo "  $ACTIVATE"
    exit 1
fi

source "$ACTIVATE"

export MPLBACKEND=Agg

# Prevent BLAS/OpenMP oversubscription.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

ONCODRIVECLUSTL="$(command -v oncodriveclustl)"

VERSION="$(
    python - <<'PY'
from importlib.metadata import version
print(version("oncodriveclustl"))
PY
)"

if [[ "$VERSION" != "$EXPECTED_VERSION" ]]; then
    echo "ERROR: Expected OncodriveCLUSTL $EXPECTED_VERSION, found $VERSION."
    exit 1
fi

# ------------------------------------------------------------
# Check required files
# ------------------------------------------------------------

for FILE in "$INPUT" "$REGIONS"; do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Required file is missing or empty:"
        echo "  $FILE"
        exit 1
    fi

    gzip -t "$FILE"
done

# ------------------------------------------------------------
# Validate mutation input
# ------------------------------------------------------------

INPUT_HEADER="$(zgrep -m1 '^' "$INPUT")"
EXPECTED_INPUT_HEADER=$'CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE'

if [[ "$INPUT_HEADER" != "$EXPECTED_INPUT_HEADER" ]]; then
    echo "ERROR: Unexpected mutation-input header."
    echo "Observed: $INPUT_HEADER"
    echo "Expected: $EXPECTED_INPUT_HEADER"
    exit 1
fi

read -r SNV_COUNT NON_SNV_COUNT SAMPLE_COUNT < <(
    zcat "$INPUT" |
    awk -F'\t' '
        NR > 1 {
            samples[$5] = 1
            if ($3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/ && $3 != $4) {
                snv++
            } else {
                non_snv++
            }
        }
        END {
            for (sample in samples) {
                n_samples++
            }
            print snv+0, non_snv+0, n_samples+0
        }
    '
)
if [[ "$SNV_COUNT" -eq 0 ]]; then
    echo "ERROR: No valid SNVs were found."
    exit 1
fi

if [[ "$NON_SNV_COUNT" -ne 0 ]]; then
    echo "ERROR: Input contains $NON_SNV_COUNT non-SNV records."
    exit 1
fi

# ------------------------------------------------------------
# Validate region header
# ------------------------------------------------------------

REGION_HEADER="$(zgrep -m1 '^' "$REGIONS")"

for COLUMN in CHROMOSOME START END ELEMENT SYMBOL STRAND; do
    if ! printf '%s\n' "$REGION_HEADER" |
         tr '\t' '\n' |
         grep -Fxq "$COLUMN"; then

        echo "ERROR: Region file is missing column: $COLUMN"
        echo "Header: $REGION_HEADER"
        exit 1
    fi
done

# ------------------------------------------------------------
# Confirm GRCh38 reference access
# ------------------------------------------------------------

python - <<'PY'
from bgreference import refseq

base = refseq("hg38", "1", 10001, 1)

if not base:
    raise RuntimeError("Unable to access the hg38 reference.")

print(f"PASS: hg38 reference available; chr1:10001 = {base}")
PY

# ------------------------------------------------------------
# Prevent accidental result pileup
# ------------------------------------------------------------

if [[ -d "$RUN_OUTPUT" ]] &&
   find "$RUN_OUTPUT" -mindepth 1 -print -quit | grep -q .; then

    echo "ERROR: The result directory already contains files:"
    echo "  $RUN_OUTPUT"
    echo
    echo "Review or remove it before intentionally rerunning."
    exit 1
fi

mkdir -p "$RUN_OUTPUT"

# ------------------------------------------------------------
# Run OncodriveCLUSTL
# ------------------------------------------------------------

echo "============================================================"
echo "Running OncodriveCLUSTL"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Version: $VERSION"
echo "SNVs: $SNV_COUNT"
echo "Samples: $SAMPLE_COUNT"
echo "Input: $INPUT"
echo "Regions: $REGIONS"
echo "Output: $RUN_OUTPUT"
echo "============================================================"

"$ONCODRIVECLUSTL" \
    --input-file "$INPUT" \
    --regions-file "$REGIONS" \
    --output-directory "$RUN_OUTPUT" \
    --genome hg38 \
    --element-mutations 2 \
    --cluster-mutations 2 \
    --smooth-window 11 \
    --cluster-window 11 \
    --kmer 3 \
    --n-simulations 1000 \
    --simulation-mode region_restricted \
    --simulation-window 31 \
    --signature-calculation region_normalized \
    --cores "${SLURM_CPUS_PER_TASK:-4}" \
    --seed 123 \
    --concatenate \
    --log-level info

# ------------------------------------------------------------
# Check native outputs
# ------------------------------------------------------------

for FILE in \
    "${RUN_OUTPUT}/elements_results.txt" \
    "${RUN_OUTPUT}/clusters_results.tsv" \
    "${RUN_OUTPUT}/results.log"; do

    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Expected output is missing or empty:"
        echo "  $FILE"
        exit 1
    fi
done

echo "============================================================"
echo "OncodriveCLUSTL completed successfully."
echo "Finished: $(date)"
echo "Results: $RUN_OUTPUT"
echo "============================================================"