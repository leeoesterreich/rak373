#!/usr/bin/env bash
#SBATCH --job-name=CBaSE_pILC
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/Driver_Gene_Analysis/CBaSE/logs/CBaSE_pILC_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/Driver_Gene_Analysis/CBaSE/logs/CBaSE_pILC_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00

set -euo pipefail

# ============================================================
# Standalone CBaSE v1.2 analysis for the pILC WES cohort
#
# Primary analysis:
#   VCF mode = 1
#   genome   = hg38
#   context  = trinucleotide (CMODE=3)
#   model    = compare all six models (MODEL=0)
#
# CBaSE runs in node-local temporary space. Only the canonical
# q-value result is retained permanently; auxiliary outputs are
# summarized in the Slurm log and removed with temporary space.
# ============================================================

PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/Driver_Gene_Analysis/CBaSE"
TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/CBaSE"
SOURCE_DIR="${TOOL_DIR}/source/CBaSE_v1.2"
ENV_DIR="${TOOL_DIR}/env_cbase"

CBASE_SCRIPT="${SOURCE_DIR}/CBaSE_v1.2.py"
AUXILIARY_DIR="${SOURCE_DIR}/Auxiliary"
PYTHON="${ENV_DIR}/bin/python"

INPUT="${PROJECT_DIR}/input/pILC_WES_cbase_input.vcf"
RESULT_DIR="${PROJECT_DIR}/results"
OUTPUT="${RESULT_DIR}/pILC_WES_CBaSE_q_values.tsv"
SIGNIFICANT_OUTPUT="${RESULT_DIR}/pILC_WES_CBaSE_positive_selection_q_lt_0.05.tsv"

OUTNAME="pILC_WES"
VCF_MODE=1
BUILD="hg38"
CONTEXT_MODE=3
MODEL=0

EXPECTED_MUTATIONS=37662
EXPECTED_SAMPLES=31

mkdir -p \
    "${PROJECT_DIR}/input" \
    "${RESULT_DIR}" \
    "${PROJECT_DIR}/logs"

module purge
module load anaconda3

# Ensure CBaSE and all internally launched auxiliary scripts
# use the dedicated CBaSE Python environment.
export VIRTUAL_ENV="${ENV_DIR}"
export PATH="${ENV_DIR}/bin:${PATH}"
export PYTHONNOUSERSITE=1
hash -r

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "============================================================"
echo "Standalone CBaSE v1.2"
echo "============================================================"
echo "Job ID:   ${SLURM_JOB_ID:-NA}"
echo "Host:     $(hostname)"
echo "Start:    $(date)"
echo "Input:    ${INPUT}"
echo "Output:   ${OUTPUT}"
echo "Filtered: ${SIGNIFICANT_OUTPUT}"
echo "Build:    ${BUILD}"
echo "Context:  ${CONTEXT_MODE}-mer"
echo "Model:    ${MODEL} (compare all models)"
echo "============================================================"

# ------------------------------------------------------------
# Validate installation and input
# ------------------------------------------------------------

for FILE in "${CBASE_SCRIPT}" "${PYTHON}" "${INPUT}"; do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required file is missing or empty:"
        echo "       ${FILE}"
        exit 1
    fi
done

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: CBaSE Python executable is not executable:"
    echo "       ${PYTHON}"
    exit 1
fi

if [[ ! -d "${AUXILIARY_DIR}" ]]; then
    echo "ERROR: CBaSE Auxiliary directory is missing:"
    echo "       ${AUXILIARY_DIR}"
    exit 1
fi

EXPECTED_HEADER=$'#CHROM\tPOS\tID\tREF\tALT\tSAMPLE_ID'
OBSERVED_HEADER="$(sed -n '1p' "${INPUT}")"

if [[ "${OBSERVED_HEADER}" != "${EXPECTED_HEADER}" ]]; then
    echo "ERROR: Unexpected CBaSE input header."
    echo "Expected: ${EXPECTED_HEADER}"
    echo "Observed: ${OBSERVED_HEADER}"
    exit 1
fi

MUTATION_COUNT="$(awk 'NR > 1 {n++} END {print n+0}' "${INPUT}")"
SAMPLE_COUNT="$(
    awk -F'\t' 'NR > 1 {samples[$6]=1} END {print length(samples)}' \
        "${INPUT}"
)"

if [[ "${MUTATION_COUNT}" -ne "${EXPECTED_MUTATIONS}" ]]; then
    echo "ERROR: Unexpected CBaSE mutation count."
    echo "Expected: ${EXPECTED_MUTATIONS}"
    echo "Observed: ${MUTATION_COUNT}"
    exit 1
fi

if [[ "${SAMPLE_COUNT}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Unexpected CBaSE sample count."
    echo "Expected: ${EXPECTED_SAMPLES}"
    echo "Observed: ${SAMPLE_COUNT}"
    exit 1
fi

echo "Input mutations: ${MUTATION_COUNT}"
echo "Samples:         ${SAMPLE_COUNT}"

echo
echo "Python environment validation:"
echo "Dedicated Python: ${PYTHON}"
echo "Generic python:   $(command -v python)"
echo "Generic python3:  $(command -v python3)"

EXPECTED_PYTHON="$(readlink -f "${PYTHON}")"
OBSERVED_PYTHON="$(readlink -f "$(command -v python)")"

if [[ "${OBSERVED_PYTHON}" != "${EXPECTED_PYTHON}" ]]; then
    echo "ERROR: Generic python does not resolve to env_cbase."
    echo "Expected: ${EXPECTED_PYTHON}"
    echo "Observed: ${OBSERVED_PYTHON}"
    exit 1
fi

"${PYTHON}" - <<'PYTHON'
import sys
import mpmath
import numpy
import scipy

print("Python executable:", sys.executable)
print("NumPy:", numpy.__version__)
print("SciPy:", scipy.__version__)
print("mpmath:", mpmath.__version__)
print("CBaSE Python dependencies: PASSED")
PYTHON

# Resolve the installed mpmath package once. The two CBaSE
# auxiliary scripts will load this package directly by absolute
# path, independent of whichever Python executable CBaSE uses
# internally for its child processes.
MPMATH_PACKAGE_DIR="$(
    "${PYTHON}" - <<'PYTHON'
from pathlib import Path
import mpmath

print(Path(mpmath.__file__).resolve().parent)
PYTHON
)"

MPMATH_INIT="${MPMATH_PACKAGE_DIR}/__init__.py"

if [[ ! -s "${MPMATH_INIT}" ]]; then
    echo "ERROR: Installed mpmath package could not be resolved:"
    echo "       ${MPMATH_INIT}"
    exit 1
fi

echo "mpmath package:  ${MPMATH_PACKAGE_DIR}"


# ------------------------------------------------------------
# Create significant positive-selection result
# ------------------------------------------------------------

create_significant_output() {
    local source_file="$1"
    local destination_file="$2"
    local temporary_file="${destination_file}.tmp.${SLURM_JOB_ID:-$$}"

    "${PYTHON}" - \
        "${source_file}" \
        "${temporary_file}" <<'PYTHON'
import csv
import sys

source_path = sys.argv[1]
output_path = sys.argv[2]

with open(source_path, newline="") as source_handle:
    reader = csv.DictReader(source_handle, delimiter="\t")

    if not reader.fieldnames:
        raise RuntimeError("CBaSE result has no header.")

    required = {"gene", "q_phi_pos"}
    missing = required.difference(reader.fieldnames)

    if missing:
        raise RuntimeError(
            "CBaSE result is missing required columns: "
            + ", ".join(sorted(missing))
        )

    significant_rows = []

    for row in reader:
        try:
            q_value = float(row["q_phi_pos"])
        except (TypeError, ValueError):
            continue

        if q_value < 0.05:
            significant_rows.append((q_value, row))

significant_rows.sort(key=lambda item: item[0])

with open(output_path, "w", newline="") as output_handle:
    writer = csv.DictWriter(
        output_handle,
        fieldnames=reader.fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    for _, row in significant_rows:
        writer.writerow(row)

print(f"Significant positive-selection genes: {len(significant_rows)}")
PYTHON

    if [[ ! -s "${temporary_file}" ]]; then
        echo "ERROR: Significant positive-selection result was not generated."
        rm -f "${temporary_file}"
        exit 1
    fi

    if ! head -n 1 "${temporary_file}" | grep -q $'gene.*q_phi_pos'; then
        echo "ERROR: Significant positive-selection output header is invalid."
        head -n 3 "${temporary_file}"
        rm -f "${temporary_file}"
        exit 1
    fi

    mv "${temporary_file}" "${destination_file}"
}

# ------------------------------------------------------------
# Protect an existing completed result
# ------------------------------------------------------------

if [[ -s "${OUTPUT}" ]]; then
    if head -n 1 "${OUTPUT}" | grep -q $'gene.*q_phi_pos'; then
        echo
        echo "Completed CBaSE output already exists:"
        echo "${OUTPUT}"

        create_significant_output             "${OUTPUT}"             "${SIGNIFICANT_OUTPUT}"

        SIGNIFICANT_ROWS="$(
            awk 'NR > 1 {n++} END {print n+0}'                 "${SIGNIFICANT_OUTPUT}"
        )"

        echo
        echo "Significant positive-selection result prepared:"
        echo "${SIGNIFICANT_OUTPUT}"
        echo "Genes with q_phi_pos < 0.05: ${SIGNIFICANT_ROWS}"
        echo "No CBaSE rerun was performed."
        exit 0
    fi

    echo "ERROR: Existing output is not a valid completed CBaSE result:"
    echo "       ${OUTPUT}"
    echo "Review or remove only this file before rerunning."
    exit 1
fi

# ------------------------------------------------------------
# Isolated node-local run directory
# ------------------------------------------------------------

TMP_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TMP_WORK="${TMP_ROOT}/cbase_${SLURM_JOB_ID:-$$}"

mkdir -p \
    "${TMP_WORK}/Output" \
    "${TMP_WORK}/Auxiliary"

cleanup() {
    rm -rf "${TMP_WORK}"
}

trap cleanup EXIT

# ------------------------------------------------------------
# Prepare top-level CBaSE files
# ------------------------------------------------------------

# CBaSE expects to run beside its distributed resources.
# Link all top-level entries except Output and Auxiliary.
# The main script is copied because we also report its exact
# auxiliary-launch commands for troubleshooting.
while IFS= read -r -d '' ENTRY; do
    NAME="$(basename "${ENTRY}")"

    case "${NAME}" in
        Output|Auxiliary|__pycache__)
            continue
            ;;
        CBaSE_v1.2.py)
            cp "${ENTRY}" "${TMP_WORK}/${NAME}"
            ;;
        *)
            ln -s "${ENTRY}" "${TMP_WORK}/${NAME}"
            ;;
    esac
done < <(
    find "${SOURCE_DIR}" \
        -mindepth 1 \
        -maxdepth 1 \
        -print0
)

# ------------------------------------------------------------
# Prepare node-local Auxiliary directory
# ------------------------------------------------------------

# Copy only Python scripts because two scripts need a node-local
# import adjustment. Link all static reference resources directly
# from the canonical installation to avoid duplication.
while IFS= read -r -d '' ENTRY; do
    NAME="$(basename "${ENTRY}")"

    if [[ "${NAME}" == "__pycache__" ]]; then
        continue
    fi

    if [[ -f "${ENTRY}" && "${NAME}" == *.py ]]; then
        cp "${ENTRY}" "${TMP_WORK}/Auxiliary/${NAME}"
    else
        ln -s "${ENTRY}" "${TMP_WORK}/Auxiliary/${NAME}"
    fi
done < <(
    find "${AUXILIARY_DIR}" \
        -mindepth 1 \
        -maxdepth 1 \
        -print0
)

# ------------------------------------------------------------
# Patch the two mpmath-dependent auxiliary scripts
# ------------------------------------------------------------

# CBaSE v1.2 launches these scripts internally. On this cluster,
# that child interpreter does not inherit the dedicated virtual
# environment reliably. Replace only "import mpmath as mp" in the
# temporary copies with an absolute import of the validated mpmath
# package. The canonical source files remain unchanged.

for AUX_SCRIPT_NAME in \
    "CBaSE_v1.2_parameters.py" \
    "CBaSE_v1.2_qvalues.py"
do
    AUX_SCRIPT="${TMP_WORK}/Auxiliary/${AUX_SCRIPT_NAME}"

    if [[ ! -s "${AUX_SCRIPT}" ]]; then
        echo "ERROR: Required CBaSE auxiliary script is missing:"
        echo "       ${AUX_SCRIPT}"
        exit 1
    fi

    "${PYTHON}" - \
        "${AUX_SCRIPT}" \
        "${MPMATH_PACKAGE_DIR}" <<'PYTHON'
from pathlib import Path
import sys

script_path = Path(sys.argv[1])
mpmath_dir = Path(sys.argv[2]).resolve()
mpmath_init = mpmath_dir / "__init__.py"

text = script_path.read_text()
needle = "import mpmath as mp"

if text.count(needle) != 1:
    raise RuntimeError(
        f"Expected exactly one '{needle}' statement in {script_path}; "
        f"observed {text.count(needle)}."
    )

replacement = f"""# Load mpmath directly from the validated CBaSE environment.
# This node-local adjustment avoids dependence on the Python
# interpreter selected internally by CBaSE.
import importlib.util as _cbase_importlib_util
import sys as _cbase_sys

_cbase_mpmath_dir = {str(mpmath_dir)!r}
_cbase_mpmath_init = {str(mpmath_init)!r}

_cbase_spec = _cbase_importlib_util.spec_from_file_location(
    "mpmath",
    _cbase_mpmath_init,
    submodule_search_locations=[_cbase_mpmath_dir],
)

if _cbase_spec is None or _cbase_spec.loader is None:
    raise ImportError(
        "Unable to create an import specification for mpmath."
    )

mp = _cbase_importlib_util.module_from_spec(_cbase_spec)
_cbase_sys.modules["mpmath"] = mp
_cbase_spec.loader.exec_module(mp)
"""

script_path.write_text(
    text.replace(needle, replacement, 1)
)
PYTHON

done

# ------------------------------------------------------------
# Validate the temporary patched scripts
# ------------------------------------------------------------

if [[ ! -s "${TMP_WORK}/CBaSE_v1.2.py" ]]; then
    echo "ERROR: Temporary CBaSE main script was not prepared."
    exit 1
fi

if [[ ! -d "${TMP_WORK}/Auxiliary" ]]; then
    echo "ERROR: Temporary CBaSE Auxiliary directory was not prepared."
    exit 1
fi

for AUX_SCRIPT_NAME in \
    "CBaSE_v1.2_parameters.py" \
    "CBaSE_v1.2_qvalues.py"
do
    AUX_SCRIPT="${TMP_WORK}/Auxiliary/${AUX_SCRIPT_NAME}"

    if grep -q '^import mpmath as mp$' "${AUX_SCRIPT}"; then
        echo "ERROR: mpmath import was not patched in:"
        echo "       ${AUX_SCRIPT}"
        exit 1
    fi

    if ! grep -q '_cbase_spec.loader.exec_module(mp)' "${AUX_SCRIPT}"; then
        echo "ERROR: Absolute mpmath loader is missing from:"
        echo "       ${AUX_SCRIPT}"
        exit 1
    fi

    "${PYTHON}" -m py_compile "${AUX_SCRIPT}"
done

echo
echo "CBaSE internal auxiliary launch lines:"

grep -nE \
    'CBaSE_v1\.2_(parameters|qvalues)\.py' \
    "${TMP_WORK}/CBaSE_v1.2.py" || true

echo
echo "Node-local auxiliary patch validation: PASSED"


# ------------------------------------------------------------
# Run CBaSE
# ------------------------------------------------------------

echo
echo "Running CBaSE..."

cd "${TMP_WORK}"

"${PYTHON}" \
    "./CBaSE_v1.2.py" \
    "${INPUT}" \
    "${VCF_MODE}" \
    "${BUILD}" \
    "${CONTEXT_MODE}" \
    "${MODEL}" \
    "${OUTNAME}"

CBASE_OUTPUT="${TMP_WORK}/Output/q_values_${OUTNAME}.txt"
MODEL_OUTPUT="${TMP_WORK}/Output/used_params_and_model_${OUTNAME}.txt"

if [[ ! -s "${CBASE_OUTPUT}" ]]; then
    echo "ERROR: CBaSE q-value output was not generated:"
    echo "       ${CBASE_OUTPUT}"
    exit 1
fi

# CBaSE writes the selected-model summary before the tabular
# gene-result header. Locate the real header instead of assuming
# that it is the first line.
HEADER_LINE="$(
    awk -F'\t' '
        $1 == "gene" {
            for (i = 1; i <= NF; i++) {
                if ($i == "q_phi_pos") {
                    print NR
                    exit
                }
            }
        }
    ' "${CBASE_OUTPUT}"
)"

if [[ -z "${HEADER_LINE}" ]]; then
    echo "ERROR: CBaSE output does not contain a gene header with q_phi_pos."
    echo "Output preview:"
    sed -n '1,5p' "${CBASE_OUTPUT}"
    exit 1
fi

echo
echo "CBaSE result header found at line: ${HEADER_LINE}"

# Preserve any model/parameter preamble in the Slurm log, but
# write only the clean gene-level table to the canonical result.
if [[ "${HEADER_LINE}" -gt 1 ]]; then
    echo
    echo "Chosen model and parameter summary:"
    sed -n "1,$((HEADER_LINE - 1))p" "${CBASE_OUTPUT}"
elif [[ -s "${MODEL_OUTPUT}" ]]; then
    echo
    echo "Chosen model and parameter summary:"
    cat "${MODEL_OUTPUT}"
fi

TEMP_RESULT="${OUTPUT}.tmp.${SLURM_JOB_ID:-$$}"

tail -n +"${HEADER_LINE}" "${CBASE_OUTPUT}" > "${TEMP_RESULT}"

if [[ ! -s "${TEMP_RESULT}" ]]; then
    echo "ERROR: Clean CBaSE result could not be generated."
    rm -f "${TEMP_RESULT}"
    exit 1
fi

if ! awk -F'\t' '
    NR == 1 {
        if ($1 != "gene") {
            exit 1
        }

        found = 0
        for (i = 1; i <= NF; i++) {
            if ($i == "q_phi_pos") {
                found = 1
            }
        }

        exit(found ? 0 : 1)
    }
' "${TEMP_RESULT}"; then
    echo "ERROR: Clean CBaSE result header is invalid."
    head -n 3 "${TEMP_RESULT}"
    rm -f "${TEMP_RESULT}"
    exit 1
fi

mv "${TEMP_RESULT}" "${OUTPUT}"

RESULT_ROWS="$(awk 'NR > 1 {n++} END {print n+0}' "${OUTPUT}")"

create_significant_output     "${OUTPUT}"     "${SIGNIFICANT_OUTPUT}"

SIGNIFICANT_ROWS="$(
    awk 'NR > 1 {n++} END {print n+0}'         "${SIGNIFICANT_OUTPUT}"
)"

# ------------------------------------------------------------
# Summarize significant positive-selection genes in the log
# ------------------------------------------------------------

echo
echo "Genes with combined positive-selection q-value < 0.05:"

"${PYTHON}" - "${OUTPUT}" <<'PYTHON'
import csv
import sys

path = sys.argv[1]

with open(path, newline="") as handle:
    sample = handle.read(4096)
    handle.seek(0)

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t ,")
    except csv.Error:
        dialect = csv.excel_tab

    reader = csv.DictReader(handle, dialect=dialect)

    if not reader.fieldnames:
        raise RuntimeError("CBaSE result has no header.")

    required = {"gene", "q_phi_pos"}
    missing = required.difference(reader.fieldnames)
    if missing:
        raise RuntimeError(
            "CBaSE result is missing required columns: "
            + ", ".join(sorted(missing))
        )

    rows = []
    for row in reader:
        try:
            q_value = float(row["q_phi_pos"])
        except (TypeError, ValueError):
            continue

        if q_value < 0.05:
            rows.append((q_value, row))

    rows.sort(key=lambda item: item[0])

    if not rows:
        print("None")
    else:
        preferred = [
            "gene",
            "q_phi_pos",
            "p_phi_pos",
            "phi_pos_or_p(m=0|s)*p(k=0|s)",
            "m_obs",
            "k_obs",
            "s_obs",
            "dm/ds",
            "dk/ds",
            "d(m+k)/ds",
        ]
        columns = [column for column in preferred if column in reader.fieldnames]
        print("\t".join(columns))
        for _, row in rows:
            print("\t".join(row.get(column, "") for column in columns))

    print(f"Significant genes: {len(rows)}")
PYTHON

echo
echo "============================================================"
echo "CBaSE completed successfully"
echo "============================================================"
echo "Input mutations: ${MUTATION_COUNT}"
echo "Samples:         ${SAMPLE_COUNT}"
echo "Result rows:     ${RESULT_ROWS}"
echo "Significant:     ${SIGNIFICANT_ROWS} genes with q_phi_pos < 0.05"
echo "Output:          ${OUTPUT}"
echo "Filtered output: ${SIGNIFICANT_OUTPUT}"
echo "End:             $(date)"
echo "============================================================"

echo
echo "Result preview:"
sed -n '1,10p' "${OUTPUT}"

echo
echo "Significant positive-selection result:"
sed -n '1,10p' "${SIGNIFICANT_OUTPUT}"