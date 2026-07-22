#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Standalone CBaSE v1.2 installation
#
# Official distribution:
#   https://weghornlab.org/Software/CBaSE_v1.2.zip
#
# The GitHub repository provides documentation; the executable
# and auxiliary reference files are distributed in the ZIP.
#
# Permanent installation:
#   tools/CBaSE/source/CBaSE_v1.2
#   tools/CBaSE/env_cbase
#
# The downloaded ZIP is temporary and is deleted automatically.
# ============================================================

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/CBaSE"
SOURCE_DIR="${TOOL_DIR}/source/CBaSE_v1.2"
ENV_DIR="${TOOL_DIR}/env_cbase"

DOWNLOAD_URLS=(
    "https://weghornlab.org/Software/CBaSE_v1.2.zip"
    "https://weghornlab.net/Software/CBaSE_v1.2.zip"
)

mkdir -p "${TOOL_DIR}/source"

module purge
module load anaconda3

BASE_PYTHON="$(command -v python3 || command -v python || true)"

if [[ -z "${BASE_PYTHON}" ]]; then
    echo "ERROR: Python was not found after loading anaconda3."
    exit 1
fi

echo "============================================================"
echo "Installing standalone CBaSE v1.2"
echo "============================================================"
echo "Tool directory: ${TOOL_DIR}"
echo "Source:         ${SOURCE_DIR}"
echo "Environment:    ${ENV_DIR}"
echo "Base Python:    ${BASE_PYTHON}"
"${BASE_PYTHON}" --version
echo "============================================================"

# ------------------------------------------------------------
# Download and install the official source only when absent
# ------------------------------------------------------------

if [[ -s "${SOURCE_DIR}/CBaSE_v1.2.py" ]] && \
   [[ -d "${SOURCE_DIR}/Auxiliary" ]]; then
    echo
    echo "Reusing existing CBaSE source and auxiliary resources:"
    echo "${SOURCE_DIR}"
else
    if [[ -e "${SOURCE_DIR}" ]]; then
        echo "ERROR: Source path exists but is incomplete:"
        echo "       ${SOURCE_DIR}"
        echo "Remove only this incomplete directory and rerun."
        exit 1
    fi

    TMP_ROOT="${TMPDIR:-/tmp}"
    TMP_WORK="$(mktemp -d "${TMP_ROOT}/cbase_install.XXXXXX")"

    cleanup() {
        rm -rf "${TMP_WORK}"
    }

    trap cleanup EXIT

    ZIP_FILE="${TMP_WORK}/CBaSE_v1.2.zip"
    EXTRACT_DIR="${TMP_WORK}/extract"
    mkdir -p "${EXTRACT_DIR}"

    DOWNLOADED=0

    for URL in "${DOWNLOAD_URLS[@]}"; do
        echo
        echo "Downloading: ${URL}"

        if curl \
            --fail \
            --location \
            --retry 3 \
            --retry-delay 5 \
            --connect-timeout 30 \
            --output "${ZIP_FILE}" \
            "${URL}"
        then
            DOWNLOADED=1
            break
        fi
    done

    if [[ "${DOWNLOADED}" -ne 1 ]] || [[ ! -s "${ZIP_FILE}" ]]; then
        echo "ERROR: CBaSE v1.2 could not be downloaded from the official URLs."
        exit 1
    fi

    if ! unzip -tq "${ZIP_FILE}" >/dev/null; then
        echo "ERROR: Downloaded CBaSE archive is not a valid ZIP file."
        exit 1
    fi

    unzip -q "${ZIP_FILE}" -d "${EXTRACT_DIR}"

    CBASE_SCRIPT="$(
        find "${EXTRACT_DIR}" \
            -type f \
            -name 'CBaSE_v1.2.py' \
            -print \
            -quit
    )"

    if [[ -z "${CBASE_SCRIPT}" ]]; then
        echo "ERROR: CBaSE_v1.2.py was not found in the downloaded archive."
        exit 1
    fi

    EXTRACTED_ROOT="$(dirname "${CBASE_SCRIPT}")"

    if [[ ! -d "${EXTRACTED_ROOT}/Auxiliary" ]]; then
        echo "ERROR: The official Auxiliary directory was not found beside:"
        echo "       ${CBASE_SCRIPT}"
        exit 1
    fi

    mkdir -p "$(dirname "${SOURCE_DIR}")"
    mv "${EXTRACTED_ROOT}" "${SOURCE_DIR}"

    echo
    echo "Installed official CBaSE source:"
    echo "${SOURCE_DIR}"
fi

# ------------------------------------------------------------
# Create a small dedicated Python environment
# ------------------------------------------------------------

if [[ -x "${ENV_DIR}/bin/python" ]]; then
    echo
    echo "Reusing existing CBaSE Python environment:"
    echo "${ENV_DIR}"
else
    if [[ -e "${ENV_DIR}" ]]; then
        echo "Removing incomplete CBaSE environment:"
        echo "${ENV_DIR}"
        rm -rf "${ENV_DIR}"
    fi

    echo
    echo "Creating CBaSE Python environment..."

    "${BASE_PYTHON}" -m venv "${ENV_DIR}"
fi

# Avoid retaining a second pip download cache.
export PIP_NO_CACHE_DIR=1

# Upgrade the environment packaging tools.
"${ENV_DIR}/bin/python" -m pip install \
    --upgrade \
    pip \
    setuptools \
    wheel

# Install or repair all CBaSE numerical and statistical dependencies.
"${ENV_DIR}/bin/python" -m pip install \
    numpy \
    scipy \
    mpmath

# ------------------------------------------------------------
# Validate source syntax, imports, and auxiliary resources
# ------------------------------------------------------------

echo
echo "Validating CBaSE installation..."

"${ENV_DIR}/bin/python" -m py_compile \
    "${SOURCE_DIR}/CBaSE_v1.2.py"

"${ENV_DIR}/bin/python" - <<'PYTHON'
import numpy
import scipy

print("NumPy:", numpy.__version__)
print("SciPy:", scipy.__version__)
print("Required Python imports: PASSED")
PYTHON

AUX_FILE_COUNT="$(
    find "${SOURCE_DIR}/Auxiliary" -type f | awk 'END {print NR+0}'
)"

if [[ "${AUX_FILE_COUNT}" -eq 0 ]]; then
    echo "ERROR: CBaSE Auxiliary directory contains no files."
    exit 1
fi

echo "Auxiliary files: ${AUX_FILE_COUNT}"

echo
echo "============================================================"
echo "CBaSE installation completed successfully"
echo "============================================================"
echo "Source:      ${SOURCE_DIR}"
echo "Environment: ${ENV_DIR}"
echo "Script:      ${SOURCE_DIR}/CBaSE_v1.2.py"
echo "Auxiliary:   ${SOURCE_DIR}/Auxiliary"
echo "============================================================"
