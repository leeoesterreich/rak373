#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Standalone SMRegions installation
#
# Official method:
#   clone repository
#   install cloned source with pip
#
# Uses a Python virtual environment, not conda create.
# Therefore, no Anaconda channel or ToS acceptance is needed.
# ============================================================

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/SMRegions"
SOURCE_DIR="${TOOL_DIR}/source/smregions"
ENV_DIR="${TOOL_DIR}/env_smregions"

REPOSITORY="https://github.com/bbglab/smregions.git"

mkdir -p "${TOOL_DIR}/source"

echo "============================================================"
echo "Installing standalone SMRegions"
echo "Source:      ${SOURCE_DIR}"
echo "Environment: ${ENV_DIR}"
echo "============================================================"

# ------------------------------------------------------------
# Load Python through the existing cluster module
# ------------------------------------------------------------

module purge
module load anaconda3

BASE_PYTHON="$(command -v python3 || command -v python)"

if [[ -z "${BASE_PYTHON}" ]]; then
    echo "ERROR: Python was not found after loading anaconda3."
    exit 1
fi

echo "Base Python: ${BASE_PYTHON}"
"${BASE_PYTHON}" --version

# ------------------------------------------------------------
# Clone only when the canonical source clone is absent
# ------------------------------------------------------------

if [[ -d "${SOURCE_DIR}/.git" ]]; then
    echo
    echo "Reusing existing SMRegions source clone:"
    echo "${SOURCE_DIR}"
else
    if [[ -e "${SOURCE_DIR}" ]]; then
        echo "ERROR: Source path exists but is not a Git repository:"
        echo "       ${SOURCE_DIR}"
        exit 1
    fi

    git clone "${REPOSITORY}" "${SOURCE_DIR}"
fi

SOURCE_COMMIT="$(
    git -C "${SOURCE_DIR}" rev-parse HEAD
)"

echo "Source commit: ${SOURCE_COMMIT}"

# ------------------------------------------------------------
# Reuse a valid completed installation
# ------------------------------------------------------------

if [[ -x "${ENV_DIR}/bin/smregions" ]]; then
    echo
    echo "SMRegions is already installed."

    "${ENV_DIR}/bin/smregions" --help >/dev/null

    echo "Executable: ${ENV_DIR}/bin/smregions"
    exit 0
fi

# Remove only an incomplete environment created by the failed
# installation attempt.
if [[ -d "${ENV_DIR}" ]]; then
    echo
    echo "Removing incomplete environment:"
    echo "${ENV_DIR}"

    rm -rf "${ENV_DIR}"
fi

# ------------------------------------------------------------
# Create a standard Python virtual environment
# ------------------------------------------------------------

echo
echo "Creating Python virtual environment..."

"${BASE_PYTHON}" -m venv "${ENV_DIR}"

# Avoid retaining a second pip download cache.
export PIP_NO_CACHE_DIR=1

"${ENV_DIR}/bin/python" -m pip install \
    --upgrade \
    pip \
    setuptools \
    wheel

# ------------------------------------------------------------
# Install from the canonical cloned repository
# ------------------------------------------------------------

echo
echo "Installing SMRegions from the cloned source..."

"${ENV_DIR}/bin/python" -m pip install \
    "${SOURCE_DIR}"

# ------------------------------------------------------------
# Validate installation
# ------------------------------------------------------------

echo
echo "Validating installation..."

"${ENV_DIR}/bin/smregions" --help >/dev/null
"${ENV_DIR}/bin/bgsignature" --help >/dev/null

"${ENV_DIR}/bin/python" - <<'PYTHON'
import smregions
import bgreference
import bgsignature
import numpy
import pandas
import scipy
import statsmodels

print("SMRegions version:", smregions.__version__)
print("Required Python imports: PASSED")
PYTHON

echo
echo "============================================================"
echo "SMRegions installation completed successfully"
echo "============================================================"
echo "Source:      ${SOURCE_DIR}"
echo "Commit:      ${SOURCE_COMMIT}"
echo "Environment: ${ENV_DIR}"
echo "Executable:  ${ENV_DIR}/bin/smregions"
echo "============================================================"