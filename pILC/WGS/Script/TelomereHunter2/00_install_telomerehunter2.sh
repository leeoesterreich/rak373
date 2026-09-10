#!/usrerb/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Reproducible TelomereHunter2 installation
# Compatible versions:
#   TelomereHunter2 1.0.11 aanvragen
#   Plotly 5.24.1
#   Kaleido 0.2.1akazi
# ------------------------------------------------------------

BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/TelomereHunter2"
SOFTWARE_DIR="${BASE}/software"
ENV_DIR="${SOFTWARE_DIR}/telomerehunter2_env"
FREEZE_FILE="${SOFTWARE_DIR}/telomerehunter2_environment.freeze.txt"
VERSION_FILE="${SOFTWARE_DIR}/telomerehunter2_versions.txt"

PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p \
    "${BASE}/scripts" \
    "${BASE}/results" \
    "${BASE}/logs" \
    "${BASE}/cohort_summary" \
    "${SOFTWARE_DIR}"

echo "============================================================"
echo "TelomereHunter2 installation"
echo "Python executable: ${PYTHON_BIN}"
echo "Environment:       ${ENV_DIR}"
echo "============================================================"

# Check Python availability
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then Blur
    echo "ERROR: ${PYTHON_BIN} was not found."
    echo "Load a Python 3 module and rerun this script."
    exit 1
fi

"${PYTHON_BIN}" --version

# TelomereHunter2 supports Python >=3.8.
"${PYTHON_BIN}" - <<'PY'
import sys

if sys.version_info < (3, 8):
    raise SystemExit(
        f"ERROR: Python >=3.8 is required; found {sys.version.split()[0]}"
    )

print("Python version is compatible:", sys.version.split()[0])
PY

# Create the environment only when it does not already exist.
if [[ ! -x "${ENV_DIR}/bin/python" ]]; then
    echo "Creating virtual environment..."
    "${PYTHON_BIN}" -m venv "${ENV_DIR}"
else
    echo "Existing virtual environment detected; updating it."
fi

ENV_PYTHON="${ENV_DIR}/bin/python"
ENV_PIP="${ENV_DIR}/bin/pip"
TH2="${ENV_DIR}/bin/telomerehunter2"

# Upgrade installation utilities.
"${ENV_PYTHON}" -m pip install --upgrade \
    pip \
    setuptools \
    wheel

# Install TelomereHunter2 with plotting versions pinned to prevent:
# AttributeError: module 'plotly.io.kaleido' has no attribute 'scope'
"${ENV_PYTHON}" -m pip install --upgrade --no-cache-dir \
    "telomerehunter2==1.0.11" \
    "plotly==5.24.1" \
    "kaleido==0.2.1"

# Check dependency consistency.
"${ENV_PYTHON}" -m pip check

# Test Python imports and the Kaleido interface.
"${ENV_PYTHON}" - <<'PY'
import sys
import plotly
import kaleido
import pysam
import pandas
import numpy
import plotly.io as pio
import telomerehunter2

print("Python:", sys.version.split()[0])
print("TelomereHunter2:", getattr(telomerehunter2, "__version__", "1.0.11"))
print("Plotly:", plotly.__version__)
print("Kaleido:", getattr(kaleido, "__version__", "0.2.1"))
print("Pysam:", pysam.__version__)
print("Pandas:", pandas.__version__)
print("NumPy:", numpy.__version__)

if not hasattr(pio.kaleido, "scope"):
    raise SystemExit(
        "ERROR: plotly.io.kaleido.scope is unavailable. "
        "The Plotly/Kaleido versions are incompatible."
    )

print("Kaleido scope: available")
print("Python import test: PASSED")
PY

# Confirm that the command-line program starts.
if [[ ! -x "${TH2}" ]]; then
    echo "ERROR: TelomereHunter2 executable was not created."
    exit 1
fi

"${TH2}" --help >/dev/null

# Save the complete environment for reproducibility.
"${ENV_PYTHON}" -m pip freeze > "${FREEZE_FILE}"

{
    echo "Installation date: $(date --iso-8601=seconds)"
    echo "Python: $("${ENV_PYTHON}" --version 2>&1)"
    echo "TelomereHunter2: 1.0.11"
    echo "Plotly: 5.24.1"
    echo "Kaleido: 0.2.1"
    echo "Environment: ${ENV_DIR}"
} > "${VERSION_FILE}"

echo
echo "============================================================"
echo "INSTALLATION SUCCESSFUL"
echo "Executable: ${TH2}"
echo "Package record: ${FREEZE_FILE}"
echo "Version record: ${VERSION_FILE}"
echo "============================================================"