#!/usr/bin/env bash

set -Eeuo pipefail
trap 'echo "ERROR: Installation failed at line ${LINENO}." >&2' ERR

# ============================================================
# Standalone OncodriveCLUSTL 1.2.0 installation
# Pitt CRC / Lmod + Anaconda
#
# This creates an environment independent of IntOGen and
# independent of the standalone OncodriveFML environment.
#
# OncodriveCLUSTL 1.2.0 requires Python >=3.10 and <3.13.
# ============================================================

# ------------------------------------------------------------
# User-editable settings
# ------------------------------------------------------------

ANACONDA_MODULE="anaconda3/2023.09-0-python_3.11.5"

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveCLUSTL"
ENV_DIR="${TOOL_DIR}/env_oncodriveclustl_1.2.0"
BGDATA_DIR="${TOOL_DIR}/bgdata"

ONCODRIVECLUSTL_VERSION="1.2.0"
PYTHON_VERSION="3.10"

LOG_DIR="${TOOL_DIR}/installation_logs"
RECORD_DIR="${TOOL_DIR}/installation_records"

# ------------------------------------------------------------
# Prepare directories and logging
# ------------------------------------------------------------

mkdir -p \
    "$TOOL_DIR" \
    "$BGDATA_DIR" \
    "$LOG_DIR" \
    "$RECORD_DIR"

RUN_DATE="$(date +%Y%m%d_%H%M%S)"
INSTALL_LOG="${LOG_DIR}/install_oncodriveclustl_${RUN_DATE}.log"

exec > >(tee -a "$INSTALL_LOG") 2>&1

 echo "============================================================"
 echo "Standalone OncodriveCLUSTL installation"
 echo "Started: $(date)"
 echo "Host: $(hostname)"
 echo "Tool directory: $TOOL_DIR"
 echo "Environment: $ENV_DIR"
 echo "BGData directory: $BGDATA_DIR"
 echo "Installation log: $INSTALL_LOG"
 echo "============================================================"

# Prevent accidental imports from ~/.local and prevent GUI use.
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg

# ------------------------------------------------------------
# Initialize the module system
# ------------------------------------------------------------

if ! type module >/dev/null 2>&1; then
    if [[ -r /etc/profile.d/modules.sh ]]; then
        # shellcheck source=/dev/null
        source /etc/profile.d/modules.sh
    elif [[ -r /usr/share/lmod/lmod/init/bash ]]; then
        # shellcheck source=/dev/null
        source /usr/share/lmod/lmod/init/bash
    else
        echo "ERROR: The cluster module system could not be initialized."
        exit 1
    fi
fi

# ------------------------------------------------------------
# Load Anaconda and initialize Conda
# ------------------------------------------------------------

 echo
 echo "Loading module: $ANACONDA_MODULE"
module load "$ANACONDA_MODULE"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found after loading $ANACONDA_MODULE."
    exit 1
fi

CONDA_BASE="$(conda info --base)"

if [[ ! -r "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: Conda initialization script was not found:"
    echo "  ${CONDA_BASE}/etc/profile.d/conda.sh"
    exit 1
fi

# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"

 echo "Conda executable: $(command -v conda)"
 echo "Conda version: $(conda --version)"
 echo "Conda base: $CONDA_BASE"

# ------------------------------------------------------------
# Create the independent environment when absent
# ------------------------------------------------------------

if [[ ! -x "${ENV_DIR}/bin/python" ]]; then
    echo
    echo "Creating Conda environment..."

    conda create \
        --override-channels \
        --prefix "$ENV_DIR" \
        --channel conda-forge \
        "python=${PYTHON_VERSION}" \
        pip \
        setuptools \
        wheel \
        -y
else
    echo
    echo "Existing environment detected; creation skipped:"
    echo "  $ENV_DIR"
fi

# ------------------------------------------------------------
# Activate and validate the environment
# ------------------------------------------------------------

conda activate "$ENV_DIR"
hash -r

if [[ "${CONDA_PREFIX:-}" != "$ENV_DIR" ]]; then
    echo "ERROR: Unexpected active Conda environment:"
    echo "  ${CONDA_PREFIX:-not set}"
    exit 1
fi

 echo
 echo "Active environment: $CONDA_PREFIX"
 echo "Python executable: $(command -v python)"
python --version

if [[ "$(command -v python)" != "${ENV_DIR}/bin/python" ]]; then
    echo "ERROR: Python is not resolving from the intended environment."
    exit 1
fi

python - <<'PY'
import sys
if not ((3, 10) <= sys.version_info[:2] < (3, 13)):
    raise RuntimeError(
        f"OncodriveCLUSTL 1.2.0 requires Python >=3.10,<3.13; found {sys.version}"
    )
print("PASS: Python version is compatible with OncodriveCLUSTL 1.2.0.")
PY

# ------------------------------------------------------------
# Install or verify OncodriveCLUSTL
# ------------------------------------------------------------

INSTALLED_VERSION="$({
    python - <<'PY' 2>/dev/null || true
try:
    from importlib.metadata import version
    print(version("oncodriveclustl"))
except Exception:
    pass
PY
} | tail -n 1)"

if [[ "$INSTALLED_VERSION" == "$ONCODRIVECLUSTL_VERSION" ]]; then
    echo
    echo "OncodriveCLUSTL ${ONCODRIVECLUSTL_VERSION} is already installed."
else
    echo
    echo "Installing OncodriveCLUSTL ${ONCODRIVECLUSTL_VERSION}..."

    python -m pip install \
        --no-cache-dir \
        --upgrade \
        "oncodriveclustl==${ONCODRIVECLUSTL_VERSION}"
fi

hash -r

# ------------------------------------------------------------
# Verify executable and package version
# ------------------------------------------------------------

if ! command -v oncodriveclustl >/dev/null 2>&1; then
    echo "ERROR: oncodriveclustl executable was not found."
    exit 1
fi

if [[ "$(command -v oncodriveclustl)" != "${ENV_DIR}/bin/oncodriveclustl" ]]; then
    echo "ERROR: OncodriveCLUSTL is not resolving from the intended environment."
    echo "Found: $(command -v oncodriveclustl)"
    exit 1
fi

DETECTED_VERSION="$(python - <<'PY'
from importlib.metadata import version
print(version("oncodriveclustl"))
PY
)"

 echo
 echo "OncodriveCLUSTL executable: $(command -v oncodriveclustl)"
 echo "OncodriveCLUSTL package version: $DETECTED_VERSION"

if [[ "$DETECTED_VERSION" != "$ONCODRIVECLUSTL_VERSION" ]]; then
    echo "ERROR: Expected OncodriveCLUSTL ${ONCODRIVECLUSTL_VERSION}, but found:"
    echo "  $DETECTED_VERSION"
    exit 1
fi

oncodriveclustl --help >/dev/null

# ------------------------------------------------------------
# Configure environment variables on activation
# ------------------------------------------------------------

ACTIVATE_DIR="${ENV_DIR}/etc/conda/activate.d"
DEACTIVATE_DIR="${ENV_DIR}/etc/conda/deactivate.d"

mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat > "${ACTIVATE_DIR}/oncodriveclustl_environment.sh" <<EOF2
#!/usr/bin/env bash
export BGDATA_LOCAL="${BGDATA_DIR}"
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
EOF2

cat > "${DEACTIVATE_DIR}/oncodriveclustl_environment.sh" <<'EOF2'
#!/usr/bin/env bash
unset BGDATA_LOCAL
unset PYTHONNOUSERSITE
unset MPLBACKEND
unset OMP_NUM_THREADS
unset MKL_NUM_THREADS
unset OPENBLAS_NUM_THREADS
unset NUMEXPR_NUM_THREADS
EOF2

chmod 755 \
    "${ACTIVATE_DIR}/oncodriveclustl_environment.sh" \
    "${DEACTIVATE_DIR}/oncodriveclustl_environment.sh"

# Apply settings to the current shell.
export BGDATA_LOCAL="$BGDATA_DIR"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

 echo
 echo "BGDATA_LOCAL=$BGDATA_LOCAL"

# ------------------------------------------------------------
# Preload and test the GRCh38 reference
# ------------------------------------------------------------
#
# OncodriveCLUSTL uses bgreference. The first access may download
# the GRCh38 reference into BGDATA_LOCAL. Running this during the
# installation avoids relying on internet access from a compute node.
# ------------------------------------------------------------

 echo
 echo "Testing/preloading the GRCh38 reference..."

python - <<'PY'
from bgreference import refseq

base = refseq("hg38", "1", 10001, 1)
if not base or base.upper() not in {"A", "C", "G", "T", "N"}:
    raise RuntimeError(f"Unexpected hg38 reference result: {base!r}")
print(f"PASS: hg38 chr1:10001 reference base = {base}")
PY

# ------------------------------------------------------------
# Create a reusable activation helper
# ------------------------------------------------------------

ACTIVATION_HELPER="${TOOL_DIR}/activate_oncodriveclustl.sh"

cat > "$ACTIVATION_HELPER" <<EOF2
#!/usr/bin/env bash

set -euo pipefail

if ! type module >/dev/null 2>&1; then
    if [[ -r /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
    elif [[ -r /usr/share/lmod/lmod/init/bash ]]; then
        source /usr/share/lmod/lmod/init/bash
    else
        echo "ERROR: The cluster module system could not be initialized." >&2
        return 1 2>/dev/null || exit 1
    fi
fi

module load "${ANACONDA_MODULE}"

CONDA_BASE="\$(conda info --base)"
source "\${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${ENV_DIR}"

export BGDATA_LOCAL="${BGDATA_DIR}"
export PYTHONNOUSERSITE=1
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

VERSION="\$(python -c 'from importlib.metadata import version; print(version("oncodriveclustl"))')"

echo "Environment: \$CONDA_PREFIX"
echo "OncodriveCLUSTL: \$(command -v oncodriveclustl)"
echo "OncodriveCLUSTL version: \$VERSION"
echo "BGDATA_LOCAL=\$BGDATA_LOCAL"
EOF2

chmod 755 "$ACTIVATION_HELPER"

# ------------------------------------------------------------
# Save reproducibility and diagnostic records
# ------------------------------------------------------------

oncodriveclustl --help \
    > "${RECORD_DIR}/oncodriveclustl_${ONCODRIVECLUSTL_VERSION}_help.txt"

python -m pip show \
    oncodriveclustl \
    bgparsers \
    bgreference \
    bgsignature \
    > "${RECORD_DIR}/oncodriveclustl_${ONCODRIVECLUSTL_VERSION}_package_details.txt"

python -m pip freeze \
    > "${RECORD_DIR}/pip_freeze_oncodriveclustl_${ONCODRIVECLUSTL_VERSION}.txt"

conda list --prefix "$ENV_DIR" \
    > "${RECORD_DIR}/conda_list_oncodriveclustl_${ONCODRIVECLUSTL_VERSION}.txt"

conda list --explicit --prefix "$ENV_DIR" \
    > "${RECORD_DIR}/conda_explicit_oncodriveclustl_${ONCODRIVECLUSTL_VERSION}.txt"

conda env export --prefix "$ENV_DIR" \
    > "${RECORD_DIR}/conda_environment_oncodriveclustl_${ONCODRIVECLUSTL_VERSION}.yml"

{
    echo "Installation date: $(date)"
    echo "Host: $(hostname)"
    echo "Anaconda module: $ANACONDA_MODULE"
    echo "Conda environment: $ENV_DIR"
    echo "BGDATA_LOCAL: $BGDATA_DIR"
    echo "Python: $(python --version 2>&1)"
    echo "OncodriveCLUSTL executable: $(command -v oncodriveclustl)"
    echo "OncodriveCLUSTL version: $DETECTED_VERSION"
} > "${RECORD_DIR}/installation_summary.txt"

PIP_CHECK_FILE="${RECORD_DIR}/pip_check_oncodriveclustl_${ONCODRIVECLUSTL_VERSION}.txt"

if python -m pip check > "$PIP_CHECK_FILE" 2>&1; then
    echo
    echo "pip check: no dependency problems reported."
else
    echo
    echo "WARNING: pip check reported dependency problems. Review:"
    echo "  $PIP_CHECK_FILE"
    exit 1
fi

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

 echo
 echo "============================================================"
 echo "OncodriveCLUSTL installation completed successfully."
 echo "Finished: $(date)"
 echo
 echo "Environment:"
 echo "  $ENV_DIR"
 echo
 echo "Executable:"
 echo "  $(command -v oncodriveclustl)"
 echo
 echo "Version:"
 echo "  $DETECTED_VERSION"
 echo
 echo "BGData directory:"
 echo "  $BGDATA_DIR"
 echo
 echo "Activation helper:"
 echo "  $ACTIVATION_HELPER"
 echo
 echo "Installation records:"
 echo "  $RECORD_DIR"
 echo
 echo "Installation log:"
 echo "  $INSTALL_LOG"
 echo "============================================================"
