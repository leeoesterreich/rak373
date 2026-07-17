#!/usr/bin/env bash

set -Eeuo pipefail

trap 'echo "ERROR: Installation failed at line ${LINENO}." >&2' ERR

# ============================================================
# Standalone OncodriveFML 2.5.0 installation
# Pitt CRC / Lmod + Anaconda
#
# This script:
#   1. Loads the cluster Anaconda module.
#   2. Creates an independent Python 3.10 Conda environment.
#   3. Installs htslib/tabix and OncodriveFML 2.5.0.
#   4. Creates a dedicated BGData directory.
#   5. Creates a reusable environment-activation helper.
#   6. Saves version, help, package, and environment records.
#
# Safe to rerun: an existing environment is reused and verified.
# This installs the software only. GRCh38 configuration, coding
# elements, reference genome, and score resources are separate.
# ============================================================

# ------------------------------------------------------------
# User-editable settings
# ------------------------------------------------------------

ANACONDA_MODULE="anaconda3/2023.09-0-python_3.11.5"

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveFML"
ENV_DIR="${TOOL_DIR}/env_oncodrivefml_2.5.0"
BGDATA_DIR="${TOOL_DIR}/bgdata"

ONCODRIVEFML_VERSION="2.5.0"
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
INSTALL_LOG="${LOG_DIR}/install_oncodrivefml_${RUN_DATE}.log"

exec > >(tee -a "$INSTALL_LOG") 2>&1

echo "============================================================"
echo "Standalone OncodriveFML installation"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Tool directory: $TOOL_DIR"
echo "Environment: $ENV_DIR"
echo "BGData directory: $BGDATA_DIR"
echo "Installation log: $INSTALL_LOG"
echo "============================================================"

# Avoid accidentally importing packages from ~/.local.
export PYTHONNOUSERSITE=1

# ------------------------------------------------------------
# Make the Lmod 'module' command available
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
        --channel bioconda \
        "python=${PYTHON_VERSION}" \
        pip \
        htslib \
        zlib \
        -y
else
    echo
    echo "Existing environment detected; creation skipped:"
    echo "  $ENV_DIR"
fi

# ------------------------------------------------------------
# Activate and validate the base environment
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

if ! command -v tabix >/dev/null 2>&1; then
    echo "ERROR: tabix is missing from the environment."
    exit 1
fi

echo "Tabix executable: $(command -v tabix)"
tabix --version

# ------------------------------------------------------------
# Install or verify OncodriveFML
# ------------------------------------------------------------

INSTALLED_VERSION="$(
    python - <<'PY' 2>/dev/null || true
try:
    from importlib.metadata import version
    print(version("oncodrivefml"))
except Exception:
    pass
PY
)"

if [[ "$INSTALLED_VERSION" == "$ONCODRIVEFML_VERSION" ]]; then
    echo
    echo "OncodriveFML ${ONCODRIVEFML_VERSION} is already installed."
else
    echo
    echo "Installing OncodriveFML ${ONCODRIVEFML_VERSION}..."

    python -m pip install \
        --no-cache-dir \
        --upgrade \
        "oncodrivefml==${ONCODRIVEFML_VERSION}"
fi

hash -r

# ------------------------------------------------------------
# Verify executable and version
# ------------------------------------------------------------

if ! command -v oncodrivefml >/dev/null 2>&1; then
    echo "ERROR: oncodrivefml executable was not found after installation."
    exit 1
fi

if [[ "$(command -v oncodrivefml)" != "${ENV_DIR}/bin/oncodrivefml" ]]; then
    echo "ERROR: OncodriveFML is not resolving from the intended environment."
    echo "Found: $(command -v oncodrivefml)"
    exit 1
fi

DETECTED_VERSION="$(oncodrivefml --version 2>&1)"

echo
echo "OncodriveFML executable: $(command -v oncodrivefml)"
echo "OncodriveFML version: $DETECTED_VERSION"

if [[ "$DETECTED_VERSION" != *"$ONCODRIVEFML_VERSION"* ]]; then
    echo "ERROR: Expected OncodriveFML ${ONCODRIVEFML_VERSION}, but found:"
    echo "  $DETECTED_VERSION"
    exit 1
fi

if ! command -v bg-data >/dev/null 2>&1; then
    echo "ERROR: bg-data executable was not installed."
    exit 1
fi

echo "BGData executable: $(command -v bg-data)"
bg-data --version || true

# ------------------------------------------------------------
# Configure BGDATA_LOCAL on future environment activation
# ------------------------------------------------------------

ACTIVATE_DIR="${ENV_DIR}/etc/conda/activate.d"
DEACTIVATE_DIR="${ENV_DIR}/etc/conda/deactivate.d"

mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

cat > "${ACTIVATE_DIR}/oncodrivefml_environment.sh" <<EOF
#!/usr/bin/env bash
export BGDATA_LOCAL="${BGDATA_DIR}"
export PYTHONNOUSERSITE=1
EOF

cat > "${DEACTIVATE_DIR}/oncodrivefml_environment.sh" <<'EOF'
#!/usr/bin/env bash
unset BGDATA_LOCAL
unset PYTHONNOUSERSITE
EOF

chmod 755 \
    "${ACTIVATE_DIR}/oncodrivefml_environment.sh" \
    "${DEACTIVATE_DIR}/oncodrivefml_environment.sh"

# Apply it to the current installation shell.
export BGDATA_LOCAL="$BGDATA_DIR"

echo
echo "BGDATA_LOCAL=$BGDATA_LOCAL"

# ------------------------------------------------------------
# Create a reusable activation helper
# ------------------------------------------------------------

ACTIVATION_HELPER="${TOOL_DIR}/activate_oncodrivefml.sh"

cat > "$ACTIVATION_HELPER" <<EOF
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

echo "Environment: \$CONDA_PREFIX"
echo "OncodriveFML: \$(command -v oncodrivefml)"
oncodrivefml --version
echo "Tabix: \$(command -v tabix)"
tabix --version
echo "BGDATA_LOCAL=\$BGDATA_LOCAL"
EOF

chmod 755 "$ACTIVATION_HELPER"

# ------------------------------------------------------------
# Save reproducibility and diagnostic records
# ------------------------------------------------------------

oncodrivefml --help \
    > "${RECORD_DIR}/oncodrivefml_${ONCODRIVEFML_VERSION}_help.txt"

python -m pip show \
    oncodrivefml \
    bgdata \
    bgreference \
    bgsignature \
    pytabix \
    > "${RECORD_DIR}/oncodrivefml_${ONCODRIVEFML_VERSION}_package_details.txt"

python -m pip freeze \
    > "${RECORD_DIR}/pip_freeze_oncodrivefml_${ONCODRIVEFML_VERSION}.txt"

conda list --prefix "$ENV_DIR" \
    > "${RECORD_DIR}/conda_list_oncodrivefml_${ONCODRIVEFML_VERSION}.txt"

conda list --explicit --prefix "$ENV_DIR" \
    > "${RECORD_DIR}/conda_explicit_oncodrivefml_${ONCODRIVEFML_VERSION}.txt"

conda env export --prefix "$ENV_DIR" \
    > "${RECORD_DIR}/conda_environment_oncodrivefml_${ONCODRIVEFML_VERSION}.yml"

{
    echo "Installation date: $(date)"
    echo "Host: $(hostname)"
    echo "Anaconda module: $ANACONDA_MODULE"
    echo "Conda environment: $ENV_DIR"
    echo "BGDATA_LOCAL: $BGDATA_DIR"
    echo "Python: $(python --version 2>&1)"
    echo "OncodriveFML executable: $(command -v oncodrivefml)"
    echo "OncodriveFML version: $(oncodrivefml --version 2>&1)"
    echo "Tabix executable: $(command -v tabix)"
    echo "Tabix version: $(tabix --version 2>&1 | head -n 1)"
    echo "BGData executable: $(command -v bg-data)"
} > "${RECORD_DIR}/installation_summary.txt"

# pip check may report unrelated packages already present in an
# environment. Save the report, but do not invalidate a working
# OncodriveFML installation solely because of those warnings.
PIP_CHECK_FILE="${RECORD_DIR}/pip_check_oncodrivefml_${ONCODRIVEFML_VERSION}.txt"

if python -m pip check > "$PIP_CHECK_FILE" 2>&1; then
    echo
    echo "pip check: no dependency problems reported."
else
    echo
    echo "WARNING: pip check reported dependency warnings."
    echo "Review:"
    echo "  $PIP_CHECK_FILE"
    echo "The core OncodriveFML executable and version checks passed."
fi

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "============================================================"
echo "OncodriveFML installation completed successfully."
echo "Finished: $(date)"
echo
echo "Environment:"
echo "  $ENV_DIR"
echo
echo "Executable:"
echo "  $(command -v oncodrivefml)"
echo
echo "Version:"
echo "  $(oncodrivefml --version 2>&1)"
echo
echo "Tabix:"
echo "  $(command -v tabix)"
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
echo
echo "To activate in a future login shell, run:"
echo "  source \"$ACTIVATION_HELPER\""
echo
echo "NOTE:"
echo "  This confirms software installation only."
echo "  GRCh38 configuration, genomic elements, reference sequence,"
echo "  and compatible functional-impact score data must still be"
echo "  prepared before running the pILC WES analysis."
echo "============================================================"
