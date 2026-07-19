#!/usr/bin/env bash

#SBATCH --job-name=SMRegions_pILC
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/SMRegions/logs/SMRegions_pILC_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00

set -euo pipefail

# ============================================================
# Standalone SMRegions
#
# Uses:
#   - one SMRegions-specific nonsynonymous SNV input
#   - combined PCGR MAF for cohort-signature generation
#   - existing canonical reference resources by direct path
#
# No reference or mutation files are copied.
# No IntOGen workflow output or intermediate file is used.
# ============================================================

PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/SMRegions"

ENV_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/SMRegions/env_smregions"

mkdir -p \
    "${PROJECT_DIR}/config" \
    "${PROJECT_DIR}/input" \
    "${PROJECT_DIR}/results" \
    "${PROJECT_DIR}/logs"

# ------------------------------------------------------------
# Environment and software
# ------------------------------------------------------------

module purge
module load anaconda3

SMREGIONS="${ENV_DIR}/bin/smregions"
BGSIGNATURE="${ENV_DIR}/bin/bgsignature"
PYTHON="${ENV_DIR}/bin/python"

# ------------------------------------------------------------
# Mutation inputs
# ------------------------------------------------------------

# Nonsynonymous coding SNVs prepared for SMRegions
SMREGIONS_INPUT="${PROJECT_DIR}/input/pILC_WES_smregions_input.tsv.gz"

# Canonical combined PCGR MAF
# A temporary all-SNV input is created only when the cohort
# signature does not already exist.
COMBINED_MAF="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Tier_1_common_all_samples.combined.maf"

# ------------------------------------------------------------
# Existing canonical reference resources
# ------------------------------------------------------------

CDS_REGIONS="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveFML/resources/GRCh38/regions/cds.regions.gz"

SMREGIONS_RESOURCES="/ix1/alee/LO_LAB/Personal/Rahul/tools/SMRegions/resources"

CDS_COUNTS="${SMREGIONS_RESOURCES}/cds.counts.gz"

PFAM_REGIONS="${SMREGIONS_RESOURCES}/regions_pfam.tsv"

BGDATA_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveCLUSTL/bgdata"

# CDS_COUNTS and PFAM_REGIONS are static reference resources.
# They are not cohort-specific IntOGen workflow outputs.

# ------------------------------------------------------------
# SMRegions files
# ------------------------------------------------------------

CONFIG="${PROJECT_DIR}/config/smregions.conf"

SIGNATURE="${PROJECT_DIR}/input/pILC_WES_cds_signature.json"

OUTPUT="${PROJECT_DIR}/results/pILC_WES_smregions_results.tsv.gz"

EXPECTED_SMREGIONS_MUTATIONS=22412
EXPECTED_ALL_SNVS=37662
EXPECTED_SAMPLES=31

echo "============================================================"
echo "Standalone SMRegions"
echo "============================================================"
echo "Job ID:  ${SLURM_JOB_ID:-NA}"
echo "Host:    $(hostname)"
echo "Start:   $(date)"
echo "Input:   ${SMREGIONS_INPUT}"
echo "Output:  ${OUTPUT}"
echo "============================================================"

# ------------------------------------------------------------
# Validate software
# ------------------------------------------------------------

for EXECUTABLE in \
    "${SMREGIONS}" \
    "${BGSIGNATURE}" \
    "${PYTHON}"
do
    if [[ ! -x "${EXECUTABLE}" ]]; then
        echo "ERROR: Executable was not found:"
        echo "       ${EXECUTABLE}"
        exit 1
    fi
done

"${SMREGIONS}" --version

# ------------------------------------------------------------
# Validate inputs and references
# ------------------------------------------------------------

for FILE in \
    "${SMREGIONS_INPUT}" \
    "${COMBINED_MAF}" \
    "${CDS_REGIONS}" \
    "${CDS_COUNTS}" \
    "${PFAM_REGIONS}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required file is missing or empty:"
        echo "       ${FILE}"
        exit 1
    fi
done

if [[ ! -d "${BGDATA_DIR}" ]]; then
    echo "ERROR: BGDATA directory was not found:"
    echo "       ${BGDATA_DIR}"
    exit 1
fi

# ------------------------------------------------------------
# Validate SMRegions-specific mutation input
# ------------------------------------------------------------

if ! gzip -t "${SMREGIONS_INPUT}"; then
    echo "ERROR: SMRegions input is not a valid gzip file:"
    echo "       ${SMREGIONS_INPUT}"
    exit 1
fi

SMREGIONS_HEADER=$(
    gzip -cd "${SMREGIONS_INPUT}" |
    sed -n '1p'
)

EXPECTED_HEADER=$'CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE'

if [[ "${SMREGIONS_HEADER}" != "${EXPECTED_HEADER}" ]]; then
    echo "ERROR: Unexpected SMRegions input header."
    echo "Expected: ${EXPECTED_HEADER}"
    echo "Observed: ${SMREGIONS_HEADER}"
    exit 1
fi

SMREGIONS_MUTATIONS=$(
    gzip -cd "${SMREGIONS_INPUT}" |
    awk 'NR > 1 {n++} END {print n+0}'
)

SMREGIONS_SAMPLES=$(
    gzip -cd "${SMREGIONS_INPUT}" |
    awk -F'\t' '
        NR > 1 {
            samples[$5] = 1
        }
        END {
            print length(samples)
        }
    '
)

if [[ "${SMREGIONS_MUTATIONS}" -ne "${EXPECTED_SMREGIONS_MUTATIONS}" ]]; then
    echo "ERROR: Unexpected SMRegions mutation count."
    echo "Expected: ${EXPECTED_SMREGIONS_MUTATIONS}"
    echo "Observed: ${SMREGIONS_MUTATIONS}"
    exit 1
fi

if [[ "${SMREGIONS_SAMPLES}" -ne "${EXPECTED_SAMPLES}" ]]; then
    echo "ERROR: Unexpected SMRegions sample count."
    echo "Expected: ${EXPECTED_SAMPLES}"
    echo "Observed: ${SMREGIONS_SAMPLES}"
    exit 1
fi

echo
echo "SMRegions input: ${SMREGIONS_MUTATIONS} mutations"
echo "Samples:         ${SMREGIONS_SAMPLES}"
echo "Signature source: combined PCGR MAF"

# ------------------------------------------------------------
# Create configuration once
# ------------------------------------------------------------

if [[ ! -s "${CONFIG}" ]]; then
    cat > "${CONFIG}" <<'EOF'
reference_genome = 'hg38'
muts_min = 3
sampling = 1000
sampling_chunk = 100
EOF
fi

echo
echo "Configuration:"
cat "${CONFIG}"

# ------------------------------------------------------------
# Configure existing bgdata cache
# ------------------------------------------------------------

export BGDATA_LOCAL="${BGDATA_DIR}"
export BGDATA_OFFLINE="TRUE"

# ------------------------------------------------------------
# Generate cohort signature once
# ------------------------------------------------------------

if [[ ! -s "${SIGNATURE}" ]]; then

    echo
    echo "Preparing temporary all-SNV signature input..."

    TMP_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"

    TMP_WORK="${TMP_ROOT}/smregions_${SLURM_JOB_ID:-$$}"

    mkdir -p "${TMP_WORK}"

    ALL_SNV_INPUT="${TMP_WORK}/pILC_WES_all_SNVs.tsv.gz"

    cleanup() {
        rm -rf "${TMP_WORK}"
    }

    trap cleanup EXIT

    # --------------------------------------------------------
    # Extract all simple SNVs from the plain-text combined MAF
    # --------------------------------------------------------

    "${PYTHON}" - \
        "${COMBINED_MAF}" \
        "${ALL_SNV_INPUT}" <<'PYTHON'

import csv
import gzip
import os
import sys

maf_path = sys.argv[1]
output_path = sys.argv[2]

valid_bases = {"A", "C", "G", "T"}

seen = set()
samples = set()

input_rows = 0
written = 0

excluded_non_snv = 0
excluded_mitochondrial = 0
excluded_invalid = 0
excluded_duplicate = 0


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "rt")


def find_column(header, candidates, description):
    for candidate in candidates:
        if candidate in header:
            return candidate

    raise RuntimeError(
        f"Could not find the {description} column. "
        f"Tried: {', '.join(candidates)}"
    )


with open_text(maf_path) as handle:

    header = None

    # Skip MAF metadata/comment lines and locate header.
    for line in handle:
        line = line.rstrip("\n\r")

        if not line or line.startswith("#"):
            continue

        header = line.split("\t")
        break

    if header is None:
        raise RuntimeError(
            "No header was found in the combined MAF."
        )

    chromosome_col = find_column(
        header,
        ["Chromosome", "CHROMOSOME", "chr"],
        "chromosome",
    )

    position_col = find_column(
        header,
        ["Start_Position", "POSITION", "Position", "pos"],
        "position",
    )

    ref_col = find_column(
        header,
        ["Reference_Allele", "REF", "ref"],
        "reference allele",
    )

    alt_col = find_column(
        header,
        [
            "Tumor_Seq_Allele2",
            "ALT",
            "alt",
            "mut",
        ],
        "alternate allele",
    )

    sample_col = find_column(
        header,
        [
            "Tumor_Sample_Barcode",
            "SAMPLE",
            "Sample",
            "sampleID",
        ],
        "sample",
    )

    reader = csv.DictReader(
        handle,
        fieldnames=header,
        delimiter="\t",
    )

    with gzip.open(output_path, "wt") as output_handle:

        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "CHROMOSOME",
            "POSITION",
            "REF",
            "ALT",
            "SAMPLE",
        ])

        for row in reader:

            input_rows += 1

            chromosome = str(
                row.get(chromosome_col, "")
            ).strip()

            if chromosome.lower().startswith("chr"):
                chromosome = chromosome[3:]

            if chromosome.upper() in {"M", "MT"}:
                excluded_mitochondrial += 1
                continue

            position_raw = str(
                row.get(position_col, "")
            ).strip()

            ref = str(
                row.get(ref_col, "")
            ).strip().upper()

            alt = str(
                row.get(alt_col, "")
            ).strip().upper()

            sample = str(
                row.get(sample_col, "")
            ).strip()

            # Retain all simple SNVs regardless of functional
            # classification.
            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in valid_bases
                or alt not in valid_bases
                or ref == alt
            ):
                excluded_non_snv += 1
                continue

            try:
                position = int(float(position_raw))
            except (TypeError, ValueError):
                excluded_invalid += 1
                continue

            if (
                position < 1
                or not chromosome
                or not sample
            ):
                excluded_invalid += 1
                continue

            key = (
                chromosome,
                position,
                ref,
                alt,
                sample,
            )

            if key in seen:
                excluded_duplicate += 1
                continue

            seen.add(key)
            samples.add(sample)

            writer.writerow([
                chromosome,
                position,
                ref,
                alt,
                sample,
            ])

            written += 1


if written == 0:
    if os.path.exists(output_path):
        os.remove(output_path)

    raise RuntimeError(
        "No valid SNVs were generated."
    )


print(f"Combined MAF rows:       {input_rows}")
print(f"All-SNV mutations:       {written}")
print(f"Unique samples:          {len(samples)}")
print(f"Excluded non-SNV rows:   {excluded_non_snv}")
print(f"Excluded mitochondrial:  {excluded_mitochondrial}")
print(f"Excluded invalid rows:   {excluded_invalid}")
print(f"Excluded duplicates:     {excluded_duplicate}")

PYTHON

    # --------------------------------------------------------
    # Validate temporary all-SNV input
    # --------------------------------------------------------

    if ! gzip -t "${ALL_SNV_INPUT}"; then
        echo "ERROR: Temporary all-SNV input is not valid gzip."
        exit 1
    fi

    ALL_SNV_COUNT=$(
        gzip -cd "${ALL_SNV_INPUT}" |
        awk 'NR > 1 {n++} END {print n+0}'
    )

    ALL_SNV_SAMPLES=$(
        gzip -cd "${ALL_SNV_INPUT}" |
        awk -F'\t' '
            NR > 1 {
                samples[$5] = 1
            }
            END {
                print length(samples)
            }
        '
    )

    echo
    echo "Temporary signature input:"
    echo "Mutations: ${ALL_SNV_COUNT}"
    echo "Samples:   ${ALL_SNV_SAMPLES}"

    if [[ "${ALL_SNV_COUNT}" -ne "${EXPECTED_ALL_SNVS}" ]]; then
        echo "ERROR: Unexpected all-SNV mutation count."
        echo "Expected: ${EXPECTED_ALL_SNVS}"
        echo "Observed: ${ALL_SNV_COUNT}"
        exit 1
    fi

    if [[ "${ALL_SNV_SAMPLES}" -ne "${EXPECTED_SAMPLES}" ]]; then
        echo "ERROR: Unexpected all-SNV sample count."
        echo "Expected: ${EXPECTED_SAMPLES}"
        echo "Observed: ${ALL_SNV_SAMPLES}"
        exit 1
    fi

    # --------------------------------------------------------
    # Generate cohort trinucleotide signature
    # --------------------------------------------------------

    echo
    echo "Generating cohort trinucleotide signature..."

    "${BGSIGNATURE}" normalize \
        -m "${ALL_SNV_INPUT}" \
        -r "${CDS_REGIONS}" \
        --normalize "${CDS_COUNTS}" \
        -s 3 \
        -g hg38 \
        --collapse \
        --cores "${SLURM_CPUS_PER_TASK:-8}" \
        -o "${SIGNATURE}"

else

    echo
    echo "Reusing existing signature:"
    echo "${SIGNATURE}"

fi

# ------------------------------------------------------------
# Validate signature
# ------------------------------------------------------------

if [[ ! -s "${SIGNATURE}" ]]; then
    echo "ERROR: Signature generation failed."
    exit 1
fi

"${PYTHON}" - "${SIGNATURE}" <<'PYTHON'

import json
import sys

signature_path = sys.argv[1]

with open(signature_path) as handle:
    signature = json.load(handle)

if not signature:
    raise RuntimeError(
        "Signature JSON is empty."
    )

print("Signature validation: PASSED")

PYTHON

# ------------------------------------------------------------
# Protect existing completed output
# ------------------------------------------------------------

if [[ -s "${OUTPUT}" ]]; then

    if gzip -t "${OUTPUT}" 2>/dev/null; then
        echo
        echo "Completed output already exists:"
        echo "${OUTPUT}"
        echo "No rerun was performed."
        exit 0
    fi

    echo
    echo "Removing invalid partial output:"
    echo "${OUTPUT}"

    rm -f "${OUTPUT}"

elif [[ -e "${OUTPUT}" ]]; then

    # Remove an empty file left by an interrupted run.
    rm -f "${OUTPUT}"

fi

# ------------------------------------------------------------
# Run SMRegions
# ------------------------------------------------------------

echo
echo "Running SMRegions..."

"${SMREGIONS}" \
    -m "${SMREGIONS_INPUT}" \
    -e "${CDS_REGIONS}" \
    -r "${PFAM_REGIONS}" \
    -s "${SIGNATURE}" \
    -c "${CONFIG}" \
    --cores "${SLURM_CPUS_PER_TASK:-8}" \
    --seed 12345 \
    -o "${OUTPUT}"

# ------------------------------------------------------------
# Validate final output
# ------------------------------------------------------------

if [[ ! -s "${OUTPUT}" ]]; then
    echo "ERROR: SMRegions output was not generated."
    exit 1
fi

if ! gzip -t "${OUTPUT}"; then
    echo "ERROR: SMRegions output is not valid gzip."
    exit 1
fi

RESULT_ROWS=$(
    gzip -cd "${OUTPUT}" |
    awk 'NR > 1 {n++} END {print n+0}'
)

echo
echo "============================================================"
echo "SMRegions completed successfully"
echo "============================================================"
echo "Input mutations: ${SMREGIONS_MUTATIONS}"
echo "Samples:         ${SMREGIONS_SAMPLES}"
echo "Result rows:     ${RESULT_ROWS}"
echo "Output:          ${OUTPUT}"
echo "End:             $(date)"
echo "============================================================"

echo
echo "Result preview:"

gzip -cd "${OUTPUT}" |
sed -n '1,10p'