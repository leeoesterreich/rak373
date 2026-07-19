#!/usr/bin/env bash

#SBATCH --job-name=Build_SMRegions_Resources
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/SMRegions/logs/Build_SMRegions_Resources_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

# ============================================================
# Build standalone GRCh38 resources required by SMRegions
#
# Permanent resources created/reused:
#   tools/SMRegions/resources/cds_biomart.tsv
#   tools/SMRegions/resources/pfam_biomart.tsv.gz
#   tools/SMRegions/resources/cds.counts.gz
#   tools/SMRegions/resources/regions_pfam.tsv
#
# Reused by direct path (not copied):
#   - canonical GRCh38 CDS regions
#   - existing BGData genome cache
#
# Build-only standalone prerequisites:
#   tools/SMRegions/build_resources/transvar.simg
#   tools/SMRegions/build_resources/transvar/
#
# No IntOGen workflow or IntOGen dataset path is used.
# ============================================================

PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/SMRegions"

TOOL_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/SMRegions"
ENV_DIR="${TOOL_DIR}/env_smregions"
RESOURCE_DIR="${TOOL_DIR}/resources"
BUILD_RESOURCE_DIR="${TOOL_DIR}/build_resources"

# Installed with standalone SMRegions.
PYTHON="${ENV_DIR}/bin/python"
BGSIGNATURE="${ENV_DIR}/bin/bgsignature"

# Permanent standalone SMRegions resources.
CDS_BIOMART="${RESOURCE_DIR}/cds_biomart.tsv"
PFAM_BIOMART="${RESOURCE_DIR}/pfam_biomart.tsv.gz"
CDS_COUNTS="${RESOURCE_DIR}/cds.counts.gz"
PFAM_OUTPUT="${RESOURCE_DIR}/regions_pfam.tsv"

# Reuse existing validated canonical resources without copying.
CDS_REGIONS="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveFML/resources/GRCh38/regions/cds.regions.gz"
BGDATA_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveCLUSTL/bgdata"

# Standalone build-only TransVar resources.
TRANSVAR_CONTAINER="${BUILD_RESOURCE_DIR}/transvar.simg"
TRANSVAR_DATA="${BUILD_RESOURCE_DIR}/transvar"

# Fixed Ensembl archive used for the validated resources.
BIOMART_URL="https://jan2024.archive.ensembl.org/biomart/martservice"

CORES="${SLURM_CPUS_PER_TASK:-8}"

EXPECTED_PFAM_HEADER=$'CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT_ID\tSEGMENT\tSYMBOL'

mkdir -p \
    "${RESOURCE_DIR}" \
    "${BUILD_RESOURCE_DIR}" \
    "${PROJECT_DIR}/logs"

module purge
module load singularity/4.3.2

echo "============================================================"
echo "Building standalone SMRegions resources"
echo "============================================================"
echo "Job ID:       ${SLURM_JOB_ID:-NA}"
echo "Host:         $(hostname)"
echo "Start:        $(date)"
echo "Resources:    ${RESOURCE_DIR}"
echo "Build inputs: ${BUILD_RESOURCE_DIR}"
echo "============================================================"

# ------------------------------------------------------------
# Validation functions
# ------------------------------------------------------------

validate_cds_biomart() {
    local file="$1"
    local rows
    local invalid

    [[ -s "${file}" ]] || return 1

    rows=$(
        awk 'END {print NR+0}' "${file}"
    )

    invalid=$(
        awk -F'\t' '
            NF != 13 ||
            $1 == "" ||
            $2 == "" ||
            $4 == "" ||
            $5 !~ /^[0-9]+$/ ||
            $6 !~ /^[0-9]+$/ ||
            $10 !~ /^-?1$/ ||
            $11 == "" {
                n++
            }
            END {
                print n+0
            }
        ' "${file}"
    )

    [[ "${rows}" -ge 100000 && "${invalid}" -eq 0 ]]
}

validate_pfam_biomart() {
    local file="$1"
    local rows
    local invalid

    [[ -s "${file}" ]] || return 1
    gzip -t "${file}" 2>/dev/null || return 1

    rows=$(
        gzip -cd "${file}" |
        awk 'END {print NR+0}'
    )

    invalid=$(
        gzip -cd "${file}" |
        awk -F'\t' '
            NF != 5 ||
            $1 == "" ||
            $2 == "" ||
            $3 !~ /^[0-9]+$/ ||
            $4 !~ /^[0-9]+$/ ||
            $5 == "" {
                n++
            }
            END {
                print n+0
            }
        '
    )

    [[ "${rows}" -ge 1000 && "${invalid}" -eq 0 ]]
}

validate_cds_counts() {
    local file="$1"

    [[ -s "${file}" ]] || return 1
    gzip -t "${file}" 2>/dev/null || return 1

    "${PYTHON}" - "${file}" <<'PYTHON' >/dev/null
import gzip
import json
import sys

path = sys.argv[1]

with gzip.open(path, "rt") as handle:
    counts = json.load(handle)

valid_bases = {"A", "C", "G", "T"}

if len(counts) != 64:
    raise RuntimeError(
        f"Expected 64 trinucleotide keys; observed {len(counts)}"
    )

for kmer, value in counts.items():
    if len(kmer) != 3 or set(kmer) - valid_bases:
        raise RuntimeError(f"Invalid trinucleotide key: {kmer}")

    if not isinstance(value, int) or value <= 0:
        raise RuntimeError(
            f"Invalid trinucleotide count for {kmer}: {value}"
        )
PYTHON
}

validate_pfam_output() {
    local file="$1"
    local header
    local rows
    local invalid

    [[ -s "${file}" ]] || return 1

    header=$(
        sed -n '1p' "${file}"
    )

    [[ "${header}" == "${EXPECTED_PFAM_HEADER}" ]] || return 1

    rows=$(
        awk 'NR > 1 {n++} END {print n+0}' "${file}"
    )

    invalid=$(
        awk -F'\t' '
            NR > 1 {
                if (
                    NF != 7 ||
                    $1 == "" ||
                    $2 !~ /^[0-9]+$/ ||
                    $3 !~ /^[0-9]+$/ ||
                    $2 > $3 ||
                    $5 == "" ||
                    $6 == "" ||
                    $7 == ""
                ) {
                    n++
                }
            }
            END {
                print n+0
            }
        ' "${file}"
    )

    [[ "${rows}" -ge 1000 && "${invalid}" -eq 0 ]]
}

# ------------------------------------------------------------
# Validate installed software and canonical resources
# ------------------------------------------------------------

for executable in \
    "${PYTHON}" \
    "${BGSIGNATURE}"
do
    if [[ ! -x "${executable}" ]]; then
        echo "ERROR: Required executable was not found:"
        echo "       ${executable}"
        exit 1
    fi
done

if [[ ! -s "${CDS_REGIONS}" ]]; then
    echo "ERROR: Canonical CDS regions file is missing:"
    echo "       ${CDS_REGIONS}"
    exit 1
fi

if [[ ! -d "${BGDATA_DIR}" ]]; then
    echo "ERROR: BGData directory is missing:"
    echo "       ${BGDATA_DIR}"
    exit 1
fi

export BGDATA_LOCAL="${BGDATA_DIR}"
export BGDATA_OFFLINE="TRUE"

# ------------------------------------------------------------
# Download the Ensembl 111 CDS BioMart table
# ------------------------------------------------------------

if validate_cds_biomart "${CDS_BIOMART}"; then

    CDS_BIOMART_ROWS=$(
        awk 'END {print NR+0}' "${CDS_BIOMART}"
    )

    echo
    echo "Reusing validated CDS BioMart table:"
    echo "${CDS_BIOMART}"
    echo "Rows: ${CDS_BIOMART_ROWS}"

elif [[ -e "${CDS_BIOMART}" ]]; then

    echo "ERROR: Existing CDS BioMart file is incomplete or malformed:"
    echo "       ${CDS_BIOMART}"
    echo "Remove it intentionally before rebuilding."
    exit 1

else

    echo
    echo "Downloading Ensembl 111 MANE CDS annotations..."

    TMP_CDS="${CDS_BIOMART}.tmp"

    rm -f "${TMP_CDS}"

    CDS_QUERY=$(
        cat <<'XML'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="chromosome_name" value="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
    <Filter name="biotype" value="protein_coding"/>
    <Filter name="mane_select" excluded="0"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="ensembl_peptide_id"/>
    <Attribute name="chromosome_name"/>
    <Attribute name="genomic_coding_start"/>
    <Attribute name="genomic_coding_end"/>
    <Attribute name="cds_start"/>
    <Attribute name="cds_end"/>
    <Attribute name="cds_length"/>
    <Attribute name="strand"/>
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="exon_chrom_start"/>
    <Attribute name="exon_chrom_end"/>
  </Dataset>
</Query>
XML
    )

    curl \
        --fail \
        --location \
        --silent \
        --show-error \
        --get \
        --data-urlencode "query=${CDS_QUERY}" \
        "${BIOMART_URL}" |
    tail -n +2 |
    awk -F'\t' '$5 != ""' \
        > "${TMP_CDS}"

    if ! validate_cds_biomart "${TMP_CDS}"; then
        echo "ERROR: Downloaded CDS BioMart table failed validation."
        rm -f "${TMP_CDS}"
        exit 1
    fi

    mv "${TMP_CDS}" "${CDS_BIOMART}"

    CDS_BIOMART_ROWS=$(
        awk 'END {print NR+0}' "${CDS_BIOMART}"
    )

    echo "CDS BioMart rows: ${CDS_BIOMART_ROWS}"

fi

# ------------------------------------------------------------
# Download the Ensembl 111 Pfam BioMart table
# ------------------------------------------------------------

if validate_pfam_biomart "${PFAM_BIOMART}"; then

    PFAM_BIOMART_ROWS=$(
        gzip -cd "${PFAM_BIOMART}" |
        awk 'END {print NR+0}'
    )

    echo
    echo "Reusing validated Pfam BioMart table:"
    echo "${PFAM_BIOMART}"
    echo "Rows: ${PFAM_BIOMART_ROWS}"

elif [[ -e "${PFAM_BIOMART}" ]]; then

    echo "ERROR: Existing Pfam BioMart file is incomplete or malformed:"
    echo "       ${PFAM_BIOMART}"
    echo "Remove it intentionally before rebuilding."
    exit 1

else

    echo
    echo "Downloading Ensembl 111 MANE/Pfam annotations..."

    TMP_PFAM="${PFAM_BIOMART}.tmp"

    rm -f "${TMP_PFAM}"

    PFAM_QUERY=$(
        cat <<'XML'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="chromosome_name" value="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
    <Filter name="biotype" value="protein_coding"/>
    <Filter name="mane_select" excluded="0"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="ensembl_transcript_id"/>
    <Attribute name="pfam_start"/>
    <Attribute name="pfam_end"/>
    <Attribute name="pfam"/>
  </Dataset>
</Query>
XML
    )

    curl \
        --fail \
        --location \
        --silent \
        --show-error \
        --get \
        --data-urlencode "query=${PFAM_QUERY}" \
        "${BIOMART_URL}" |
    tail -n +2 |
    awk -F'\t' '$4 != "" && $5 != ""' |
    gzip -c > "${TMP_PFAM}"

    if ! validate_pfam_biomart "${TMP_PFAM}"; then
        echo "ERROR: Downloaded Pfam BioMart table failed validation."
        rm -f "${TMP_PFAM}"
        exit 1
    fi

    mv "${TMP_PFAM}" "${PFAM_BIOMART}"

    PFAM_BIOMART_ROWS=$(
        gzip -cd "${PFAM_BIOMART}" |
        awk 'END {print NR+0}'
    )

    echo "Pfam BioMart rows: ${PFAM_BIOMART_ROWS}"

fi

# ------------------------------------------------------------
# Generate GRCh38 CDS trinucleotide counts
# ------------------------------------------------------------

if validate_cds_counts "${CDS_COUNTS}"; then

    echo
    echo "Reusing validated CDS trinucleotide counts:"
    echo "${CDS_COUNTS}"

elif [[ -e "${CDS_COUNTS}" ]]; then

    echo "ERROR: Existing CDS counts file is incomplete or malformed:"
    echo "       ${CDS_COUNTS}"
    echo "Remove it intentionally before rebuilding."
    exit 1

else

    echo
    echo "Generating GRCh38 CDS trinucleotide counts..."

    TMP_COUNTS="${CDS_COUNTS}.tmp.gz"

    rm -f "${TMP_COUNTS}"

    "${BGSIGNATURE}" count \
        -r "${CDS_REGIONS}" \
        -s 3 \
        -g hg38 \
        --collapse \
        --cores "${CORES}" \
        -o "${TMP_COUNTS}"

    if ! validate_cds_counts "${TMP_COUNTS}"; then
        echo "ERROR: Generated CDS counts failed validation."
        rm -f "${TMP_COUNTS}"
        exit 1
    fi

    mv "${TMP_COUNTS}" "${CDS_COUNTS}"

    echo "CDS trinucleotide counts generated:"
    echo "${CDS_COUNTS}"

fi

# ------------------------------------------------------------
# Reuse a validated completed Pfam genomic resource
# ------------------------------------------------------------

if validate_pfam_output "${PFAM_OUTPUT}"; then

    PFAM_ROWS=$(
        awk 'NR > 1 {n++} END {print n+0}' "${PFAM_OUTPUT}"
    )

    echo
    echo "Reusing validated GRCh38 Pfam genomic regions:"
    echo "${PFAM_OUTPUT}"
    echo "Mapped rows: ${PFAM_ROWS}"

    echo
    echo "============================================================"
    echo "All standalone SMRegions resources are ready"
    echo "============================================================"
    echo "CDS BioMart:   ${CDS_BIOMART}"
    echo "Pfam BioMart:  ${PFAM_BIOMART}"
    echo "CDS counts:    ${CDS_COUNTS}"
    echo "Pfam regions:  ${PFAM_OUTPUT}"
    echo "End:           $(date)"
    echo "============================================================"

    exit 0

elif [[ -e "${PFAM_OUTPUT}" ]]; then

    echo "ERROR: Existing Pfam genomic regions file is incomplete or malformed:"
    echo "       ${PFAM_OUTPUT}"
    echo "Remove it intentionally before rebuilding."
    exit 1

fi

# ------------------------------------------------------------
# Validate standalone TransVar build resources
# ------------------------------------------------------------

if [[ ! -s "${TRANSVAR_CONTAINER}" ]]; then
    echo "ERROR: Standalone TransVar container is missing:"
    echo "       ${TRANSVAR_CONTAINER}"
    exit 1
fi

if [[ ! -d "${TRANSVAR_DATA}" ]]; then
    echo "ERROR: Standalone TransVar data directory is missing:"
    echo "       ${TRANSVAR_DATA}"
    exit 1
fi

for file in \
    "${TRANSVAR_DATA}/Homo_sapiens.GRCh38.fa" \
    "${TRANSVAR_DATA}/Homo_sapiens.GRCh38.fa.fai" \
    "${TRANSVAR_DATA}/Homo_sapiens.GRCh38.111.gtf.gz" \
    "${TRANSVAR_DATA}/Homo_sapiens.GRCh38.111.gtf.gz.transvardb"
do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Required standalone TransVar resource is missing:"
        echo "       ${file}"
        exit 1
    fi
done

singularity run \
    -c \
    -B "${TRANSVAR_DATA}:/data" \
    "${TRANSVAR_CONTAINER}" \
    --help >/dev/null

echo
echo "Standalone TransVar build resources: PASSED"

# ------------------------------------------------------------
# Node-local temporary workspace
# ------------------------------------------------------------

TMP_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TMP_WORK="${TMP_ROOT}/smregions_resources_${SLURM_JOB_ID:-$$}"

mkdir -p "${TMP_WORK}"

cleanup() {
    rm -rf "${TMP_WORK}"
}

trap cleanup EXIT

TRANSCRIPTS_FILE="${TMP_WORK}/MANE_TRANSCRIPTS.tsv"
PFAM_TABLE="${TMP_WORK}/pfam_biomart.tsv"
WORKER="${TMP_WORK}/map_pfam_domain.sh"
TEMP_OUTPUT="${TMP_WORK}/regions_pfam.tsv"

# ------------------------------------------------------------
# Prepare transcript-to-gene/HUGO lookup
# ------------------------------------------------------------

awk -F'\t' '
    BEGIN {
        OFS = "\t"
    }
    {
        print $1, $11, $2
    }
' "${CDS_BIOMART}" |
sort -u > "${TRANSCRIPTS_FILE}"

if [[ ! -s "${TRANSCRIPTS_FILE}" ]]; then
    echo "ERROR: Transcript lookup could not be generated."
    exit 1
fi

# ------------------------------------------------------------
# Temporary TransVar worker
# ------------------------------------------------------------

cat > "${WORKER}" <<'WORKER_SCRIPT'
#!/usr/bin/env bash

set -euo pipefail

TRANSCRIPT="$1"
PFAM_DOMAIN="$2"
AA_START="$3"
AA_END="$4"
CONTAINER="$5"
DATA="$6"
TRANSCRIPTS_FILE="$7"

RANGE1=$(
    singularity run \
        -c \
        -B "${DATA}:/data" \
        "${CONTAINER}" \
        panno \
        --ensembl \
        -i "${TRANSCRIPT}:${AA_START}" |
    tail -n +2 |
    cut -f5 |
    cut -d'/' -f1 |
    sed 's/:g\./_/g' |
    tr '_' '\t' |
    sed -n '1p'
)

RANGE2=$(
    singularity run \
        -c \
        -B "${DATA}:/data" \
        "${CONTAINER}" \
        panno \
        --ensembl \
        -i "${TRANSCRIPT}:${AA_END}" |
    tail -n +2 |
    cut -f5 |
    cut -d'/' -f1 |
    sed 's/:g\./_/g' |
    tr '_' '\t' |
    sed -n '1p'
)

GENE_TRANSCRIPT_HUGO=$(
    awk -F'\t' -v transcript="${TRANSCRIPT}" '
        $2 == transcript {
            print
            exit
        }
    ' "${TRANSCRIPTS_FILE}"
)

HUGO=$(
    printf '%s\n' "${GENE_TRANSCRIPT_HUGO}" |
    cut -f3
)

CHROMOSOME=$(
    printf '%s\n' "${RANGE1}" |
    cut -f1 |
    sed 's/^chr//'
)

POSITIONS=$(
    printf '%s\t%s\n' "${RANGE1}" "${RANGE2}" |
    cut -f2,3,5,6
)

MIN_POSITION=$(
    printf '%s\n' "${POSITIONS}" |
    tr '\t' '\n' |
    awk '/^[0-9]+$/' |
    sort -n |
    sed -n '1p'
)

MAX_POSITION=$(
    printf '%s\n' "${POSITIONS}" |
    tr '\t' '\n' |
    awk '/^[0-9]+$/' |
    sort -nr |
    sed -n '1p'
)

if [[ -z "${CHROMOSOME}" ]] ||
   [[ -z "${MIN_POSITION}" ]] ||
   [[ -z "${MAX_POSITION}" ]] ||
   [[ -z "${HUGO}" ]]
then
    exit 0
fi

ELEMENT_ID="${TRANSCRIPT}:${PFAM_DOMAIN}:${AA_START}:${AA_END}"

printf '%s\t%s\t%s\t-\t%s\t%s\t%s\n' \
    "${CHROMOSOME}" \
    "${MIN_POSITION}" \
    "${MAX_POSITION}" \
    "${ELEMENT_ID}" \
    "${ELEMENT_ID}" \
    "${HUGO}"
WORKER_SCRIPT

chmod +x "${WORKER}"

# ------------------------------------------------------------
# Map Pfam amino-acid coordinates to GRCh38
# ------------------------------------------------------------

gzip -cd "${PFAM_BIOMART}" > "${PFAM_TABLE}"

TOTAL_ROWS=$(
    awk 'END {print NR+0}' "${PFAM_TABLE}"
)

CHUNK_SIZE=$(
    awk -v rows="${TOTAL_ROWS}" -v cores="${CORES}" '
        BEGIN {
            print int((rows + cores - 1) / cores)
        }
    '
)

split \
    -l "${CHUNK_SIZE}" \
    -d \
    -a 3 \
    "${PFAM_TABLE}" \
    "${TMP_WORK}/chunk."

printf '%s\n' "${EXPECTED_PFAM_HEADER}" > "${TEMP_OUTPUT}"

mapfile -t CHUNKS < <(
    find "${TMP_WORK}" \
        -maxdepth 1 \
        -type f \
        -name 'chunk.*' |
    sort
)

OUTPUT_PARTS=()

echo
echo "Mapping ${TOTAL_ROWS} Pfam rows using ${#CHUNKS[@]} workers..."

for chunk in "${CHUNKS[@]}"; do

    part="${chunk}.out"
    OUTPUT_PARTS+=("${part}")

    (
        while IFS=$'\t' read -r \
            ENSEMBL_GENE \
            TRANSCRIPT \
            AA_START \
            AA_END \
            PFAM_DOMAIN
        do
            if [[ -z "${TRANSCRIPT}" ]] ||
               [[ -z "${AA_START}" ]] ||
               [[ -z "${AA_END}" ]] ||
               [[ -z "${PFAM_DOMAIN}" ]]
            then
                continue
            fi

            "${WORKER}" \
                "${TRANSCRIPT}" \
                "${PFAM_DOMAIN}" \
                "${AA_START}" \
                "${AA_END}" \
                "${TRANSVAR_CONTAINER}" \
                "${TRANSVAR_DATA}" \
                "${TRANSCRIPTS_FILE}" ||
            true

        done < "${chunk}"
    ) > "${part}" &

done

wait

for part in "${OUTPUT_PARTS[@]}"; do
    cat "${part}" >> "${TEMP_OUTPUT}"
done

# ------------------------------------------------------------
# Validate and store the completed Pfam resource
# ------------------------------------------------------------

if ! validate_pfam_output "${TEMP_OUTPUT}"; then
    echo "ERROR: Generated Pfam genomic regions failed validation."
    exit 1
fi

PFAM_ROWS=$(
    awk 'NR > 1 {n++} END {print n+0}' "${TEMP_OUTPUT}"
)

mv "${TEMP_OUTPUT}" "${PFAM_OUTPUT}"

echo
echo "============================================================"
echo "Standalone SMRegions resources built successfully"
echo "============================================================"
echo "CDS BioMart rows:  ${CDS_BIOMART_ROWS}"
echo "Pfam BioMart rows: ${PFAM_BIOMART_ROWS}"
echo "Mapped Pfam rows:  ${PFAM_ROWS}"
echo "CDS BioMart:       ${CDS_BIOMART}"
echo "Pfam BioMart:      ${PFAM_BIOMART}"
echo "CDS counts:        ${CDS_COUNTS}"
echo "Pfam regions:      ${PFAM_OUTPUT}"
echo "End:               $(date)"
echo "============================================================"

echo
echo "Pfam preview:"
sed -n '1,5p' "${PFAM_OUTPUT}"
