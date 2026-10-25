#!/bin/bash
# Build nf-core/rnafusion samplesheets (paired-end) for ONE batch of FASTQ files.
# Usage: bash make_samplesheet.sh            (uses the defaults below)
#        bash make_samplesheet.sh <FASTQ_DIR> <STRAND_SUMMARY> <OUT_NAME>
#        PREFIX=CHP_ bash make_samplesheet.sh      (add a prefix to every sample name)
# Edit the three defaults for each new batch. Never overwrites an existing sheet
# (set FORCE=1 to allow it).
# Outputs: <OUT_NAME>.csv (all samples) and <OUT_NAME>_test.csv (first 2 samples)

BASE=${BASE:-/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion}
FASTQ_DIR=${1:-/ix1/alee/LO_LAB/General/Lab_Data/20240315_Organoid_RNAseq/Organoid_mRNA_03_20_2024_Children_Hosp_RAW}
STRAND_SUMMARY=${2:-${BASE}/strand_check/strand_summary.tsv}
OUT_NAME=${3:-samplesheet_fusion}
# Sample-name prefix, e.g. CHP_ . Needed when names are only digits (108, 120 ...): nf-core/rnafusion
# reads them as integers and rejects them, because sample names must be text.
PREFIX=${PREFIX:-}
R1_SUFFIX=_R1_001.fastq.gz
R2_SUFFIX=_R2_001.fastq.gz

SHEET=${BASE}/${OUT_NAME}.csv
TEST=${BASE}/${OUT_NAME}_test.csv
OTHER_SHEET=${BASE}/samplesheet_fusion.csv      # earlier batch, used for the overlap check

# ---- Safety checks ----
[[ -d ${FASTQ_DIR} ]]       || { echo "ERROR: folder not found: ${FASTQ_DIR}" >&2; exit 1; }
[[ -f ${STRAND_SUMMARY} ]]  || { echo "ERROR: strand summary not found: ${STRAND_SUMMARY}" >&2; exit 1; }
if [[ -e ${SHEET} && "${FORCE:-0}" != "1" ]]; then
    echo "ERROR: ${SHEET} already exists. Choose another OUT_NAME or run with FORCE=1." >&2; exit 1
fi
mkdir -p ${BASE}/logs
shopt -s nullglob
R1_FILES=(${FASTQ_DIR}/*${R1_SUFFIX})
EXPECTED=${#R1_FILES[@]}
[[ ${EXPECTED} -gt 0 ]] || { echo "ERROR: no *${R1_SUFFIX} files in ${FASTQ_DIR}" >&2; exit 1; }

# ---- Header (rnafusion needs exactly these first 4 columns) ----
echo "sample,fastq_1,fastq_2,strandedness" > ${SHEET}

# ---- One row per R1 file ----
NUMERIC=0
for R1 in "${R1_FILES[@]}"; do
    R2=${R1/%${R1_SUFFIX}/${R2_SUFFIX}}
    if [[ ! -f "${R2}" ]]; then
        echo "WARNING: R2 missing for ${R1}" >&2; continue
    fi
    # Raw name: drop _S<number> and the suffix, hyphens -> underscores (matches the strand summary)
    RAW=$(basename "${R1}" ${R1_SUFFIX} | sed -E 's/_S[0-9]+$//' | tr '-' '_')
    SAMPLE=${PREFIX}${RAW}                       # name used in the samplesheet
    if [[ ${SAMPLE} =~ ^[0-9]+$ ]]; then
        echo "ERROR: sample name '${SAMPLE}' is only digits. Run with PREFIX=XXX_ (e.g. PREFIX=CHP_)." >&2
        NUMERIC=$((NUMERIC+1))
    fi
    # Strandedness from the Salmon summary (last column, looked up by raw name); unclear or missing -> auto
    CALL=$(awk -F'\t' -v s="${RAW}" '$1==s {print $NF}' ${STRAND_SUMMARY})
    case "${CALL}" in
        reverse) STRAND=reverse ;;
        forward) STRAND=forward ;;
        *)       STRAND=auto; echo "NOTE: ${SAMPLE} call='${CALL:-none}' -> auto" >&2 ;;
    esac
    echo "${SAMPLE},${R1},${R2},${STRAND}" >> ${SHEET}
done

# ---- Test sheet: header + first 2 samples ----
head -n 3 ${SHEET} > ${TEST}

# ---- Validation (printed only) ----
N=$(tail -n +2 ${SHEET} | wc -l)
echo "Samples in samplesheet: ${N} (R1 files found: ${EXPECTED})"
echo "Strandedness values:"
tail -n +2 ${SHEET} | cut -d, -f4 | sort | uniq -c

echo "Duplicated sample names (should print nothing):"
DUPS=$(tail -n +2 ${SHEET} | cut -d, -f1 | sort | uniq -d); echo "${DUPS}"

echo "R2 files without an R1 (should print nothing):"
for R2 in ${FASTQ_DIR}/*${R2_SUFFIX}; do
    [[ -f "${R2/%${R2_SUFFIX}/${R1_SUFFIX}}" ]] || echo "${R2}"
done

echo "FASTQ paths that do not exist (should print nothing):"
tail -n +2 ${SHEET} | cut -d, -f2,3 | tr ',' '\n' | while read -r P; do
    [[ -f "${P}" ]] || echo "${P}"
done

echo "Names already in ${OTHER_SHEET} (should print nothing):"
if [[ -f ${OTHER_SHEET} && ${OTHER_SHEET} != ${SHEET} ]]; then
    comm -12 <(tail -n +2 ${SHEET} | cut -d, -f1 | sort) \
             <(tail -n +2 ${OTHER_SHEET} | cut -d, -f1 | sort)
fi

echo "Test sheet:"; cat ${TEST}

# ---- Hard stop if something is wrong ----
if [[ "${N}" -ne "${EXPECTED}" || -n "${DUPS}" || "${NUMERIC}" -gt 0 ]]; then
    echo "ERROR: count mismatch, duplicate names or digit-only names. Do NOT submit the pipeline." >&2; exit 1
fi
echo "Samplesheets OK: ${SHEET} and ${TEST}"
