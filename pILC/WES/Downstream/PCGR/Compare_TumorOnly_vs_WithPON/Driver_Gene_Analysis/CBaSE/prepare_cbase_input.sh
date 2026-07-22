#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Prepare standalone CBaSE v1.2 input from the canonical PCGR MAF
#
# CBaSE VCF-mode input columns:
#   CHROM POS ID REF ALT SAMPLE_ID
#
# Retained:
#   - all simple somatic SNVs, including synonymous variants
#   - canonical chromosomes 1-22, X, and Y
#
# Excluded:
#   - indels and MNVs
#   - mitochondrial/noncanonical contigs
#   - invalid rows
#   - exact sample-level duplicates
#
# Synonymous mutations must remain because CBaSE uses them to
# estimate the neutral mutation background.
# ============================================================

PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/Driver_Gene_Analysis/CBaSE"

COMBINED_MAF="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Tier_1_common_all_samples.combined.maf"

INPUT_DIR="${PROJECT_DIR}/input"
OUTPUT="${INPUT_DIR}/pILC_WES_cbase_input.vcf"

EXPECTED_MUTATIONS=37662
EXPECTED_SAMPLES=31

PYTHON="${PYTHON:-python3}"

mkdir -p "${INPUT_DIR}"

if [[ ! -s "${COMBINED_MAF}" ]]; then
    echo "ERROR: Combined PCGR MAF is missing or empty:"
    echo "       ${COMBINED_MAF}"
    exit 1
fi

if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "ERROR: Python executable was not found: ${PYTHON}"
    exit 1
fi

if [[ -s "${OUTPUT}" ]]; then
    HEADER="$(sed -n '1p' "${OUTPUT}")"
    EXPECTED_HEADER=$'#CHROM\tPOS\tID\tREF\tALT\tSAMPLE_ID'

    if [[ "${HEADER}" == "${EXPECTED_HEADER}" ]]; then
        MUTATION_COUNT="$(awk 'NR > 1 {n++} END {print n+0}' "${OUTPUT}")"
        SAMPLE_COUNT="$(
            awk -F'\t' 'NR > 1 {samples[$6]=1} END {print length(samples)}' \
                "${OUTPUT}"
        )"

        if [[ "${MUTATION_COUNT}" -eq "${EXPECTED_MUTATIONS}" ]] && \
           [[ "${SAMPLE_COUNT}" -eq "${EXPECTED_SAMPLES}" ]]; then
            echo "Valid completed CBaSE input already exists:"
            echo "${OUTPUT}"
            echo "Mutations: ${MUTATION_COUNT}"
            echo "Samples:   ${SAMPLE_COUNT}"
            exit 0
        fi
    fi

    echo "ERROR: Existing CBaSE input is incomplete or unexpected:"
    echo "       ${OUTPUT}"
    echo "Review or remove only this file before rerunning."
    exit 1
fi

echo "============================================================"
echo "Preparing standalone CBaSE input"
echo "============================================================"
echo "Source MAF: ${COMBINED_MAF}"
echo "Output:     ${OUTPUT}"
echo "============================================================"

TEMP_OUTPUT="${OUTPUT}.tmp.$$"

cleanup() {
    rm -f "${TEMP_OUTPUT}"
}

trap cleanup EXIT

"${PYTHON}" - "${COMBINED_MAF}" "${TEMP_OUTPUT}" <<'PYTHON_SCRIPT'
import csv
import gzip
import sys

maf_path = sys.argv[1]
output_path = sys.argv[2]

valid_bases = {"A", "C", "G", "T"}
valid_chromosomes = {str(i) for i in range(1, 23)} | {"X", "Y"}

seen = set()
samples = set()

stats = {
    "input_maf_rows": 0,
    "written_snvs": 0,
    "excluded_non_snv": 0,
    "excluded_noncanonical_chromosome": 0,
    "excluded_invalid": 0,
    "excluded_duplicate": 0,
}


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

    for line in handle:
        line = line.rstrip("\n\r")
        if not line or line.startswith("#"):
            continue
        header = line.split("\t")
        break

    if header is None:
        raise RuntimeError("No header was found in the combined MAF.")

    chromosome_col = find_column(
        header, ["Chromosome", "CHROMOSOME", "chr"], "chromosome"
    )
    position_col = find_column(
        header, ["Start_Position", "POSITION", "Position", "pos"], "position"
    )
    ref_col = find_column(
        header, ["Reference_Allele", "REF", "ref"], "reference allele"
    )
    alt_col = find_column(
        header,
        ["Tumor_Seq_Allele2", "ALT", "alt", "mut"],
        "alternate allele",
    )
    sample_col = find_column(
        header,
        ["Tumor_Sample_Barcode", "SAMPLE", "Sample", "sampleID"],
        "sample",
    )

    reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")

    with open(output_path, "wt", newline="") as output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "SAMPLE_ID"])

        for row in reader:
            stats["input_maf_rows"] += 1

            chromosome = str(row.get(chromosome_col, "")).strip()
            if chromosome.lower().startswith("chr"):
                chromosome = chromosome[3:]
            chromosome = chromosome.upper()

            if chromosome not in valid_chromosomes:
                stats["excluded_noncanonical_chromosome"] += 1
                continue

            ref = str(row.get(ref_col, "")).strip().upper()
            alt = str(row.get(alt_col, "")).strip().upper()
            sample = str(row.get(sample_col, "")).strip()
            position_raw = str(row.get(position_col, "")).strip()

            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in valid_bases
                or alt not in valid_bases
                or ref == alt
            ):
                stats["excluded_non_snv"] += 1
                continue

            try:
                position = int(float(position_raw))
            except (TypeError, ValueError):
                stats["excluded_invalid"] += 1
                continue

            if position < 1 or not sample:
                stats["excluded_invalid"] += 1
                continue

            key = (chromosome, position, ref, alt, sample)
            if key in seen:
                stats["excluded_duplicate"] += 1
                continue

            seen.add(key)
            samples.add(sample)
            writer.writerow([chromosome, position, ".", ref, alt, sample])
            stats["written_snvs"] += 1

if stats["written_snvs"] == 0:
    raise RuntimeError("No valid SNVs were generated for CBaSE.")

for metric, value in stats.items():
    print(f"{metric}: {value}")
print(f"unique_samples: {len(samples)}")
PYTHON_SCRIPT

EXPECTED_HEADER=$'#CHROM\tPOS\tID\tREF\tALT\tSAMPLE_ID'
OBSERVED_HEADER="$(sed -n '1p' "${TEMP_OUTPUT}")"

if [[ "${OBSERVED_HEADER}" != "${EXPECTED_HEADER}" ]]; then
    echo "ERROR: Unexpected CBaSE input header."
    echo "Observed: ${OBSERVED_HEADER}"
    exit 1
fi

MUTATION_COUNT="$(awk 'NR > 1 {n++} END {print n+0}' "${TEMP_OUTPUT}")"
SAMPLE_COUNT="$(
    awk -F'\t' 'NR > 1 {samples[$6]=1} END {print length(samples)}' \
        "${TEMP_OUTPUT}"
)"

if [[ "${MUTATION_COUNT}" -ne "${EXPECTED_MUTATIONS}" ]]; then
    echo "ERROR: Unexpected CBaSE SNV count."
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

mv "${TEMP_OUTPUT}" "${OUTPUT}"
trap - EXIT

echo
echo "============================================================"
echo "CBaSE input prepared successfully"
echo "============================================================"
echo "Mutations: ${MUTATION_COUNT}"
echo "Samples:   ${SAMPLE_COUNT}"
echo "Output:    ${OUTPUT}"
echo "============================================================"

echo
echo "Input preview:"
sed -n '1,10p' "${OUTPUT}"
