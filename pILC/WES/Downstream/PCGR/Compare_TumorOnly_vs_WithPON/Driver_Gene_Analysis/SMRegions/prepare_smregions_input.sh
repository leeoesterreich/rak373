#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Prepare standalone SMRegions mutation input
#
# Source:
#   Combined PCGR MAF
#
# Output columns:
#   CHROMOSOME POSITION REF ALT SAMPLE
#
# Retained:
#   - simple SNVs only
#   - nonsynonymous protein-coding mutations
#
# Excluded:
#   - synonymous/Silent mutations
#   - indels and MNVs
#   - noncoding mutations
#   - mitochondrial mutations
#   - exact sample-level duplicates
#
# No IntOGen input or output is used.
# ============================================================

BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/SMRegions"

MAF="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Tier_1_common_all_samples.combined.maf"

INPUT_DIR="${BASE}/input"

OUTPUT="${INPUT_DIR}/pILC_WES_smregions_input.tsv.gz"
QC_OUTPUT="${INPUT_DIR}/pILC_WES_smregions_input_QC.tsv"

# Use another Python executable by exporting PYTHON before running,
# for example: export PYTHON=/path/to/python
PYTHON="${PYTHON:-python3}"

mkdir -p "${INPUT_DIR}"

# ------------------------------------------------------------
# Validate source and Python
# ------------------------------------------------------------

if [[ ! -s "${MAF}" ]]; then
    echo "ERROR: Combined PCGR MAF is missing or empty:"
    echo "       ${MAF}"
    exit 1
fi

if ! command -v "${PYTHON}" >/dev/null 2>&1; then
    echo "ERROR: Python executable was not found:"
    echo "       ${PYTHON}"
    exit 1
fi

echo "============================================================"
echo "Preparing standalone SMRegions input"
echo "============================================================"
echo "Source MAF: ${MAF}"
echo "Output:     ${OUTPUT}"
echo "QC file:    ${QC_OUTPUT}"
echo "Python:     $(command -v "${PYTHON}")"
echo "============================================================"

# ------------------------------------------------------------
# Prepare mutation file
# ------------------------------------------------------------

"${PYTHON}" - "${MAF}" "${OUTPUT}" "${QC_OUTPUT}" <<'PYTHON_SCRIPT'
import csv
import gzip
import os
import sys
from collections import Counter

maf_path = sys.argv[1]
output_path = sys.argv[2]
qc_path = sys.argv[3]

VALID_BASES = {"A", "C", "G", "T"}

# Protein-coding nonsynonymous classes relevant to an
# amino-acid/domain enrichment analysis.
#
# Indel classifications are listed here for documentation but
# are automatically removed by the single-base SNV requirement.
NONSYNONYMOUS_CLASSES = {
    "missense_mutation",
    "nonsense_mutation",
    "nonstop_mutation",
    "translation_start_site",
    "start_codon_snp",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "protein_altering_variant",
}

SYNONYMOUS_CLASSES = {
    "silent",
    "synonymous",
    "synonymous_variant",
}

stats = Counter()
samples = set()
seen = set()
observed_classifications = Counter()


def open_text(path):
    """Open a plain-text or gzip-compressed file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def locate_column(header, candidates, description):
    """Return the first matching MAF column."""
    for candidate in candidates:
        if candidate in header:
            return candidate

    raise RuntimeError(
        f"Unable to find the {description} column. "
        f"Tried: {', '.join(candidates)}"
    )


def normalize_chromosome(value):
    """Convert chr1/CHR1 to 1 and chrM/MT to M."""
    chromosome = str(value).strip()

    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]

    if chromosome.upper() in {"M", "MT"}:
        return "M"

    return chromosome


with open_text(maf_path) as handle:

    # --------------------------------------------------------
    # Find the MAF header, skipping comment/metadata lines.
    # --------------------------------------------------------

    header = None

    for line in handle:
        line = line.rstrip("\n\r")

        if not line:
            continue

        if line.startswith("#"):
            continue

        header = line.split("\t")
        break

    if header is None:
        raise RuntimeError("No header was found in the combined MAF.")

    chromosome_col = locate_column(
        header,
        ["Chromosome", "CHROMOSOME", "chr"],
        "chromosome",
    )

    position_col = locate_column(
        header,
        ["Start_Position", "POSITION", "Position", "pos"],
        "position",
    )

    ref_col = locate_column(
        header,
        ["Reference_Allele", "REF", "ref"],
        "reference allele",
    )

    alt_col = locate_column(
        header,
        [
            "Tumor_Seq_Allele2",
            "ALT",
            "alt",
            "mut",
        ],
        "alternate allele",
    )

    sample_col = locate_column(
        header,
        [
            "Tumor_Sample_Barcode",
            "SAMPLE",
            "Sample",
            "sampleID",
        ],
        "sample",
    )

    classification_col = locate_column(
        header,
        [
            "Variant_Classification",
            "Consequence",
            "CONSEQUENCE",
        ],
        "variant classification",
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

            stats["input_maf_rows"] += 1

            chromosome = normalize_chromosome(
                row.get(chromosome_col, "")
            )

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

            classification_original = str(
                row.get(classification_col, "")
            ).strip()

            classification = classification_original.lower()

            observed_classifications[classification_original] += 1

            # ------------------------------------------------
            # Required values
            # ------------------------------------------------

            if (
                not chromosome
                or not position_raw
                or not ref
                or not alt
                or not sample
            ):
                stats["excluded_missing_required_value"] += 1
                continue

            # ------------------------------------------------
            # Exclude mitochondrial variants
            # ------------------------------------------------

            if chromosome == "M":
                stats["excluded_mitochondrial"] += 1
                continue

            # ------------------------------------------------
            # Retain simple SNVs only
            # ------------------------------------------------

            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in VALID_BASES
                or alt not in VALID_BASES
            ):
                stats["excluded_non_snv"] += 1
                continue

            if ref == alt:
                stats["excluded_ref_equals_alt"] += 1
                continue

            # ------------------------------------------------
            # Parse position
            # ------------------------------------------------

            try:
                position = int(float(position_raw))
            except (TypeError, ValueError):
                stats["excluded_invalid_position"] += 1
                continue

            if position < 1:
                stats["excluded_invalid_position"] += 1
                continue

            # ------------------------------------------------
            # Exclude synonymous variants
            # ------------------------------------------------

            if classification in SYNONYMOUS_CLASSES:
                stats["excluded_synonymous"] += 1
                continue

            # ------------------------------------------------
            # Retain protein-changing coding substitutions
            # ------------------------------------------------

            if classification not in NONSYNONYMOUS_CLASSES:
                stats["excluded_other_classification"] += 1
                continue

            # ------------------------------------------------
            # Remove exact sample-level duplicates
            # ------------------------------------------------

            key = (
                chromosome,
                position,
                ref,
                alt,
                sample,
            )

            if key in seen:
                stats["excluded_exact_duplicate"] += 1
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

            stats["smregions_input_mutations"] += 1


stats["unique_samples"] = len(samples)

if stats["smregions_input_mutations"] == 0:

    if os.path.exists(output_path):
        os.remove(output_path)

    raise RuntimeError(
        "No eligible nonsynonymous coding SNVs were written. "
        "Review the Variant_Classification values in the QC file "
        "or combined MAF."
    )


# ------------------------------------------------------------
# Write QC report
# ------------------------------------------------------------

qc_metrics = [
    "input_maf_rows",
    "smregions_input_mutations",
    "unique_samples",
    "excluded_synonymous",
    "excluded_non_snv",
    "excluded_other_classification",
    "excluded_mitochondrial",
    "excluded_ref_equals_alt",
    "excluded_missing_required_value",
    "excluded_invalid_position",
    "excluded_exact_duplicate",
]

with open(qc_path, "wt") as qc_handle:

    qc_handle.write("Section\tMetric\tValue\n")

    for metric in qc_metrics:
        qc_handle.write(
            f"summary\t{metric}\t{stats[metric]}\n"
        )

    for classification, count in sorted(
        observed_classifications.items(),
        key=lambda item: (-item[1], item[0]),
    ):
        classification_label = (
            classification if classification else "MISSING"
        )

        qc_handle.write(
            "observed_variant_classification"
            f"\t{classification_label}\t{count}\n"
        )


print("SMRegions input generation completed successfully.")

for metric in qc_metrics:
    print(f"{metric}: {stats[metric]}")
PYTHON_SCRIPT

# ------------------------------------------------------------
# Validate final file
# ------------------------------------------------------------

if [[ ! -s "${OUTPUT}" ]]; then
    echo "ERROR: SMRegions input was not created."
    exit 1
fi

HEADER=$(zcat "${OUTPUT}" | head -n 1)

EXPECTED_HEADER=$'CHROMOSOME\tPOSITION\tREF\tALT\tSAMPLE'

if [[ "${HEADER}" != "${EXPECTED_HEADER}" ]]; then
    echo "ERROR: Unexpected SMRegions input header."
    echo "Observed: ${HEADER}"
    echo "Expected: ${EXPECTED_HEADER}"
    exit 1
fi

MUTATION_COUNT=$(zcat "${OUTPUT}" | awk 'NR > 1 {n++} END {print n+0}')
SAMPLE_COUNT=$(zcat "${OUTPUT}" |
    awk -F'\t' 'NR > 1 {samples[$5]=1}
        END {print length(samples)}')

echo
echo "============================================================"
echo "SMRegions input prepared successfully"
echo "============================================================"
echo "Input mutations: ${MUTATION_COUNT}"
echo "Unique samples:  ${SAMPLE_COUNT}"
echo "Output file:     ${OUTPUT}"
echo "QC file:         ${QC_OUTPUT}"
echo "============================================================"

echo
echo "Input preview:"
zcat "${OUTPUT}" | head -n 10

echo
echo "QC summary:"
awk -F'\t' '
    NR == 1 || $1 == "summary" {
        print
    }
' "${QC_OUTPUT}" | column -t -s $'\t'