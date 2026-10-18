#!/bin/bash
#SBATCH --job-name=Check_integrity
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm_log/Check_%j.out
#SBATCH --error=slurm_log/Check_%j.err

set -uo pipefail

# -------------------- Directories --------------------
BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/Chelsea_21ILC_RNASeq_data"

FASTQ="${BASE}/Raw_Fastq"
META="${BASE}/metadata"
ACCESSIONS="${BASE}/SRR_accessions.txt"

mkdir -p "$META"

INTEGRITY_TSV="${META}/gzip_integrity.tsv"
READ_COUNTS_TSV="${META}/read_counts.tsv"
SUMMARY_TSV="${META}/integrity_summary.tsv"

echo "========================================"
echo "Chelsea 21 ILC RNA-seq: FASTQ QC"
echo "Started: $(date)"
echo "FASTQ directory: $FASTQ"
echo "========================================"

# -------------------- Validate input --------------------
if [[ ! -d "$FASTQ" ]]; then
    echo "ERROR: FASTQ directory does not exist."
    exit 1
fi

if [[ ! -s "$ACCESSIONS" ]]; then
    echo "ERROR: SRR accession list missing or empty."
    exit 1
fi

mapfile -t SAMPLES < <(grep -E '^SRR[0-9]+$' "$ACCESSIONS")

if [[ "${#SAMPLES[@]}" -ne 21 ]]; then
    echo "ERROR: Expected 21 SRR accessions."
    exit 1
fi

# -------------------- Initialize tables --------------------
printf "SRR\tRead\tStatus\tFile\n" > "$INTEGRITY_TSV"

printf "SRR\tR1_reads\tR2_reads\tpaired_match\n" \
    > "$READ_COUNTS_TSV"

printf "SRR\tR1_gzip\tR2_gzip\tR1_format\tR2_format\tPair_status\n" \
    > "$SUMMARY_TSV"

# -------------------- Function: Check FASTQ --------------------
check_fastq() {

    local file="$1"
    local srr="$2"
    local read="$3"

    if [[ ! -s "$file" ]]; then
        printf "%s\t%s\tMISSING\t%s\n" \
            "$srr" "$read" "$file" >> "$INTEGRITY_TSV"
        echo "MISSING"
        return
    fi

    if gzip -t "$file" 2>/dev/null; then
        printf "%s\t%s\tOK\t%s\n" \
            "$srr" "$read" "$file" >> "$INTEGRITY_TSV"
        echo "OK"
    else
        printf "%s\t%s\tCORRUPT\t%s\n" \
            "$srr" "$read" "$file" >> "$INTEGRITY_TSV"
        echo "CORRUPT"
    fi
}

# -------------------- Function: Count reads --------------------
count_reads() {

    local file="$1"

    if [[ ! -s "$file" ]]; then
        echo "NA"
        return
    fi

    local lines

    if ! lines=$(gzip -cd "$file" 2>/dev/null | wc -l); then
        echo "ERROR"
        return
    fi

    if (( lines % 4 != 0 )); then
        echo "INVALID"
        return
    fi

    echo $((lines / 4))
}

# -------------------- Process samples --------------------
for SRR in "${SAMPLES[@]}"; do

    echo "----------------------------------------"
    echo "Processing: $SRR"

    R1="${FASTQ}/${SRR}_1.fastq.gz"
    R2="${FASTQ}/${SRR}_2.fastq.gz"

    # 1. gzip checks
    G1=$(check_fastq "$R1" "$SRR" "R1")
    G2=$(check_fastq "$R2" "$SRR" "R2")

    # 2. Read counts
    if [[ "$G1" == "OK" ]]; then
        C1=$(count_reads "$R1")
    else
        C1="NA"
    fi

    if [[ "$G2" == "OK" ]]; then
        C2=$(count_reads "$R2")
    else
        C2="NA"
    fi

    # 3. Format status
    F1="INVALID"
    F2="INVALID"

    [[ "$C1" =~ ^[0-9]+$ ]] && F1="OK"
    [[ "$C2" =~ ^[0-9]+$ ]] && F2="OK"

    # 4. Paired read count comparison
    PAIR="FAIL"

    if [[ "$F1" == "OK" && "$F2" == "OK" ]]; then
        if [[ "$C1" -eq "$C2" ]]; then
            PAIR="YES"
        else
            PAIR="NO"
        fi
    fi

    # 5. Save results
    printf "%s\t%s\t%s\t%s\n" \
        "$SRR" "$C1" "$C2" "$PAIR" \
        >> "$READ_COUNTS_TSV"

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SRR" "$G1" "$G2" "$F1" "$F2" "$PAIR" \
        >> "$SUMMARY_TSV"

    echo "R1 reads: $C1"
    echo "R2 reads: $C2"
    echo "Pair status: $PAIR"

done

# -------------------- Overall summary --------------------
echo "========================================"
echo "FINAL SUMMARY"
echo "========================================"

TOTAL_SAMPLES="${#SAMPLES[@]}"

TOTAL_FASTQ=$(find "$FASTQ" -maxdepth 1 \
    -type f -name "*.fastq.gz" | wc -l)

PASS_SAMPLES=$(awk -F'\t' \
    'NR>1 && $6=="YES" {n++} END{print n+0}' \
    "$SUMMARY_TSV")

echo "Expected samples: 21"
echo "Samples checked: $TOTAL_SAMPLES"
echo "Expected FASTQ files: 42"
echo "FASTQ files found: $TOTAL_FASTQ"
echo "Samples passing basic checks: $PASS_SAMPLES"

echo ""
echo "Results saved to:"
echo "$INTEGRITY_TSV"
echo "$READ_COUNTS_TSV"
echo "$SUMMARY_TSV"

echo "Finished: $(date)"

if [[ "$PASS_SAMPLES" -ne 21 || "$TOTAL_FASTQ" -ne 42 ]]; then
    echo "WARNING: One or more integrity checks failed."
    exit 1
fi

echo "All basic integrity checks passed."
