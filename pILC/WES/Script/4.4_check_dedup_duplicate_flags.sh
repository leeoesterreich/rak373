#!/bin/bash
#SBATCH --job-name=check_dedup_dupflags
#SBATCH --output=slurm_log/dedup_check/dedup_check_%A_%a.out
#SBATCH --error=slurm_log/dedup_check/dedup_check_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-47

set -euo pipefail

module purge
module load samtools/1.22.1

DEDUP_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4.1_Dedup"
OUT_DIR="${DEDUP_DIR}/Duplicate_flag_check"

mkdir -p "$OUT_DIR/per_sample" "slurm_log/dedup_check"

# Only TP sample folders
mapfile -t sample_dirs < <(
    find "$DEDUP_DIR" -mindepth 1 -maxdepth 1 -type d -name "TP*" | sort
)

sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "${sample_dir}" || ! -d "$sample_dir" ]]; then
    echo "[ERROR] No valid sample directory found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

SAMPLE=$(basename "$sample_dir")
BAM="${sample_dir}/${SAMPLE}.dedup.bam"

if [[ ! -f "$BAM" ]]; then
    echo "[ERROR] Dedup BAM not found: $BAM"
    exit 1
fi

FLAGSTAT="${OUT_DIR}/per_sample/${SAMPLE}.flagstat.txt"
TSV="${OUT_DIR}/per_sample/${SAMPLE}.duplicate_check.tsv"

echo "[INFO] Sample: $SAMPLE"
echo "[INFO] BAM: $BAM"

samtools flagstat -@ "${SLURM_CPUS_PER_TASK}" "$BAM" > "$FLAGSTAT"

TOTAL_READS=$(awk '/in total/ {print $1; exit}' "$FLAGSTAT")
DUP_READS=$(awk '$0 ~ / duplicates$/ && $0 !~ /primary duplicates/ {print $1; exit}' "$FLAGSTAT")
PRIMARY_DUP_READS=$(awk '/primary duplicates/ {print $1; exit}' "$FLAGSTAT")

DUP_READS=${DUP_READS:-0}
PRIMARY_DUP_READS=${PRIMARY_DUP_READS:-0}

if [[ "$DUP_READS" -eq 0 && "$PRIMARY_DUP_READS" -eq 0 ]]; then
    STATUS="DUPLICATES_REMOVED"
else
    STATUS="DUPLICATES_REMAIN"
fi

echo -e "Sample\tBAM\tTotal_Reads\tDuplicate_Flagged_Reads\tPrimary_Duplicate_Reads\tStatus" > "$TSV"
echo -e "${SAMPLE}\t${BAM}\t${TOTAL_READS}\t${DUP_READS}\t${PRIMARY_DUP_READS}\t${STATUS}" >> "$TSV"

echo "[DONE] $SAMPLE $STATUS"