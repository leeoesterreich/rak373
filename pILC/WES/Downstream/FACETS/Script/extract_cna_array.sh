#!/bin/bash
#SBATCH --job-name=extract_cna
#SBATCH --output=slurm_log/extract_cna/extract_cna_%A_%a.out
#SBATCH --error=slurm_log/extract_cna/extract_cna_%A_%a.err
#SBATCH --array=0-15
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

set -euo pipefail

cd /ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS

module load gcc/12.2.0 r/4.4.0

export R_LIBS=/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4
export R_LIBS_USER=/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4

FILE_LIST="hisens_rds_files.txt"
R_SCRIPT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Script/extract_cna_from_facets_rds.R"
OUTDIR="cna_segs"

mkdir -p "$OUTDIR"
mkdir -p slurm_log/extract_cna

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
RDS_FILE=$(sed -n "${LINE_NUM}p" "$FILE_LIST")

if [[ -z "$RDS_FILE" ]]; then
    echo "ERROR: No RDS file found for array task $SLURM_ARRAY_TASK_ID, line $LINE_NUM" >&2
    exit 1
fi

if [[ ! -s "$RDS_FILE" ]]; then
    echo "ERROR: RDS file missing or empty: $RDS_FILE" >&2
    exit 1
fi

BASENAME=$(basename "$RDS_FILE" _hisens.rds)
OUT_FILE="${OUTDIR}/${BASENAME}_cna_segs.tsv"

echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Line number: $LINE_NUM"
echo "RDS file: $RDS_FILE"
echo "Output file: $OUT_FILE"
echo "R_LIBS: $R_LIBS"
echo "Start time: $(date)"

Rscript "$R_SCRIPT" "$RDS_FILE" "$OUT_FILE"

echo "Finished: $BASENAME"
echo "End time: $(date)"