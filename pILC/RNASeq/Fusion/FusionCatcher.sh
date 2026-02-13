#!/usr/bin/env bash
#SBATCH --job-name=FusionCatcher
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
#SBATCH --array=0-30
#SBATCH --partition=htc
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/slurm_log/%x.%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/slurm_log/%x.%A_%a.err

set -euo pipefail

module load singularity/4.3.2

# -----------------------------
# Container
# -----------------------------
SIF="/software/rhel9/manual/install/fusioncatcher/fusioncatcher-1.33.sif"

# -----------------------------
# Inputs
# -----------------------------
SAMPLE_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/pILC_RNASeq_31samples.csv"
DB="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/FusionCatcher_db/current"

# Override broken faToTwoBit in SIF
UCSC_FATOTWOBIT="/ix1/alee/LO_LAB/Personal/Rahul/tools/ucsc/faToTwoBit"

# -----------------------------
# Outputs / working
# -----------------------------
BASE_OUTDIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/FusionCatcher"
WORK_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/FusionCatcher_work"
mkdir -p "$BASE_OUTDIR" "$WORK_BASE"

# -----------------------------
# Temp (node-local)
# -----------------------------
TMPDIR="/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

# -----------------------------
# Get sample for this task
# CSV header: Sample,R1_path,R2_path
# -----------------------------
row=$((SLURM_ARRAY_TASK_ID + 2))   # +2 because header is line 1
line="$(sed -n "${row}p" "${SAMPLE_SHEET}" || true)"
[[ -n "${line}" ]] || { echo "No line for task ${SLURM_ARRAY_TASK_ID} (row ${row})"; exit 1; }

IFS=',' read -r SAMPLE R1_PATH R2_PATH <<< "${line}"
SAMPLE="$(echo "$SAMPLE"   | tr -d '\r' | xargs)"
R1_PATH="$(echo "$R1_PATH" | tr -d '\r' | xargs)"
R2_PATH="$(echo "$R2_PATH" | tr -d '\r' | xargs)"

echo "[$(date)] TASK=${SLURM_ARRAY_TASK_ID} SAMPLE=${SAMPLE}"
echo "R1=$R1_PATH"
echo "R2=$R2_PATH"
echo "DB=$DB"
echo "TMPDIR=$TMPDIR"
echo "SIF=$SIF"
echo "UCSC_FATOTWOBIT=$UCSC_FATOTWOBIT"
echo

# -----------------------------
# Validate row completeness
# -----------------------------
[[ -n "$SAMPLE"  ]] || { echo "ERROR: empty Sample in CSV row ${row}: $line"; exit 10; }
[[ -n "$R1_PATH" ]] || { echo "ERROR: empty R1_path in CSV row ${row}: $line"; exit 11; }
[[ -n "$R2_PATH" ]] || { echo "ERROR: empty R2_path in CSV row ${row}: $line"; exit 12; }

# -----------------------------
# Fast checks
# -----------------------------
[[ -s "$R1_PATH" ]] || { echo "ERROR: R1 not found/empty: $R1_PATH"; exit 2; }
[[ -s "$R2_PATH" ]] || { echo "ERROR: R2 not found/empty: $R2_PATH"; exit 2; }
[[ -d "$DB"      ]] || { echo "ERROR: DB dir not found: $DB"; exit 3; }
[[ -s "$SIF"     ]] || { echo "ERROR: SIF not found: $SIF"; exit 4; }
[[ -x "$UCSC_FATOTWOBIT" ]] || { echo "ERROR: UCSC faToTwoBit not executable: $UCSC_FATOTWOBIT"; exit 5; }

# -----------------------------
# Prepare per-sample input directory (2 files only)
# -----------------------------
SAMPLE_INDIR="${WORK_BASE}/${SAMPLE}/input"
mkdir -p "$SAMPLE_INDIR"
ln -sfn "$R1_PATH" "${SAMPLE_INDIR}/${SAMPLE}_R1.fastq.gz"
ln -sfn "$R2_PATH" "${SAMPLE_INDIR}/${SAMPLE}_R2.fastq.gz"

# -----------------------------
# Output dir
# -----------------------------
SAMPLE_OUTDIR="${BASE_OUTDIR}/${SAMPLE}"
mkdir -p "$SAMPLE_OUTDIR"

echo "Input dir listing:"
ls -lah "$SAMPLE_INDIR"
echo

# -----------------------------
# Run FusionCatcher inside container
# -----------------------------
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "[$(date)] Running FusionCatcher..."
singularity exec --cleanenv \
  --bind /ix1:/ix1 \
  --bind /tmp:/tmp \
  --bind "$UCSC_FATOTWOBIT":/usr/local/bin/faToTwoBit \
  "$SIF" \
  fusioncatcher.py \
    -d "$DB" \
    -i "$SAMPLE_INDIR" \
    -o "$SAMPLE_OUTDIR" \
    -T "$TMPDIR" \
    -p "$THREADS" \
    -Z

echo "[$(date)] DONE: ${SAMPLE}"
