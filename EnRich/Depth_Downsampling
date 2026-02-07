#!/bin/bash
# Targeted Depth Normalization
#SBATCH --job-name=RAW_Downsample
#SBATCH --output=slurm_log/Downsample/RAW_DS_%A_%a.out
#SBATCH --error=slurm_log/Downsample/RAW_DS_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-7

set -euo pipefail

module purge
module load samtools/1.22.1

# =========================
# Paths
# =========================
SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet_012425.csv"
OUTDIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Downsampling/012425_RAW"
mkdir -p "$OUTDIR" slurm_log/Downsample

# =========================
# Panel parameters
# =========================
# Panel BED total size (bp): PIK3CA + ESR1 + GATA3 + FOXA1 + TP53/ERBB2 + CDH1
PANEL_BP=264569

# Mean read length (bp)
READ_LEN=150

# Target panel coverages
TARGET_COV_20K=20000
TARGET_COV_10K=10000

# Compute target reads
TARGET_READS_20K=$(echo "scale=2; ($TARGET_COV_20K * $PANEL_BP) / $READ_LEN" | bc)
TARGET_READS_10K=$(echo "scale=2; ($TARGET_COV_10K * $PANEL_BP) / $READ_LEN" | bc)

SEED=42

# =========================
# Parse sample sheet
# =========================
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$SHEET" | tr -d '\r')
[[ -z "$LINE" ]] && exit 0

SAMPLE=$(echo "$LINE" | cut -d',' -f1 | xargs)
INPUT_BAM=$(echo "$LINE" | cut -d',' -f2 | xargs)

if [[ ! -f "$INPUT_BAM" ]]; then
  echo "[ERROR] BAM not found: $INPUT_BAM"
  exit 1
fi

SAMPLE_DIR="${OUTDIR}/${SAMPLE}"
mkdir -p "$SAMPLE_DIR"

CSV="${SAMPLE_DIR}/${SAMPLE}.raw_downsampling_panel.csv"
echo "Sample,Total_Reads,Target_Reads_20K,Fraction_20K,Target_Reads_10K,Fraction_10K" > "$CSV"

# =========================
# Total reads in BAM
# =========================
TOTAL_READS=$(samtools view -c "$INPUT_BAM")

echo "[INFO] Sample: $SAMPLE"
echo "[INFO] Total reads: $TOTAL_READS"
echo "[INFO] Target reads (20K): $TARGET_READS_20K"
echo "[INFO] Target reads (10K): $TARGET_READS_10K"

# =========================
# Fractions
# =========================
FRAC_20K=$(echo "scale=6; $TARGET_READS_20K / $TOTAL_READS" | bc)
FRAC_10K=$(echo "scale=6; $TARGET_READS_10K / $TOTAL_READS" | bc)

echo "$SAMPLE,$TOTAL_READS,$TARGET_READS_20K,$FRAC_20K,$TARGET_READS_10K,$FRAC_10K" >> "$CSV"

# =========================
# Downsampling
# =========================
mk_seed_frac () {
  local seed="$1"
  local frac="$2"

  # normalize: remove spaces, ensure we don't have leading "."
  frac="${frac%% *}"
  frac="${frac#.}"          # if ".091" -> "091"
  frac="${frac#0.}"         # if "0.091" -> "091"
  echo "${seed}.${frac}"
}

echo "[INFO] Downsampling 20K panel depth"
if awk -v f="$FRAC_20K" 'BEGIN{exit !(f>=0.999999)}'; then
  cp -f "$INPUT_BAM" "${SAMPLE_DIR}/${SAMPLE}.20K.sorted.bam"
  [[ -f "${INPUT_BAM}.bai" ]] && cp -f "${INPUT_BAM}.bai" "${SAMPLE_DIR}/${SAMPLE}.20K.sorted.bam.bai" || samtools index "${SAMPLE_DIR}/${SAMPLE}.20K.sorted.bam"
else
  SARG=$(mk_seed_frac "$SEED" "$FRAC_20K")
  samtools view -@ "$SLURM_CPUS_PER_TASK" -b -s "$SARG" "$INPUT_BAM" > "${SAMPLE_DIR}/${SAMPLE}.20K.sorted.bam"
  samtools index "${SAMPLE_DIR}/${SAMPLE}.20K.sorted.bam"
fi

echo "[INFO] Downsampling 10K panel depth"
if awk -v f="$FRAC_10K" 'BEGIN{exit !(f>=0.999999)}'; then
  cp -f "$INPUT_BAM" "${SAMPLE_DIR}/${SAMPLE}.10K.sorted.bam"
  [[ -f "${INPUT_BAM}.bai" ]] && cp -f "${INPUT_BAM}.bai" "${SAMPLE_DIR}/${SAMPLE}.10K.sorted.bam.bai" || samtools index "${SAMPLE_DIR}/${SAMPLE}.10K.sorted.bam"
else
  SARG=$(mk_seed_frac "$SEED" "$FRAC_10K")
  samtools view -@ "$SLURM_CPUS_PER_TASK" -b -s "$SARG" "$INPUT_BAM" > "${SAMPLE_DIR}/${SAMPLE}.10K.sorted.bam"
  samtools index "${SAMPLE_DIR}/${SAMPLE}.10K.sorted.bam"
fi

echo "[DONE] Raw downsampling completed for $SAMPLE"
