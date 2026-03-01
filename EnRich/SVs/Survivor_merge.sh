#!/bin/bash
# https://github.com/fritzsedlazeck/SURVIVOR
# Merge VCF files to generate a consensus file
# Params: max_dist=1000, min_support=2, use_type=1 (require same SVTYPE to merge), use_strand=1 (require same strand orientation), disabled=0, min_size=50

#SBATCH --job-name=Survivor_merge
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 1:00:00
#SBATCH --array=0-47
#SBATCH --output=slurm_log/Survivor/Survivor_%A_%a.out
#SBATCH --error=slurm_log/Survivor/Survivor_%A_%a.err

set -euo pipefail

# Make sure SLURM log directory exists
mkdir -p slurm_log/Survivor

# -----------------------------
# Modules
# -----------------------------
module load survivor/1.0.7
module load bcftools/1.22

# -----------------------------
# Inputs
# -----------------------------
DELLY_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20260213/Batch2_021326/Delly_021326"
MANTA_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20260213/Batch2_021326/Manta_021326"
LUMPY_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20260213/Batch2_021326/Lumpy_021326"

# -----------------------------
# Output base
# -----------------------------
SURVIVOR_BASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20260213/Batch2_021326/Survivor/"
VCF_U_DIR="${SURVIVOR_BASE}/vcf_uncompressed"
LIST_DIR="${SURVIVOR_BASE}/lists"
MERGE_DIR="${SURVIVOR_BASE}/merged"
LOG_DIR="${SURVIVOR_BASE}/logs"

mkdir -p "$SURVIVOR_BASE" "$VCF_U_DIR" "$LIST_DIR" "$MERGE_DIR" "$LOG_DIR"

echo "[INFO] Output base: $SURVIVOR_BASE"
echo "[INFO] Writing uncompressed VCFs to: $VCF_U_DIR"

# -----------------------------
# Build sample list + pick 1 sample for this array task
# -----------------------------
mapfile -t SAMPLES < <(ls -d "${DELLY_DIR}"/TP* 2>/dev/null | xargs -n1 basename | sort)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "[ERROR] No TP* sample folders found under: $DELLY_DIR"
  exit 1
fi

if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#SAMPLES[@]} ]]; then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} but only ${#SAMPLES[@]} samples found."
  exit 1
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "[INFO] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} -> SAMPLE=${SAMPLE}"

# -----------------------------
# Input files
# -----------------------------
DELLY_GZ="${DELLY_DIR}/${SAMPLE}/${SAMPLE}_delly_filtered_merged.vcf.gz"
MANTA_GZ="${MANTA_DIR}/${SAMPLE}/${SAMPLE}_manta_diploidSV.vcf.gz"
LUMPY_VCF="${LUMPY_DIR}/${SAMPLE}/${SAMPLE}.lumpy_filtered.vcf"

# Output copies (uncompressed)
DELLY_VCF_U="${VCF_U_DIR}/${SAMPLE}.delly.vcf"
MANTA_VCF_U="${VCF_U_DIR}/${SAMPLE}.manta.vcf"
LUMPY_VCF_U="${VCF_U_DIR}/${SAMPLE}.lumpy.vcf"

# Existence checks
[[ -f "$DELLY_GZ"  ]] || { echo "[ERROR] Missing Delly: $DELLY_GZ"; exit 1; }
[[ -f "$MANTA_GZ"  ]] || { echo "[ERROR] Missing Manta: $MANTA_GZ"; exit 1; }
[[ -f "$LUMPY_VCF" ]] || { echo "[ERROR] Missing Lumpy: $LUMPY_VCF"; exit 1; }

# -----------------------------
# 1) Make uncompressed copies (keep .vcf.gz)
# -----------------------------
if [[ ! -s "$DELLY_VCF_U" ]]; then
  echo "[INFO] Uncompressing Delly for $SAMPLE"
  bcftools view -Ov -o "$DELLY_VCF_U" "$DELLY_GZ"
fi

if [[ ! -s "$MANTA_VCF_U" ]]; then
  echo "[INFO] Uncompressing Manta for $SAMPLE"
  bcftools view -Ov -o "$MANTA_VCF_U" "$MANTA_GZ"
fi

if [[ ! -s "$LUMPY_VCF_U" ]]; then
  echo "[INFO] Copying Lumpy for $SAMPLE"
  cp -f "$LUMPY_VCF" "$LUMPY_VCF_U"
fi

# -----------------------------
# 2) Make SURVIVOR list file
# -----------------------------
LIST_FILE="${LIST_DIR}/${SAMPLE}_vcf_list.txt"
printf "%s\n%s\n%s\n" "$DELLY_VCF_U" "$MANTA_VCF_U" "$LUMPY_VCF_U" > "$LIST_FILE"

# -----------------------------
# 3) Merge
# -----------------------------
OUT_VCF="${MERGE_DIR}/${SAMPLE}_merged.vcf"
LOG="${LOG_DIR}/${SAMPLE}.merge.log"

echo "[INFO] SURVIVOR merge ${SAMPLE}"
SURVIVOR merge "$LIST_FILE" 1000 2 1 1 0 50 "$OUT_VCF" > "$LOG" 2>&1

# -----------------------------
# 4) Extract BND-like events
#   Keep records where:
#     - INFO contains SVTYPE=BND
#     - OR ALT contains '[' or ']'
# -----------------------------
BND_DIR="${SURVIVOR_BASE}/merged_bnd"
mkdir -p "$BND_DIR"

BND_VCF="${BND_DIR}/${SAMPLE}_merged.BND.vcf"

awk 'BEGIN{OFS="\t"}
     /^#/ {print; next}
     ($8 ~ /SVTYPE=BND/) || ($8 ~ /SVTYPE=TRA/) || ($5 ~ /\[/) || ($5 ~ /\]/) {print}
' "$OUT_VCF" > "$BND_VCF"

echo "[DONE] BND-only: $BND_VCF"
echo "[DONE] SAMPLE=${SAMPLE}"
echo "[DONE] List:   $LIST_FILE"
echo "[DONE] Merged: $OUT_VCF"
echo "[DONE] Log:    $LOG"
