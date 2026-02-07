#!/bin/bash
#SBATCH --job-name=Dedup
#SBATCH --output=slurm_log/Dedup/Dedup_%A_%a.out
#SBATCH --error=slurm_log/Dedup/Dedup_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7

set -euo pipefail

module load singularity/4.3.2

SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet.csv"
RG_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/3_BQSR_010526"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Dedup_010526"

SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUTPUT_DIR" slurm_log/Dedup

log_info(){ echo "[INFO] $*"; }
log_fail(){ echo "[ERROR] $*" >&2; exit 1; }

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$SHEET" | tr -d '\r')
[[ -n "${LINE:-}" ]] || { log_info "No line for task ${SLURM_ARRAY_TASK_ID}; exiting."; exit 0; }

SAMPLE=$(echo "$LINE" | cut -d',' -f1 | xargs)
[[ -n "${SAMPLE:-}" ]] || log_fail "Failed to parse SAMPLE from: $LINE"

INPUT_BAM="${RG_DIR}/${SAMPLE}/${SAMPLE}.rg.bam"
[[ -f "$INPUT_BAM" ]] || log_fail "RG BAM not found: $INPUT_BAM"
[[ -f "${INPUT_BAM%.bam}.bai" || -f "${INPUT_BAM}.bai" ]] || log_fail "RG BAM index not found for: $INPUT_BAM"

OUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

DEDUP_BAM="${OUT_DIR}/${SAMPLE}.dedup.bam"
METRICS="${OUT_DIR}/${SAMPLE}.dedup.metrics.txt"

if [[ -f "$DEDUP_BAM" && -f "${DEDUP_BAM}.bai" ]]; then
  log_info "$SAMPLE already processed — skipping."
  exit 0
fi

log_info "Sample: $SAMPLE"
log_info "Input RG BAM: $INPUT_BAM"
log_info "Output dedup BAM: $DEDUP_BAM"

singularity exec $BIND $SIF gatk MarkDuplicates \
  -I "$INPUT_BAM" \
  -O "$DEDUP_BAM" \
  -M "$METRICS" \
  --REMOVE_DUPLICATES true \
  --CREATE_INDEX true

log_info "[DONE] $SAMPLE"
