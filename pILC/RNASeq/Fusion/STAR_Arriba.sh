#!/bin/bash
#SBATCH --job-name=STAR_ArribaMode
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -t 2-00:00
#SBATCH --array=0-30
#SBATCH --partition=htc
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/logs/%x.%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/logs/%x.%A_%a.err

set -euo pipefail

module purge
module load star/2.7.11b

SAMPLE_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/pILC_RNASeq_31samples.csv"
STAR_INDEX="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7"
GTF_FILE="${STAR_INDEX}/gencode.v47.basic.annotation.gtf"
OUT_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/Alignment_arribaSTAR"

mkdir -p "${OUT_BASE}"
mkdir -p "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/logs"

THREADS="${SLURM_CPUS_PER_TASK}"

row=$((SLURM_ARRAY_TASK_ID + 2))
line=$(sed -n "${row}p" "${SAMPLE_SHEET}" || true)
[[ -n "${line}" ]] || { echo "No line for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

IFS=',' read -r SAMPLE R1 R2 <<< "${line}"
SAMPLE=$(echo "$SAMPLE" | tr -d '\r' | xargs)
R1=$(echo "$R1" | tr -d '\r' | xargs)
R2=$(echo "$R2" | tr -d '\r' | xargs)

SAMPLE_DIR="${OUT_BASE}/${SAMPLE}"
mkdir -p "${SAMPLE_DIR}"

OUT_PREFIX="${SAMPLE_DIR}/${SAMPLE}."
LOG_FILE="${SAMPLE_DIR}/${SAMPLE}.STAR.log"

[[ -f "${R1}" ]] || { echo "R1 not found: ${R1}" | tee "${LOG_FILE}"; exit 1; }
[[ -f "${R2}" ]] || { echo "R2 not found: ${R2}" | tee -a "${LOG_FILE}"; exit 1; }

echo "[$(date)] STAR Arriba-mode: ${SAMPLE}" | tee "${LOG_FILE}"

STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" --genomeLoad NoSharedMemory \
  --sjdbGTFfile "${GTF_FILE}" \
  --readFilesIn "${R1}" "${R2}" --readFilesCommand zcat \
  --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
  --outFileNamePrefix "${OUT_PREFIX}" \
  --outFilterMultimapNmax 50 \
  --alignSplicedMateMapLminOverLmate 0.5 \
  --alignSJstitchMismatchNmax 5 -1 5 5 \
  --chimSegmentMin 10 \
  --chimOutType WithinBAM HardClip \
  --chimJunctionOverhangMin 10 \
  --chimScoreDropMax 30 \
  --chimScoreJunctionNonGTAG 0 \
  --chimScoreSeparation 1 \
  --chimSegmentReadGapMax 3 \
  --chimMultimapNmax 50 \
  >> "${LOG_FILE}" 2>&1

[[ -s "${OUT_PREFIX}Aligned.out.bam" ]] || { echo "ERROR: Missing Aligned.out.bam for ${SAMPLE}" | tee -a "${LOG_FILE}"; exit 1; }

echo "[$(date)] Done STAR Arriba-mode: ${SAMPLE}" | tee -a "${LOG_FILE}"
