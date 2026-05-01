#!/bin/bash
#SBATCH --job-name=Dedup
#SBATCH --output=slurm_log/Dedup/Dedup_%A_%a.out
#SBATCH --error=slurm_log/Dedup/Dedup_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-47

set -euo pipefail

# Load required modules
module purge
module load singularity/4.3.2 samtools/1.22.1

# Paths
ALIGN_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/3_Alignment"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4.1_Dedup"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUTPUT_DIR" "slurm_log/Dedup"

# === Get sample folder from alignment directory ===
mapfile -t sample_dirs < <(find "$ALIGN_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "Qualimap" | sort)

sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "${sample_dir}" || ! -d "${sample_dir}" ]]; then
    echo "[ERROR] No valid sample directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

SAMPLE=$(basename "$sample_dir")

INPUT_BAM="${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam"
WORKDIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$WORKDIR"

TMP_RG_BAM="${WORKDIR}/${SAMPLE}.rg.bam"
DEDUP_BAM="${WORKDIR}/${SAMPLE}.dedup.bam"
METRICS="${WORKDIR}/${SAMPLE}.dedup.metrics.txt"

# === Check input BAM ===
if [[ ! -f "$INPUT_BAM" ]]; then
    echo "[ERROR] Input BAM not found for $SAMPLE: $INPUT_BAM"
    exit 1
fi

if [[ ! -f "${INPUT_BAM}.bai" ]]; then
    echo "[INFO] BAM index missing for $SAMPLE; creating index"
    samtools index "$INPUT_BAM"
fi


# === Read Group values ===
# Since this workflow merged lanes before trimming/alignment, use sample-level RG values
RG_ID="${SAMPLE}.RG1"
RG_PU="${SAMPLE}.unit1"

# === Step 1: AddOrReplaceReadGroups ===
echo "[INFO] Adding Read Groups for $SAMPLE"
singularity exec $BIND $SIF \
    gatk AddOrReplaceReadGroups \
    -I "$INPUT_BAM" \
    -O "$TMP_RG_BAM" \
    -RGID "$RG_ID" \
    -RGSM "$SAMPLE" \
    -RGLB "${SAMPLE}_lib" \
    -RGPL ILLUMINA \
    -RGPU "$RG_PU" \
    --CREATE_INDEX true

# === Step 2: Remove Duplicates ===
echo "[INFO] Removing duplicates for $SAMPLE"
singularity exec $BIND $SIF \
    gatk MarkDuplicates \
    -I "$TMP_RG_BAM" \
    -O "$DEDUP_BAM" \
    -M "$METRICS" \
    --REMOVE_DUPLICATES true \
    --CREATE_INDEX true

# === Cleanup of intermediate RG BAM and index ===
rm -f "$TMP_RG_BAM" "${TMP_RG_BAM%.bam}.bai"

echo "[DONE] AddOrReplaceReadGroups, Remove Duplicates completed for $SAMPLE"