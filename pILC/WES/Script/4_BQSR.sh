#!/bin/bash
#SBATCH --job-name=BQSR
#SBATCH --output=slurm_log/BQSR/BQSR_%A_%a.out
#SBATCH --error=slurm_log/BQSR/BQSR_%A_%a.err
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
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4_BQSR"
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
INDELS_GOLD="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
SNPS_1000G="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUTPUT_DIR" "slurm_log/BQSR"

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
BQSR_TABLE="${WORKDIR}/${SAMPLE}.recal.table"
RECAL_BAM="${WORKDIR}/${SAMPLE}.bqsr.bam"

# === Check input BAM ===
if [[ ! -f "$INPUT_BAM" ]]; then
    echo "[ERROR] Input BAM not found for $SAMPLE: $INPUT_BAM"
    exit 1
fi

if [[ ! -f "${INPUT_BAM}.bai" ]]; then
    echo "[INFO] BAM index missing for $SAMPLE; creating index"
    samtools index "$INPUT_BAM"
fi

# === Optional skip if final outputs already exist ===
if [[ -f "$RECAL_BAM" && -f "${RECAL_BAM}.bai" && -f "$BQSR_TABLE" && -f "$METRICS" ]]; then
    echo "[SKIP] BQSR already completed for $SAMPLE"
    exit 0
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

# === Step 2: MarkDuplicates ===
echo "[INFO] Marking duplicates for $SAMPLE"
singularity exec $BIND $SIF \
    gatk MarkDuplicates \
    -I "$TMP_RG_BAM" \
    -O "$DEDUP_BAM" \
    -M "$METRICS" \
    --CREATE_INDEX true

# === Step 3: BaseRecalibrator ===
echo "[INFO] Running BaseRecalibrator for $SAMPLE"
singularity exec $BIND $SIF \
    gatk BaseRecalibrator \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --known-sites "$INDELS_GOLD" \
    --known-sites "$SNPS_1000G" \
    -O "$BQSR_TABLE"

# === Step 4: ApplyBQSR ===
echo "[INFO] Applying BQSR for $SAMPLE"
singularity exec $BIND $SIF \
    gatk ApplyBQSR \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --bqsr-recal-file "$BQSR_TABLE" \
    -O "$RECAL_BAM"

# === Index final BAM ===
samtools index "$RECAL_BAM"

# === Optional cleanup of intermediate RG BAM ===
rm -f "$TMP_RG_BAM" "${TMP_RG_BAM}.bai"

echo "[DONE] AddOrReplaceReadGroups, MarkDuplicates, BaseRecalibrator, and ApplyBQSR completed for $SAMPLE"