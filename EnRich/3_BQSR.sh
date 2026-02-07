#!/bin/bash
#SBATCH --job-name=BQSR
#SBATCH --output=slurm_log/BQSR/BQSR_%A_%a.out
#SBATCH --error=slurm_log/BQSR/BQSR_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7

module load singularity/4.3.2

SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet.csv"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/3_BQSR_010526"

REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"
INDELS_GOLD="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
SNPS_1000G="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"

SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUTPUT_DIR" slurm_log/BQSR

# 0-based array: pick the (task_id+2)th line (skip header)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$SHEET" | tr -d '\r')

if [[ -z "$LINE" ]]; then
  echo "[SKIP] No line for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

SAMPLE=$(echo "$LINE" | cut -d',' -f1 | xargs)
INPUT_BAM=$(echo "$LINE" | cut -d',' -f2 | xargs)

if [[ -z "$SAMPLE" || -z "$INPUT_BAM" ]]; then
  echo "[ERROR] Could not parse SAMPLE/BAM from line: $LINE" >&2
  exit 1
fi

if [[ ! -f "$INPUT_BAM" ]]; then
  echo "[ERROR] BAM not found: $INPUT_BAM" >&2
  exit 1
fi

WORKDIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$WORKDIR"

TMP_RG_BAM="${WORKDIR}/${SAMPLE}.rg.bam"
MKDUP_BAM="${WORKDIR}/${SAMPLE}.mkdup.bam"
METRICS="${WORKDIR}/${SAMPLE}.mkdup.metrics.txt"
BQSR_TABLE="${WORKDIR}/${SAMPLE}.recal.table"
RECAL_BAM="${WORKDIR}/${SAMPLE}.bqsr.bam"

threads="${SLURM_CPUS_PER_TASK:-1}"

# ---- Read group: simple + consistent ----
RG_ID="$SAMPLE"
RG_PU="$SAMPLE"

echo "[INFO] Sample: $SAMPLE"
echo "[INFO] Input BAM: $INPUT_BAM"

# Step 1: AddOrReplaceReadGroups
singularity exec $BIND $SIF gatk AddOrReplaceReadGroups \
  -I "$INPUT_BAM" \
  -O "$TMP_RG_BAM" \
  -RGID "$RG_ID" \
  -RGSM "$SAMPLE" \
  -RGLB "${SAMPLE}_lib" \
  -RGPL ILLUMINA \
  -RGPU "$RG_PU" \
  --CREATE_INDEX true

# Step 2: MarkDuplicates
singularity exec $BIND $SIF gatk MarkDuplicates \
  -I "$TMP_RG_BAM" \
  -O "$MKDUP_BAM" \
  -M "$METRICS" \
  --CREATE_INDEX true

# Step 3: BaseRecalibrator
singularity exec $BIND $SIF gatk BaseRecalibrator \
  -R "$REF" \
  -I "$MKDUP_BAM" \
  --known-sites "$INDELS_GOLD" \
  --known-sites "$SNPS_1000G" \
  -O "$BQSR_TABLE"

# Step 4: ApplyBQSR
singularity exec $BIND $SIF gatk ApplyBQSR \
  -R "$REF" \
  -I "$MKDUP_BAM" \
  --bqsr-recal-file "$BQSR_TABLE" \
  -O "$RECAL_BAM"

echo "[DONE] BQSR completed for $SAMPLE"
