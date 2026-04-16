#!/bin/bash
#SBATCH --job-name=BQSR
#SBATCH --output=slurm_log/BQSR/BQSR_%A_%a.out
#SBATCH --error=slurm_log/BQSR/BQSR_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1

# Load required module
module load singularity/4.3.2

# Paths
SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/3_Alignment/Sample_sheet3.csv"
ALIGN_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/3_Alignment"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/4_BQSR"
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
INDELS_GOLD="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
SNPS_1000G="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

# === Get sample info from sheet ===
LINE=$(tail -n +2 "$SHEET" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
SAMPLE=$(echo "$LINE" | cut -d',' -f1)
FASTQ_R1=$(echo "$LINE" | cut -d',' -f2)

INPUT_BAM="${ALIGN_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam"
WORKDIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$WORKDIR"

DEDUP_BAM="${WORKDIR}/${SAMPLE}.dedup.bam"
METRICS="${WORKDIR}/${SAMPLE}.dedup.metrics.txt"
BQSR_TABLE="${WORKDIR}/${SAMPLE}.recal.table"
RECAL_BAM="${WORKDIR}/${SAMPLE}.bqsr.bam"

# === Extract Read Group ===
echo "[INFO] Extracting Read Group for $SAMPLE"
if [[ "$FASTQ_R1" == *.gz ]]; then
    FIRST_READ=$(zcat "$FASTQ_R1" | head -1)
else
    FIRST_READ=$(head -1 "$FASTQ_R1")
fi

READ_NAME=${FIRST_READ#@}
IFS=':' read -ra PARTS <<< "$READ_NAME"
if [[ ${#PARTS[@]} -ge 4 ]]; then
    RG_ID="${PARTS[2]}.L${PARTS[3]}"
    RG_PU="${PARTS[2]}.L${PARTS[3]}"
else
    RG_ID="$SAMPLE"
    RG_PU="$SAMPLE"
fi

# === Step 1: AddOrReplaceReadGroups ===
TMP_RG_BAM="${WORKDIR}/${SAMPLE}.rg.bam"
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
echo "[INFO] Marking duplicates"
singularity exec $BIND $SIF \
    gatk MarkDuplicates \
    -I "$TMP_RG_BAM" \
    -O "$DEDUP_BAM" \
    -M "$METRICS" \
    --CREATE_INDEX true

# === Step 3: BaseRecalibrator ===
echo "[INFO] Running BQSR"
singularity exec $BIND $SIF \
    gatk BaseRecalibrator \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --known-sites "$INDELS_GOLD" \
    --known-sites "$SNPS_1000G" \
    -O "$BQSR_TABLE"

singularity exec $BIND $SIF \
    gatk ApplyBQSR \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --bqsr-recal-file "$BQSR_TABLE" \
    -O "$RECAL_BAM"

echo "[DONE] AddOrReplaceReadGroups MarkDuplicates BaseRecalibrator have completed for $SAMPLE"
