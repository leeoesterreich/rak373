#!/bin/bash

# This script is used to call somatic variants in a tumor only mode

#SBATCH --job-name=mutect2_tumoronly
#SBATCH --output=slurm_log/mutect2/tumor_only/mutect2_tumoronly_%A_%a.out
#SBATCH --error=slurm_log/mutect2/tumor_only/mutect2_tumoronly_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-7

module load singularity/4.3.2

log_info()    { echo "[INFO]    $*"; }
log_fail()    { echo "[ERROR]   $*" >&2; exit 1; }
log_success() { echo "[SUCCESS] $*"; }

SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet.csv"
BQSR_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/3_BQSR_010526"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/4_Mutect2_tumoronly_010526"

SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"
INTERVALS="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Bed_file/combined_good_CDH1.sorted.interval_list"
GERMLINE_RESOURCE="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/gnomad.exomes.r2.1.1.sites.hg19.ucsc.vcf.gz"

mkdir -p "$OUTPUT_DIR" slurm_log/mutect2/tumor_only

threads="${SLURM_CPUS_PER_TASK:-1}"

# ---- Get sample from CSV (Sample,BAM_Path) ----
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$SHEET" | tr -d '\r')
[[ -n "$LINE" ]] || { log_info "No line for task ${SLURM_ARRAY_TASK_ID}; exiting."; exit 0; }

SAMPLE=$(echo "$LINE" | cut -d',' -f1 | xargs)
[[ -n "$SAMPLE" ]] || log_fail "Failed to parse SAMPLE from: $LINE"

TUMOR_BAM="${BQSR_DIR}/${SAMPLE}/${SAMPLE}.bqsr.bam"
[[ -f "$TUMOR_BAM" ]] || log_fail "Tumor BAM not found: $TUMOR_BAM"
[[ -f "${TUMOR_BAM}.bai" ]] || log_fail "Tumor BAM index (.bai) not found: ${TUMOR_BAM}.bai"

# ---- Outputs ----
OUT_DIR="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

UNFILTERED_VCF="${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz"
F1R2_TAR="${OUT_DIR}/${SAMPLE}.f1r2.tar.gz"
ORIENT_PRIORS="${OUT_DIR}/${SAMPLE}.orientation-priors.tar.gz"
PILEUP="${OUT_DIR}/${SAMPLE}.pileups.table"
CONTAM_TABLE="${OUT_DIR}/${SAMPLE}.contamination.table"
SEGMENTS_TABLE="${OUT_DIR}/${SAMPLE}.segments.table"
FILTERED_VCF="${OUT_DIR}/${SAMPLE}.filtered.vcf.gz"

[[ -f "$FILTERED_VCF" ]] && { log_success "$SAMPLE already processed — skipping."; exit 0; }

# STEP 1: Mutect2 (tumor-only)
log_info "Running Mutect2 (tumor-only) for $SAMPLE"
singularity exec $BIND $SIF gatk Mutect2 \
  -R "$REF" \
  -I "$TUMOR_BAM" \
  -L "$INTERVALS" \
  --germline-resource "$GERMLINE_RESOURCE" \
  --f1r2-tar-gz "$F1R2_TAR" \
  --native-pair-hmm-threads "$threads" \
  -max-mnp-distance 0 \
  -O "$UNFILTERED_VCF" || log_fail "Mutect2 failed"

# STEP 2: Learn Read Orientation Model
log_info "Learning orientation model"
singularity exec $BIND $SIF gatk LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIENT_PRIORS" || log_fail "LearnReadOrientationModel failed"

# STEP 3: GetPileupSummaries
log_info "Running GetPileupSummaries"
singularity exec $BIND $SIF gatk GetPileupSummaries \
  -I "$TUMOR_BAM" \
  -L "$INTERVALS" \
  -V "$GERMLINE_RESOURCE" \
  -R "$REF" \
  -O "$PILEUP" || log_fail "GetPileupSummaries failed"

# STEP 4: CalculateContamination
log_info "Running CalculateContamination"
singularity exec $BIND $SIF gatk CalculateContamination \
  -I "$PILEUP" \
  --tumor-segmentation "$SEGMENTS_TABLE" \
  -O "$CONTAM_TABLE" || log_fail "CalculateContamination failed"

# STEP 5: FilterMutectCalls
log_info "Filtering Mutect2 calls"
singularity exec $BIND $SIF gatk FilterMutectCalls \
  -R "$REF" \
  -V "$UNFILTERED_VCF" \
  --contamination-table "$CONTAM_TABLE" \
  --tumor-segmentation "$SEGMENTS_TABLE" \
  --orientation-bias-artifact-priors "$ORIENT_PRIORS" \
  -O "$FILTERED_VCF" || log_fail "FilterMutectCalls failed"

log_success "$SAMPLE complete — final VCF: $FILTERED_VCF"
