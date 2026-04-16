#!/bin/bash

# This script is used to call somatic variants in a tumor with matched normal mode

#SBATCH --job-name=mutect2_paired
#SBATCH --output=slurm_log/mutect2/paired_normal/mutect2_paired_%A_%a.out
#SBATCH --error=slurm_log/mutect2/paired_normal/mutect2_paired_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-10

module load singularity/4.3.2

# Logging functions
log_info()    { echo "[INFO]    $*"; }
log_fail()    { echo "[ERROR]   $*" >&2; exit 1; }
log_success() { echo "[SUCCESS] $*"; }

# INPUTS 
PAIRS="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/5_PON/Normal_Tumor.csv"
BQSR_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/4_BQSR"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/6_Mutect2"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

# GATK References
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/af-only-gnomad.hg38.vcf.gz"
PON="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/5_PON/pon_rak373_20251023_161850.vcf.gz"

# GET SAMPLE NAMES 
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PAIRS")
NORMAL=$(echo "$LINE" | cut -d',' -f1)
TUMOR=$(echo "$LINE" | cut -d',' -f2)

# FIND BAM FILES 
TUMOR_BAM=$(find "$BQSR_DIR" -type f -name "${TUMOR}.bqsr.bam")
NORMAL_BAM=$(find "$BQSR_DIR" -type f -name "${NORMAL}.bqsr.bam")

[[ -f "$TUMOR_BAM" ]] || log_fail "Tumor BAM not found: $TUMOR"
[[ -f "$NORMAL_BAM" ]] || log_fail "Normal BAM not found: $NORMAL"

# OUTPUT FILES 
OUT_DIR="${OUTPUT_DIR}/${TUMOR}_vs_${NORMAL}"
mkdir -p "$OUT_DIR"

UNFILTERED_VCF="${OUT_DIR}/${TUMOR}.unfiltered.vcf.gz"
F1R2_TAR="${OUT_DIR}/${TUMOR}.f1r2.tar.gz"
ORIENT_PRIORS="${OUT_DIR}/${TUMOR}.orientation-priors.tar.gz"
PILEUP="${OUT_DIR}/${TUMOR}.pileups.table"
CONTAM_TABLE="${OUT_DIR}/${TUMOR}.contamination.table"
SEGMENTS_TABLE="${OUT_DIR}/${TUMOR}.segments.table"
FILTERED_VCF="${OUT_DIR}/${TUMOR}.filtered.vcf.gz"

# SKIP IF DONE 
[[ -f "$FILTERED_VCF" ]] && { log_success "$TUMOR already processed — skipping."; exit 0; }

# STEP 1: Mutect2 
log_info "Running Mutect2 on $TUMOR vs $NORMAL"
singularity exec $BIND $SIF gatk Mutect2 \
  -R "$REF" \
  -I "$TUMOR_BAM" \
  -I "$NORMAL_BAM" \
  -tumor-sample "$TUMOR" \
  -normal "$NORMAL" \
  --germline-resource "$GERMLINE_RESOURCE" \
  --panel-of-normals "$PON" \
  --f1r2-tar-gz "$F1R2_TAR" \
  -O "$UNFILTERED_VCF" || log_fail "Mutect2 failed"

# STEP 2: Learn Read Orientation Model 
log_info "Learning orientation model"
singularity exec $BIND $SIF gatk LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIENT_PRIORS" || log_fail "Orientation model failed"

# STEP 3: GetPileupSummaries 
log_info "Running GetPileupSummaries"
singularity exec $BIND $SIF gatk GetPileupSummaries \
  -I "$TUMOR_BAM" \
  -V "$GERMLINE_RESOURCE" \
  -R "$REF" \
  -O "$PILEUP" || log_fail "GetPileupSummaries failed"

# STEP 4: CalculateContamination 
log_info "Running CalculateContamination"
singularity exec $BIND $SIF gatk CalculateContamination \
  -I "$PILEUP" \
  -tumor-segmentation "$SEGMENTS_TABLE" \
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

log_success "$TUMOR vs $NORMAL complete — final VCF: $FILTERED_VCF"
