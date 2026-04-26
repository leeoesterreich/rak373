#!/bin/bash
#SBATCH --job-name=mutect2_pon_full
#SBATCH --output=slurm_log/mutect2/with_pon/mutect2_pon_full_%A_%a.out
#SBATCH --error=slurm_log/mutect2/with_pon/mutect2_pon_full_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-30

set -euo pipefail

module purge
module load singularity/4.3.2

# Logging
log_info()    { echo -e "\033[1;34m[$(date '+%F %T')] [INFO] $1\033[0m"; }
log_fail()    { echo -e "\033[1;31m[$(date '+%F %T')] [ERROR] $1\033[0m"; exit 1; }
log_success() { echo -e "\033[1;32m[$(date '+%F %T')] [DONE] $1\033[0m"; }

# Paths
BQSR_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4_BQSR"
OUT_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/Mutect2_with_PON"

# Update this to your actual final PoN filename
PON="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/pon_rak373_20260422_152132.vcf.gz"

REGIONS="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"
COMMON_SITES="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/af-only-gnomad.hg38.vcf.gz"
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"

SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUT_ROOT" "slurm_log/mutect2/with_pon"

# Get tumor sample folders directly from BQSR output
mapfile -t TUMOR_DIRS < <(find "$BQSR_ROOT" -mindepth 1 -maxdepth 1 -type d -name "*_T" | sort)

SAMPLE_DIR="${TUMOR_DIRS[$SLURM_ARRAY_TASK_ID]:-}"

[[ -z "$SAMPLE_DIR" || ! -d "$SAMPLE_DIR" ]] && \
    log_fail "No valid tumor sample directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

SAMPLE=$(basename "$SAMPLE_DIR")

RECAL_BAM="${SAMPLE_DIR}/${SAMPLE}.bqsr.bam"
RECAL_BAI="${SAMPLE_DIR}/${SAMPLE}.bqsr.bai"

OUT_DIR="${OUT_ROOT}/${SAMPLE}"
mkdir -p "$OUT_DIR"

UNFILTERED_VCF="${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz"
UNFILTERED_STATS="${OUT_DIR}/${SAMPLE}.unfiltered.vcf.gz.stats"
F1R2_TAR="${OUT_DIR}/${SAMPLE}.f1r2.tar.gz"
ORIENT_PRIORS="${OUT_DIR}/${SAMPLE}.orientation-priors.tar.gz"
PILEUP="${OUT_DIR}/${SAMPLE}.pileups.table"
CONTAM_TABLE="${OUT_DIR}/${SAMPLE}.contamination.table"
SEGMENTS_TABLE="${OUT_DIR}/${SAMPLE}.segments.table"
FILTERED_VCF="${OUT_DIR}/${SAMPLE}.vcf.gz"

# Checks
[[ ! -f "$RECAL_BAM" ]] && log_fail "Missing input BAM: $RECAL_BAM"
[[ ! -f "$RECAL_BAI" ]] && log_fail "Missing GATK-style BAM index: $RECAL_BAI"
[[ ! -f "$PON" ]] && log_fail "Missing PoN file: $PON"
[[ ! -f "$COMMON_SITES" ]] && log_fail "Missing germline resource: $COMMON_SITES"

# Skip if already done
[[ -f "$FILTERED_VCF" ]] && \
    { log_success "$SAMPLE already processed — skipping."; exit 0; }

# Step 1: Run Mutect2
log_info "Running Mutect2 on $SAMPLE"
singularity exec $BIND $SIF gatk Mutect2 \
  -R "$REF" \
  -I "$RECAL_BAM" \
  --panel-of-normals "$PON" \
  --germline-resource "$COMMON_SITES" \
  -L "$REGIONS" \
  --max-mnp-distance 0 \
  --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
  --f1r2-tar-gz "$F1R2_TAR" \
  -O "$UNFILTERED_VCF" || log_fail "Mutect2 failed"

# Step 2: Learn Read Orientation Model
log_info "Learning orientation model"
singularity exec $BIND $SIF gatk LearnReadOrientationModel \
  -I "$F1R2_TAR" \
  -O "$ORIENT_PRIORS" || log_fail "Orientation model failed"

# Step 3: GetPileupSummaries
log_info "Running GetPileupSummaries"
singularity exec $BIND $SIF gatk GetPileupSummaries \
  -I "$RECAL_BAM" \
  -V "$COMMON_SITES" \
  -L "$REGIONS" \
  -R "$REF" \
  -O "$PILEUP" || log_fail "GetPileupSummaries failed"

# Step 4: CalculateContamination
log_info "Running CalculateContamination"
singularity exec $BIND $SIF gatk CalculateContamination \
  -I "$PILEUP" \
  --tumor-segmentation "$SEGMENTS_TABLE" \
  -O "$CONTAM_TABLE" || log_fail "CalculateContamination failed"

# Step 5: FilterMutectCalls
log_info "Filtering calls"
singularity exec $BIND $SIF gatk FilterMutectCalls \
  -R "$REF" \
  -V "$UNFILTERED_VCF" \
  --stats "$UNFILTERED_STATS" \
  --contamination-table "$CONTAM_TABLE" \
  --tumor-segmentation "$SEGMENTS_TABLE" \
  --orientation-bias-artifact-priors "$ORIENT_PRIORS" \
  -L "$REGIONS" \
  --f-score-beta 1.0 \
  -O "$FILTERED_VCF" || log_fail "FilterMutectCalls failed"

log_success "$SAMPLE complete — Filtered VCF written to $FILTERED_VCF"
