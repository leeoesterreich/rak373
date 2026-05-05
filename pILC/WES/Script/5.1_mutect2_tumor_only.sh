#!/bin/bash
#SBATCH --job-name=mutect2_tumoronly
#SBATCH --output=slurm_log/mutect2/tumor_only/mutect2_tumoronly_%A_%a.out
#SBATCH --error=slurm_log/mutect2/tumor_only/mutect2_tumoronly_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-30

set -euo pipefail

module purge
module load singularity/4.3.2

# ===============================
# Logging
# ===============================
log_info()    { echo -e "\033[1;34m[$(date '+%F %T')] [INFO] $1\033[0m"; }
log_fail()    { echo -e "\033[1;31m[$(date '+%F %T')] [ERROR] $1\033[0m" >&2; exit 1; }
log_success() { echo -e "\033[1;32m[$(date '+%F %T')] [DONE] $1\033[0m"; }

# ===============================
# Paths
# ===============================
BQSR_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4_BQSR"
OUT_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/5.1_mutect2_tumor_only"
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
GERMLINE_RESOURCE="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/af-only-gnomad.hg38.vcf.gz"
REGIONS="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"

SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUT_ROOT" "slurm_log/mutect2/tumor_only"

# ===============================
# Get tumor sample folders
# ===============================
mapfile -t TUMOR_DIRS < <(
    find "$BQSR_ROOT" -mindepth 1 -maxdepth 1 -type d -name "*_T" | sort
)

NUM_TUMORS=${#TUMOR_DIRS[@]}
log_info "Total tumor samples found: $NUM_TUMORS"

if [[ "$NUM_TUMORS" -eq 0 ]]; then
    log_fail "No tumor sample folders found under $BQSR_ROOT"
fi

if [[ "$SLURM_ARRAY_TASK_ID" -ge "$NUM_TUMORS" ]]; then
    log_info "Task ID ${SLURM_ARRAY_TASK_ID} exceeds tumor count ${NUM_TUMORS}; exiting."
    exit 0
fi

SAMPLE_DIR="${TUMOR_DIRS[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "$SAMPLE_DIR")

# ===============================
# Input/output files
# ===============================
RECAL_BAM="${SAMPLE_DIR}/${SAMPLE}.bqsr.bam"
RECAL_INDEX="${RECAL_BAM%.bam}.bai"  # sample.bqsr.bai

OUT_DIR="${OUT_ROOT}/${SAMPLE}"
mkdir -p "$OUT_DIR"

UNFILTERED_VCF="${OUT_DIR}/${SAMPLE}.tumor_only.unfiltered.vcf.gz"
FILTERED_VCF="${OUT_DIR}/${SAMPLE}.tumor_only.filtered.vcf.gz"

# ===============================
# Checks
# ===============================
[[ ! -f "$RECAL_BAM" ]] && log_fail "Missing input BAM: $RECAL_BAM"

if [[ -f "$RECAL_INDEX" ]]; then
    log_info "Found BAM index: $RECAL_INDEX"
else
    log_fail "Missing BAM index: $RECAL_INDEX"
fi

# Skip if output already exists
if [[ -f "$FILTERED_VCF" ]]; then
    log_success "$SAMPLE tumor-only filtered VCF already exists — skipping."
    exit 0
fi

# ===============================
# Step 1: Run Mutect2 tumor-only
# ===============================
log_info "Running tumor-only Mutect2 for $SAMPLE"

singularity exec $BIND "$SIF" \
    gatk Mutect2 \
    -R "$REF" \
    -I "$RECAL_BAM" \
    -tumor "$SAMPLE" \
    --germline-resource "$GERMLINE_RESOURCE" \
    --af-of-alleles-not-in-resource 0.0000025 \
    -L "$REGIONS" \
    --max-mnp-distance 0 \
    --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
    -O "$UNFILTERED_VCF"

# ===============================
# Step 2: FilterMutectCalls
# ===============================
log_info "Filtering tumor-only Mutect2 calls for $SAMPLE"

singularity exec $BIND "$SIF" \
    gatk FilterMutectCalls \
    -R "$REF" \
    -V "$UNFILTERED_VCF" \
    -O "$FILTERED_VCF"

log_success "$SAMPLE tumor-only Mutect2 complete: $FILTERED_VCF"