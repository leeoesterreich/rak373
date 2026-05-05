#!/bin/bash
#SBATCH --job-name=mutect2_pon
#SBATCH --output=slurm_log/mutect2/pon/mutect2_pon_%A_%a.out
#SBATCH --error=slurm_log/mutect2/pon/mutect2_pon_%A_%a.err
#SBATCH --time=1-00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-16

set -euo pipefail

module purge
module load singularity/4.3.2

# Constants
BQSR_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4_BQSR"
OUT_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/5.2_Mutect2_PON_normals"
REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

mkdir -p "$OUT_ROOT" "slurm_log/mutect2/pon"

# Get normal sample directories directly from BQSR output
mapfile -t NORMAL_DIRS < <(find "$BQSR_ROOT" -mindepth 1 -maxdepth 1 -type d -name "*_N" | sort)

SAMPLE_DIR="${NORMAL_DIRS[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "$SAMPLE_DIR" || ! -d "$SAMPLE_DIR" ]]; then
    log_error "No valid normal sample directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

SAMPLE=$(basename "$SAMPLE_DIR")

# Input and output paths
RECAL_BAM="${SAMPLE_DIR}/${SAMPLE}.bqsr.bam"
RECAL_BAI="${SAMPLE_DIR}/${SAMPLE}.bqsr.bai"
OUT_DIR="${OUT_ROOT}/${SAMPLE}"
OUT_VCF="${OUT_DIR}/${SAMPLE}.vcf.gz"

mkdir -pv "$OUT_DIR"

# Checks
if [[ ! -f "$RECAL_BAM" ]]; then
    log_error "Missing input BAM: $RECAL_BAM"
    exit 1
fi

if [[ ! -f "$RECAL_BAI" ]]; then
    log_error "Missing GATK-style BAM index: $RECAL_BAI"
    exit 1
fi

# Skip if output already exists
if [[ -f "$OUT_VCF" ]]; then
    log_success "$SAMPLE VCF already exists — skipping."
    exit 0
fi

# Run Mutect2 with max-mnp-distance 0 for PoN workflow
singularity exec $BIND $SIF \
    gatk Mutect2 \
    -R "$REF" \
    -I "$RECAL_BAM" \
    --max-mnp-distance 0 \
    -O "$OUT_VCF"

log_success "$SAMPLE Mutect2 (PoN mode) complete."