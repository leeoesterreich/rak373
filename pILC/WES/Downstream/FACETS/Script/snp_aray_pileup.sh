#!/bin/bash

# This script is used to generate *.snp_pileup.gz by using run snp-pileup-wrapper.R (https://github.com/mskcc/facets-suite/tree/master)


#SBATCH --job-name=snp_pileup
#SBATCH --output=slurm_log/snp_pileup_%A_%a.out
#SBATCH --error=slurm_log/snp_pileup_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=0-15  

module load r/4.4.0

CSV_PATH="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Paired_samples.csv"
RSCRIPT_PATH="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Script/snp-pileup-wrapper.R"
SNP_PILEUP_PATH="/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4/facets/extcode/snp-pileup"
VCF_FILE="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/00-common_all.vcf.gz"
PILEUP_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/snp-pileup"

mkdir -p "$PILEUP_DIR"

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 2))  # Skip header
LINE=$(sed -n "${LINE_NUM}p" "$CSV_PATH")

if [ -z "$LINE" ]; then
    echo "[$(date)] No line found for task ID $SLURM_ARRAY_TASK_ID (line $LINE_NUM). Skipping."
    exit 1
fi

IFS=',' read -r SAMPLE TUMOR_BAM NORMAL_BAM TUMOR_BAM_PATH NORMAL_BAM_PATH <<< "$LINE"

OUTFILE="${PILEUP_DIR}/${SAMPLE}.snp_pileup.gz"
if [ -f "$OUTFILE" ]; then
    echo "[$(date)] Output already exists for $SAMPLE. Skipping."
    exit 0
fi

echo "[$(date)] Processing sample: $SAMPLE"

Rscript "$RSCRIPT_PATH" \
    --snp-pileup-path "$SNP_PILEUP_PATH" \
    --vcf-file "$VCF_FILE" \
    --normal-bam "$NORMAL_BAM_PATH" \
    --tumor-bam "$TUMOR_BAM_PATH" \
    --output-prefix "$PILEUP_DIR/$SAMPLE"

echo "[$(date)] Completed: $SAMPLE"
