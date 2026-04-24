#!/bin/bash

# This script is used to call CNV by CNVkit with normal-tumor paired samples.

#SBATCH --job-name=cnvkit_CNV
#SBATCH --output=slurm_log/cnvkit_%A_%a.out
#SBATCH --error=slurm_log/cnvkit_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-9  # Adjust this based on number of sample pairs in the CSV

# Load the CNVkit module (Singularity-based)
module purge
module load cnvkit/0.9.12

# Optional: suppress graphical backend and warnings
export MPLBACKEND=Agg
export PYTHONWARNINGS="ignore"

# ==== CONFIGURATION ====
samples_csv="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/CNVKit/Paired_samples.csv"
targets_bed="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/Exome_v1_hg38_Probes_Standard.filtered.bed"
refFlat="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/refFlat.filtered.txt"
fasta="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/Human_genome_hg38_p7_chr_only.fna"
access_bed="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/CNVKit/access-5kb-mappable.hg38.clean.bed"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/CNVKit/Result/Tumor_Normal"

# ==== SETUP ====
mkdir -p "$output_dir"

# Get the current line from the sample CSV (skip header, task 0 gets line 2)
sample_line=$(tail -n +2 "$samples_csv" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Read columns from the CSV
IFS=',' read -r SAMPLE TUMOR_BAM NORMAL_BAM TUMOR_BAM_PATH NORMAL_BAM_PATH <<< "$sample_line"

# Validate input BAMs
if [[ ! -f "$NORMAL_BAM_PATH" ]]; then
    echo "ERROR: Normal BAM not found: $NORMAL_BAM_PATH"
    exit 1
fi

if [[ ! -f "$TUMOR_BAM_PATH" ]]; then
    echo "ERROR: Tumor BAM not found: $TUMOR_BAM_PATH"
    exit 1
fi

# Create per-sample output directory
SAMPLE_OUT="${output_dir}/${SAMPLE}"
mkdir -p "$SAMPLE_OUT"

# ==== RUN CNVKIT ====
echo "[$(date)] Starting CNVkit for sample: $SAMPLE"

cnvkit.py batch "$TUMOR_BAM_PATH" \
    --normal "$NORMAL_BAM_PATH" \
    --targets "$targets_bed" \
    --annotate "$refFlat" \
    --fasta "$fasta" \
    --access "$access_bed" \
    --output-reference "$SAMPLE_OUT/my_reference.cnn" \
    --output-dir "$SAMPLE_OUT" \
    --diagram --scatter

echo "[$(date)] Completed CNVkit for sample: $SAMPLE"