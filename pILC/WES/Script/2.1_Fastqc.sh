#!/bin/bash

#SBATCH --job-name=Trim_fastqc
#SBATCH --output=slurm_log/Trim_fastqc/Trim_fastqc_%A_%a.out
#SBATCH --error=slurm_log/Trim_fastqc/Trim_fastqc_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-47

set -euo pipefail

# Load FastQC module
module purge
module load fastqc/0.11.7

# Base directories
trim_base="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/2_TrimGalore"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/2_TrimGalore/Fastqc"

mkdir -p "$output_dir"
mkdir -p slurm_log/Trim_fastqc

# Get trimmed sample folders
mapfile -t sample_dirs < <(find "$trim_base" -mindepth 1 -maxdepth 1 -type d ! -name "Fastqc" | sort)

sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "${sample_dir}" || ! -d "${sample_dir}" ]]; then
    echo "[ERROR] No valid sample directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

sample_name=$(basename "$sample_dir")
sample_outdir="${output_dir}/${sample_name}"
mkdir -p "$sample_outdir"

# Detect trimmed FASTQ files
R1_trimmed=$(find "$sample_dir" -maxdepth 1 -type f -name "*_val_1.fq.gz" | sort | head -n 1)
R2_trimmed=$(find "$sample_dir" -maxdepth 1 -type f -name "*_val_2.fq.gz" | sort | head -n 1)

echo "Running FastQC for $sample_name"
echo "  R1_trimmed: ${R1_trimmed:-NOT FOUND}"
echo "  R2_trimmed: ${R2_trimmed:-NOT FOUND}"

# Check files exist
if [[ -z "${R1_trimmed}" || ! -f "${R1_trimmed}" ]]; then
    echo "[ERROR] R1 trimmed file not found for $sample_name"
    exit 1
fi

if [[ -z "${R2_trimmed}" || ! -f "${R2_trimmed}" ]]; then
    echo "[ERROR] R2 trimmed file not found for $sample_name"
    exit 1
fi

# Run FastQC on both trimmed files
fastqc -t 4 -o "$sample_outdir" "$R1_trimmed" "$R2_trimmed" \
    > "$sample_outdir/${sample_name}_trimmed_fastqc.log" 2>&1

echo "Finished FastQC for $sample_name"