#!/bin/bash

# This script runs FastQC for each sample

#SBATCH --job-name=FastQC
#SBATCH --output=slurm_log/fastqc/fastqc_%A_%a.out
#SBATCH --error=slurm_log/fastqc/fastqc_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-19

# Load modules
module purge
module load fastqc/0.11.7

# Paths
sample_sheet="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/Sample_sheet.csv"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/1_fastqc_result"
mkdir -p "$output_dir"

# Read the corresponding line (skip header using tail +2)
line=$(tail -n +2 "$sample_sheet" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Extract values
IFS=',' read -r sample_name R1_path R2_path <<< "$line"

echo "Processing sample: $sample_name"
echo "  R1: $R1_path"
echo "  R2: $R2_path"

# Run FastQC if R1 exists
if [[ -n "$R1_path" && -f "$R1_path" ]]; then
    fastqc -o "$output_dir" "$R1_path" > "$output_dir/${sample_name}_R1_fastqc.log" 2>&1
    echo "FastQC done for R1"
else
    echo "Skipping R1 for $sample_name: file not found or empty"
fi

# Run FastQC if R2 exists
if [[ -n "$R2_path" && -f "$R2_path" ]]; then
    fastqc -o "$output_dir" "$R2_path" > "$output_dir/${sample_name}_R2_fastqc.log" 2>&1
    echo "FastQC done for R2"
else
    echo "Skipping R2 for $sample_name: file not found or empty"
fi

echo "Finished FastQC for $sample_name"
