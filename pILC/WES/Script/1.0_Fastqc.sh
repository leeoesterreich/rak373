#!/bin/bash

#SBATCH --job-name=FastQC
#SBATCH --output=slurm_log/fastqc/fastqc_%A_%a.out
#SBATCH --error=slurm_log/fastqc/fastqc_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-47

# Load modules
module purge
module load fastqc/0.11.7

# Base input directory containing 48 sample folders
base_dir="/ix1/alee/LO_LAB/General/Lab_Data/20250611_HTGC_WES_pILC_Jagmohan"

# Output directory
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/1_fastqc_result"
mkdir -p "$output_dir"
mkdir -p slurm_log/fastqc

# Get list of sample folders
sample_dirs=($(find "$base_dir" -mindepth 1 -maxdepth 1 -type d | sort))

# Pick the folder for this array task
sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]}"

# Safety check
if [[ -z "$sample_dir" || ! -d "$sample_dir" ]]; then
    echo "Error: No valid sample directory found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

sample_name=$(basename "$sample_dir")

echo "Processing sample folder: $sample_name"
echo "Folder path: $sample_dir"

# Find all FASTQ files in this sample folder
fastq_files=($(find "$sample_dir" -maxdepth 1 -type f -name "*.fastq.gz" | sort))

# Check if FASTQ files exist
if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "No FASTQ files found in $sample_dir"
    exit 1
fi

echo "FASTQ files found:"
printf '%s\n' "${fastq_files[@]}"

# Optional: make sample-specific output folder
sample_outdir="$output_dir/$sample_name"
mkdir -p "$sample_outdir"

# Run FastQC on all FASTQ files in this folder
fastqc -t 8 -o "$sample_outdir" "${fastq_files[@]}" \
    > "$sample_outdir/${sample_name}_fastqc.log" 2>&1

echo "Finished FastQC for $sample_name"