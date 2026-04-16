#!/bin/bash

#SBATCH --job-name=Trim_fastqc
#SBATCH --output=slurm_log/Trim_fastqc/Trim_fastqc_TP16-M397.out
#SBATCH --error=slurm_log/Trim_fastqc/Trim_fastqc_TP16-M397.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load FastQC module
module purge
module load fastqc/0.11.7

# Set paths
sample_name="TP16-M397"
trimmed_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/2_TrimGalore/${sample_name}"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/2_TrimGalore/Fastqc"
mkdir -p "$output_dir" "slurm_log/Trim_fastqc"

# Trimmed FASTQ files
R1_trimmed="${trimmed_dir}/TP16-M397_B_OES757A8_CL_S5_R1_001_val_1.fq.gz"
R2_trimmed="${trimmed_dir}/TP16-M397_B_OES757A8_CL_S5_R2_001_val_2.fq.gz"

echo "Running FastQC for $sample_name"
echo "  R1_trimmed: $R1_trimmed"
echo "  R2_trimmed: $R2_trimmed"

# Run FastQC on R1
if [[ -f "$R1_trimmed" ]]; then
    fastqc -o "$output_dir" "$R1_trimmed" > "$output_dir/${sample_name}_R1_fastqc.log" 2>&1
    echo "FastQC done for R1"
else
    echo "R1 trimmed file not found: $R1_trimmed"
fi

# Run FastQC on R2
if [[ -f "$R2_trimmed" ]]; then
    fastqc -o "$output_dir" "$R2_trimmed" > "$output_dir/${sample_name}_R2_fastqc.log" 2>&1
    echo "FastQC done for R2"
else
    echo "R2 trimmed file not found: $R2_trimmed"
fi

echo "Finished FastQC for $sample_name"
