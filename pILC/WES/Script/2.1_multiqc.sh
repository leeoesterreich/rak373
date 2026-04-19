#!/bin/bash

#SBATCH --job-name=MultiQC
#SBATCH --output=slurm_log/multiqc/multiqc_%j.out
#SBATCH --error=slurm_log/multiqc/multiqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load module
module purge
module load multiqc/1.30

# Directories
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/2_TrimGalore/Fastqc"

echo "Running MultiQC..."
multiqc "$output_dir" -o "$output_dir"
echo "MultiQC complete."
