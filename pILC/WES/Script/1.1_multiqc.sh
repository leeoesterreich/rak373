#!/bin/bash

#SBATCH --job-name=MultiQC
#SBATCH --output=slurm_log/multiqc/multiqc_%j.out
#SBATCH --error=slurm_log/multiqc/multiqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Load module
module purge
module load multiqc/1.30

# Directories
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/1_fastqc_result/"

mkdir -p slurm_log/multiqc
mkdir -p "$output_dir"

echo "Running MultiQC on: $output_dir"

multiqc "$output_dir" -o "$output_dir" -n multiqc_fastqc_report.html

echo "MultiQC complete."