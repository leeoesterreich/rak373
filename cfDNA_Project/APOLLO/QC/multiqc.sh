#!/bin/bash

#SBATCH --job-name=MultiQC
#SBATCH --output=slurm_log/multiqc/multiqc_%j.out
#SBATCH --error=slurm_log/multiqc/multiqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# ==== Load module ====
module purge
module load multiqc/1.30

# ==== Directories ====
input_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Qualimap/"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Qualimap/"

# Ensure log dir exists
mkdir -p slurm_log/multiqc

echo "[INFO] Running MultiQC..."
multiqc "$input_dir" --outdir "$output_dir" --force

echo "[INFO] MultiQC complete. Report saved to: $output_dir/multiqc_report.html"
