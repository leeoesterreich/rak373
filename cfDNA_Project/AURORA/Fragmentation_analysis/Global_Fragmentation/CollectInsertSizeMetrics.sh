#!/bin/bash

#This script is used to derived insert size distribution and read orientation of paired-end libraries.

#SBATCH --job-name=InsertSizeMetrics
#SBATCH --output=logs/InsertSizeMetrics_%A_%a.out
#SBATCH --error=logs/InsertSizeMetrics_%A_%a.err
#SBATCH --array=0-383%50  # MaxArraySize: 500
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 1-00:00:00

# Load module
module purge
module load picard/3.0.0
module load r/4.4.0

# Paths
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Global_Fragment_analysis/Result"
mkdir -p "$output_dir"

# Read bam_samples.txt into array
mapfile -t lines < bam_samples_2.txt
line="${lines[$SLURM_ARRAY_TASK_ID]}"

# Parse sample and BAM path
sample_name=$(echo "$line" | awk '{print $1}')
bam_path=$(echo "$line" | awk '{print $2}')

# Output files
metrics="${output_dir}/${sample_name}.metrics.txt"
histogram="${output_dir}/${sample_name}.histogram.pdf"

echo "[$sample_name] Processing BAM: $bam_path"

# Run Picard CollectInsertSizeMetrics
picard CollectInsertSizeMetrics \
  -I "$bam_path" \
  -O "$metrics" \
  -H "$histogram" \
  -M 0.05

echo "[$sample_name] Done."
