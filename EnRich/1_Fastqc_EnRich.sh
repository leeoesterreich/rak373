#!/bin/bash

# This script is used to check the quality of high-throughput sequence data.

#SBATCH --job-name=FastQC
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 1-00:00:00
#SBATCH --array=0-16
#SBATCH --output=slurm_logs/FastQC_%A_%a.out
#SBATCH --error=slurm_logs/FastQC_%A_%a.err

module purge
module load fastqc/0.11.7
module load multiqc/1.30

designated_folder="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
output_dir="$designated_folder/1_fastqc_result_010526"
mkdir -p "$output_dir" slurm_logs

input_dir="/ix1/alee/LO_LAB/General/Lab_Data/20251224_SeqCenter_LBP_Enrich_Rahul/EnRICH_pilot"

# Collect only TP samples (exclude Undetermined)
FASTQ_FILES=($(ls $input_dir/TP*.fastq.gz))

fastq_file="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"
sample=$(basename "$fastq_file" .fastq.gz)

fastqc -t 8 -o "$output_dir" "$fastq_file" \
  > "$output_dir/${sample}_fastqc.log"

echo "FastQC completed for $sample"
