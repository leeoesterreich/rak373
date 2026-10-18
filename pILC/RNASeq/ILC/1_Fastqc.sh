#!/bin/bash

#This script is used to check the quality of high-throughput sequence data.

#SBATCH --job-name=FastQC
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH -t 1-00:00

# Modules required
module purge
module load  fastqc/0.11.7
module load multiqc/1.30

# Set paths
designated_folder="/ix1/alee/LO_LAB/Personal/Rahul/pILC/Chelsea_21ILC_RNASeq_data"
input_dir="$designated_folder/metadata"
output_dir="$designated_folder/1_fastqc_result"

# Create output directory
mkdir -p "$output_dir"

echo "Analyzing fastq files for each sample..."

for fastq_file in "$input_dir"/*.fastq.gz; do
    sample=$(basename "$fastq_file" .fastq.gz)
    
    # Run FastQC using 16 threads
    fastqc -t 16 -o "$output_dir" "$fastq_file" > "$output_dir/${sample}_fastqc.log"
    
    echo "FastQC analysis completed for $sample"
done

echo "All FastQC analyses completed"

# Run MultiQC to summarize FastQC outputs
echo "Running MultiQC summary..."
multiqc "$output_dir" -o "$output_dir"

echo "MultiQC analysis completed"