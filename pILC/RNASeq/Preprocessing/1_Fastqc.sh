#!/bin/bash

#This script is used to check the quality of high-throughput sequence data.

#SBATCH --job-name=FastQC
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 1-00:00

# Modules required
module purge
module load  fastqc/0.11.7
module load  multiqc/1.19

#Create new directory
designated_folder="/ix1/alee/LO_LAB/Personal/Rahul/pILC"
mkdir -p "$designated_folder/1_fastqc_result"

#Input directory
input_dir="/ix1/alee/LO_LAB/General/Lab_Data/20250609_HSSC_RNAseq_pILC_Jagmohan"

#Output directory 
output_dir="$designated_folder/1_fastqc_result"

#Loop through all FASTQ files in the input directory
echo "Analyzing fastq files for each sample..."

for fastq_file in $input_dir/*.fastq.gz; do
    # Extract the filename without extension
    sample=$(basename "$fastq_file" .fastq.gz)
    
    # Run FastQC
	fastqc -o "$output_dir" "$fastq_file" > "$output_dir/${sample}_fastqc.log"
    
    echo "FastQC analysis completed for $sample"
done

echo "All FastQC analyses completed"

#MultiQC is a reporting tool that aggregate results from bioinformatics analyses across many samples into a single report.

#Run MultiQC on the input directory
echo "Analyzing fastq files for all samples..."

multiqc "$output_dir" -o "$output_dir"

echo "MultiQC analysis completed"