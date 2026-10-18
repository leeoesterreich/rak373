#!/bin/bash

#SBATCH --job-name=Alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=16  
#SBATCH --output=slurm_log/align/align_%A_%a.out
#SBATCH --error=slurm_log/align/align_%A_%a.err        
#SBATCH -t 1-00:00
#SBATCH --array=0-19                

module purge
module load gcc/8.2.0
module load star/2.7.11b

# Directories
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
mkdir -p "$designated_path/2_alignment_file"

# Reference genome
genome_dir="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7"
gtf_file="$genome_dir/gencode.v47.basic.annotation.gtf"

# Input/Output
input_dir="/ix1/alee/LO_LAB/General/Lab_Data/20250609_HSSC_RNAseq_pILC_Jagmohan"
output_dir="$designated_path/2_alignment_file"

threads=16

# Create a list of samples
read1_files=($input_dir/*_R1_001.fastq.gz)
read1=${read1_files[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "$read1" _R1_001.fastq.gz)
read2="$input_dir/${sample}_R2_001.fastq.gz"
output_prefix="$output_dir/${sample}_"
log_file="${output_dir}/${sample}.log"

# Check inputs
if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
    echo "FASTQ files for sample $sample not found. Exiting." | tee "$log_file"
    exit 1
fi

echo "Aligning $sample with STAR..." | tee "$log_file"
STAR --genomeDir "$genome_dir" \
     --readFilesIn "$read1" "$read2" \
     --readFilesCommand zcat \
     --runThreadN "$threads" \
     --twopassMode Basic \
     --outFileNamePrefix "$output_prefix" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM \
     --quantMode TranscriptomeSAM GeneCounts \
     > "$log_file" 2>&1

echo "Sample $sample alignment complete." | tee -a "$log_file"
