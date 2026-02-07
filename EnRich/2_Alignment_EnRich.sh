#!/bin/bash

# This script is used to align the samples with the reference genome and check the quality of BAM files.

# BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome
# BWA-MEM, for high-quality queries and better performance, as it is faster and more accurate

#SBATCH --job-name=Alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH -t 1-00:00:00
#SBATCH --array=0-7
#SBATCH --output=slurm_logs/alignment_%A_%a.out
#SBATCH --error=slurm_logs/alignment_%A_%a.err

# === Load Modules ===
module purge
module load bwa/0.7.17
module load samtools/1.22.1

# === Set Up Paths ===
designated_path1="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
output_dir="$designated_path1/2_alignment_010526"
input_dir="/ix1/alee/LO_LAB/General/Lab_Data/20251224_SeqCenter_LBP_Enrich_Rahul/EnRICH_pilot"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"

# === Create Directories ===
mkdir -p "$output_dir"
mkdir -p slurm_logs

# === Get Sample for This Array Task ===
samples=($(ls "$input_dir"/TP*_R1_001.fastq.gz | xargs -n 1 basename | sed 's/_R1_001.fastq.gz//' | sort))
sample="${samples[$SLURM_ARRAY_TASK_ID]}"

read1="$input_dir/${sample}_R1_001.fastq.gz"
read2="$input_dir/${sample}_R2_001.fastq.gz"

bam_file="$output_dir/${sample}.bam"
sorted_bam="$output_dir/${sample}_sorted.bam"
bai_file="$output_dir/${sample}_sorted.bam.bai"

echo "Processing sample: $sample"
echo "Input files:"
echo "  R1: $read1"
echo "  R2: $read2"

# === Input Check ===
if [[ ! -f "$read1" || ! -f "$read2" ]]; then
    echo "FASTQ files for $sample not found. Exiting."
    exit 1
fi

# === Alignment + BAM Generation ===
echo "Running BWA MEM and piping to samtools view..."
bwa mem -t 16 "$reference" "$read1" "$read2" | samtools view -@ 16 -bS - > "$bam_file"

# === Sorting ===
echo "Sorting BAM..."
samtools sort -@ 16 "$bam_file" -o "$sorted_bam"

# === Indexing ===
echo "Indexing sorted BAM..."
samtools index "$sorted_bam" "$bai_file"

# === Cleanup ===
echo "Cleaning up intermediate BAM..."
rm -f "$bam_file"

echo "Finished processing $sample"
