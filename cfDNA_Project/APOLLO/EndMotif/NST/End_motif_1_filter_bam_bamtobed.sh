#!/bin/bash

# This script (I) is used for preprocessing and generating bedpe file.
# Filtering criteria - properly paired (‘-f 2’), exclude unmapped reads, vendor quality check failures,
# PCR duplicates, no secondary alignments (‘-F 3844’), and have a mapping quality score (MAPQ) >= 30
# Sort BAM file by read name using ‘samtools sort’ command
# Convert BAM alignments to BEDPE format
# Calculate cfDNA fragment size based on genomic coordinates of paired-end reads

#SBATCH --job-name=End_motiff_analysis
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --output=slurm_log/Endmotif/Endmotif_%A_%a.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_%A_%a.err
#SBATCH -t 1-00:00
#SBATCH --array=0-261%50  # 262 BAMs, max 50 jobs at once

# Load required modules
module load samtools/1.22.1
module load bedtools/2.31.1

# Create designated output directory
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST"
output_dir="$designated_path/End_motif/bamtobed"
mkdir -p "$output_dir"

# Path to BAM list file
bam_list="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/bam_list.txt"

# Get the BAM file corresponding to this array task
task_id=${SLURM_ARRAY_TASK_ID}

# sed is 1-based, so add 1 to the 0-based array index
bam=$(sed -n "$((task_id+1))p" "$bam_list")

# Safety check
if [[ -z "$bam" ]]; then
    echo "No BAM file found for SLURM_ARRAY_TASK_ID=${task_id}"
    exit 1
fi

# Extract sample name (remove .bam)
filename=$(basename "$bam")
sample=${filename%.bam}

echo "SLURM_ARRAY_TASK_ID=${task_id}"
echo "BAM file: ${bam}"
echo "Sample: ${sample}"

# Apply filtering criteria
echo "Applying filtering for ${sample}..."
filter_bam="${output_dir}/${sample}_filter.bam"
samtools view -bh -f 2 -F 3844 -q 30 "$bam" > "$filter_bam"
echo "Filtering complete for ${sample}."

# Summary statistics for the filtered BAM file
stat_file="${output_dir}/${sample}_flagstat.txt"
echo "Generating stat report: $stat_file"
samtools flagstat "$filter_bam" > "$stat_file"

# Sort filtered BAM by read name
echo "Sorting ${sample} by name..."
sorted_bam_name="${output_dir}/${sample}_filter_sorted.bam"
samtools sort -n "$filter_bam" -o "$sorted_bam_name"
echo "Sorting complete for ${sample}."

# Remove intermediate filtered unsorted bam
rm -f "$filter_bam"

# Convert to BEDPE format
echo "Converting ${sample} to BEDPE format..."
bedpe_file="${output_dir}/${sample}.bedpe"
bedtools bamtobed -i "$sorted_bam_name" -bedpe > "$bedpe_file"
echo "BEDPE conversion complete for ${sample}."

# Calculate insert size (fragment length) and overwrite BEDPE file
echo "Modifying BEDPE file for ${sample}..."
modified_bedpe_file="${output_dir}/${sample}_modified.bedpe"
awk 'OFS = "\t" {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $3-$2}' \
    "$bedpe_file" > "$modified_bedpe_file"
mv "$modified_bedpe_file" "$bedpe_file"
echo "Modification complete for ${sample}."

echo "Done processing ${sample}."
