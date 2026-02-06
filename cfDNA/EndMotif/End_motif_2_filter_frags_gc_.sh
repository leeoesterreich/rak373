#!/bin/bash
# This is a continuation script (II).
# Convert BEDPE file into BED format
# Sort BED file by chromosome and start position
# Filter out blacklisted and gap regions of hg38
# Calculate fragment GC content

#SBATCH --job-name=End_motif_analysis_step2
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --output=slurm_log/Endmotif/Endmotif_step2_%A_%a.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_step2_%A_%a.err
#SBATCH -t 1-00:00
#SBATCH --array=0-261%50  # 262 BEDPEs, max 50 jobs at once

# Load required modules 
module purge
module load bedtools/2.31.1

# Reference genome & filter files
reference="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/Human_genome_hg38_p7_chr_only.fna"
blacklist="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Blacklist_file/hg38-blacklist.v2.bed"
gap="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Blacklist_file/gaps.hg38.bed"

# Directory containing BEDPE files from step I
bedpe_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif/bamtobed"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif/filter_frags_gc"
mkdir -p "$output_dir"

# Collect all BEDPE files into an array
bedpe_files=( "$bedpe_dir"/*.bedpe )

# Safety: check that we actually found files
if [[ ${#bedpe_files[@]} -eq 0 || ! -e "${bedpe_files[0]}" ]]; then
    echo "No .bedpe files found in $bedpe_dir"
    exit 1
fi

echo "Found ${#bedpe_files[@]} .bedpe files in $bedpe_dir"

task_id=${SLURM_ARRAY_TASK_ID}

if (( task_id < 0 || task_id >= ${#bedpe_files[@]} )); then
    echo "No BEDPE file for SLURM_ARRAY_TASK_ID=${task_id}"
    exit 1
fi

# Select the BEDPE file for this array task
bedpe="${bedpe_files[$task_id]}"

# Extract sample name
sample=$(basename "$bedpe" .bedpe)

echo "SLURM_ARRAY_TASK_ID=${task_id}"
echo "Processing BEDPE: $bedpe"
echo "Sample: ${sample}"

# BED file format: fragment start/end interval/QNAME/frag size
# Columns from BEDPE:
#  1: chrom1, 2: start1, 6: end2, 7: name, 11: fragment length
awk 'OFS = "\t" {print $1, $2, $6, $7, $11}' "$bedpe" \
  | sort -k1,1 -k2,2n \
  | bedtools subtract -a - -b "$blacklist" -A \
  | bedtools subtract -a - -b "$gap" -A \
  | bedtools nuc -fi "$reference" -bed - \
  | awk 'OFS = "\t" {print $1, $2, $3, $4, $5, $7}' - \
  > "${output_dir}/${sample}_frags_gc.bed"

echo "${sample} complete."
