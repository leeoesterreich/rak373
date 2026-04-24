#!/bin/bash

#SBATCH --job-name=run_GISTIC
#SBATCH --output=slurm_log/run_GISTIC_%j.out
#SBATCH --error=slurm_log/run_GISTIC_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=120G

# This script runs GISTIC2.0 to identify significant copy number alterations.
# ref_gene_file can be downloaded from https://github.com/bzhanglab/GISTIC2_example/blob/master/hg38.UCSC.add_miR.160920.refgene.mat
# Template provided by Brent

# Load required modules
module purge
module load singularity/4.3.2
module load gistic/2.0.23

# --- Configuration ---
# Define file paths
seg_file="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/GISTIC/gistic_input_fixed.seg"
ref_gene_file="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/GISTIC/hg38.UCSC.add_miR.160920.refgene.mat"

# Create output directory
output_dir_base="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/GISTIC/output"
mkdir -p $output_dir_base

# --- Run GISTIC2 ---
echo "Starting GISTIC2 analysis..."

gistic2 -b $output_dir_base \
            -seg $seg_file \
            -refgene $ref_gene_file \
            -ta 0.3 \
            -td 0.3 \
            -conf 0.95 \
            -qvt 0.25 \
            -brlen 0.98 \
            -genegistic 1 \
            -armpeel 1 \
            -broad 1 \
            -savegene 1 \
            -gcm extreme \
            -v 30

# Print when gistic is finished
echo "GISTIC2 analysis completed. Results are in ${output_dir_base}" 
