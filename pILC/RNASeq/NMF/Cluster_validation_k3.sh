#!/bin/bash
# This script is used to validate the NMF cluster K3. 
#SBATCH --job-name=pILC_Cluster_validation
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm_log/validation/validation_%j.out
#SBATCH --error=slurm_log/validation/validation_%j.err
#SBATCH --mem=32G
#SBATCH --time=12:00:00

module load r/4.6.0

Rscript /ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Cluster_validation_k3.R