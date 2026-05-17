#!/bin/bash
#SBATCH --job-name=gc_loess_crpb
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=1-00:00
#SBATCH --output=slurm_log/gc_loess_%j.out
#SBATCH --error=slurm_log/gc_loess_%j.err

set -euo pipefail

cd /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Fragmentation_analysis/Regional_Fragment_analysis/

source /ihome/alee/rak373/envs/gc_loess_env/bin/activate

python -u 05_gc_loess_correct_short_long.py
