#!/bin/bash
#
#SBATCH --job-name=griffin_parallel
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 1-00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_nucleosome_profiling/griffin_controller_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_nucleosome_profiling/griffin_controller_%j.err

set -euo pipefail

module purge
module load gcc/8.2.0
module load r/4.4.0
module load python/ondemand-jupyter-python3.11

source activate /ihome/alee/rak373/.conda/envs/griffin_updated

export PATH=/ihome/alee/rak373/.conda/envs/griffin_updated/bin:$PATH

PROJECT_DIR=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_nucleosome_profiling

SNAKEFILE="${PROJECT_DIR}/griffin_nucleosome_profiling.snakefile"
SAMPLES_CONFIG="${PROJECT_DIR}/config/samples.GC.yaml"
SITES_CONFIG="${PROJECT_DIR}/config/sites.yaml"
MAIN_CONFIG="${PROJECT_DIR}/config/config.yaml"
CLUSTER_CONFIG="${PROJECT_DIR}/config/cluster_slurm.yaml"

cd "${PROJECT_DIR}"

echo "Python: $(which python)"
echo "Snakemake: $(which snakemake)"
python --version
snakemake --version

snakemake \
    --snakefile "${SNAKEFILE}" \
    --configfile "${SAMPLES_CONFIG}" \
    --configfile "${SITES_CONFIG}" \
    --configfile "${MAIN_CONFIG}" \
    --cluster-config "${CLUSTER_CONFIG}" \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete \
    --jobs 40 \
    --cluster \
    "sbatch \
        --mem={cluster.mem} \
        --time={cluster.time} \
        --cpus-per-task={cluster.ncpus} \
        --ntasks={cluster.ntasks} \
        --output={cluster.output} \
        --job-name={cluster.JobName}"