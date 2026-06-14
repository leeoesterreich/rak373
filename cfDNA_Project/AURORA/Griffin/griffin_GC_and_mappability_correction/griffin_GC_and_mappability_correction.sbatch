#!/bin/bash

#SBATCH --job-name=griffin_GC_and_mappability_correction
#SBATCH -N 1
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction_snakefile.out

module purge
module load gcc/8.2.0 r/4.4.0
module load python/ondemand-jupyter-python3.11

source activate /ihome/alee/rak373/.conda/envs/griffin_updated

# To set the default Python interpreter
export PATH=/ihome/alee/rak373/.conda/envs/griffin_updated/bin:$PATH

# Grant execute permissions to all python scripts in the scripts directory
chmod +x /ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Griffin_scripts/*.py

which python
which snakemake


# /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.snakefile

snakemake -s griffin_GC_and_mappability_correction.snakefile --cores 1 #runs one job at a time

snakemake \
  -s griffin_GC_and_mappability_correction.snakefile \
  --configfile /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/config/samples.yaml \
  --configfile /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/config/config.yaml \
  --configfile /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/config/cluster_slurm.yaml \
  --latency-wait 60 \
  --keep-going \
  --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" \
  -j 80