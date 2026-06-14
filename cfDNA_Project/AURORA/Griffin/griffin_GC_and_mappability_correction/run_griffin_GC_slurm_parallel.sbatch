#!/bin/bash
#SBATCH --job-name=griffin_GCmap_AURORA
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -t 2-00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.%j.err

set -euo pipefail

module purge
module load gcc/8.2.0 r/4.4.0
module load python/ondemand-jupyter-python3.11

source activate /ihome/alee/rak373/.conda/envs/griffin_updated
export PATH=/ihome/alee/rak373/.conda/envs/griffin_updated/bin:$PATH

WORKDIR=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_GC_and_mappability_correction
CONFIGDIR=${WORKDIR}/config
SCRIPTSDIR=/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Griffin_scripts

cd ${WORKDIR}

mkdir -p logs/snakemake logs/cluster results
chmod +x ${SCRIPTSDIR}/*.py

echo "Python:"
which python
echo "Snakemake:"
which snakemake
snakemake --version

# Number of independent Slurm jobs Snakemake is allowed to submit at once.
JOBS=40

snakemake \
  -s ${WORKDIR}/griffin_GC_and_mappability_correction.snakefile \
  --configfile ${CONFIGDIR}/samples.yaml \
  --configfile ${CONFIGDIR}/config.yaml \
  --cluster-config ${CONFIGDIR}/cluster_slurm.yaml \
  --latency-wait 120 \
  --keep-going \
  --rerun-incomplete \
  --printshellcmds \
  --cluster "sbatch --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.ncpus} --ntasks={cluster.ntasks} --output={cluster.output} --job-name={cluster.JobName}" \
  --jobs ${JOBS}
