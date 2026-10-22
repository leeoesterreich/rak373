#!/bin/bash
#
#SBATCH --job-name=griffin_nucleosome_profiling
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -t 3-00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_nucleosome_profiling/griffin_nucleosome_profiling.snakefile.out

module purge
module load gcc/8.2.0 r/4.4.0
module load python/ondemand-jupyter-python3.11

source activate /ihome/alee/rak373/.conda/envs/griffin_updated
export PATH=/ihome/alee/rak373/.conda/envs/griffin_updated/bin:$PATH

chmod +x /ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Griffin_scripts/*.py

P=/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Griffin/griffin_nucleosome_profiling
cd "$P"

# Build one target file per sample from the sample names in samples.GC.yaml
# (sample names are the lines like "  AL1Y-NB1:" under "samples:")
grep -E '^  [A-Za-z0-9_-]+:$' config/samples.GC.yaml | tr -d ' :' \
  | awk '{print "results/"$1"/"$1".GC_corrected.coverage.tsv"}' > griffin_targets.txt
echo "Targets: $(wc -l < griffin_targets.txt)"

# Release any leftover lock
snakemake -s griffin_nucleosome_profiling.snakefile \
  --configfile $P/config/samples.GC.yaml \
  --configfile $P/config/sites.yaml \
  --configfile $P/config/config.yaml \
  --unlock

# Parallel run across samples: up to 20 jobs at a time
snakemake \
  -s griffin_nucleosome_profiling.snakefile \
  --configfile $P/config/samples.GC.yaml \
  --configfile $P/config/sites.yaml \
  --configfile $P/config/config.yaml \
  --cluster-config $P/config/cluster_slurm.yaml \
  --latency-wait 60 \
  --keep-going \
  --rerun-incomplete \
  --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -e {cluster.error} -J {cluster.JobName}" \
  -j 20 \
  $(cat griffin_targets.txt)