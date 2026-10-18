#!/bin/bash
#SBATCH --job-name=RSEM_matrix
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH -t 02:00:00
#SBATCH --output=slurm_log/RSEM/RSEM_matrix_%j.out
#SBATCH --error=slurm_log/RSEM/RSEM_matrix_%j.err

# Generate the gene-level expected count matrix after the complete execution of RSEM

set -euo pipefail
module purge
module load rsem/1.3.3

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
output_dir="${designated_path}/RSEM_Count"

cd "$output_dir"

# Safety check
n=$(ls -1 *.genes.results 2>/dev/null | wc -l)
if [ "$n" -eq 0 ]; then
  echo "ERROR: No *.genes.results found in $output_dir"
  exit 1
fi

rsem-generate-data-matrix *.genes.results > RSEM_expected_counts_genes.tsv
echo "Wrote: $output_dir/RSEM_expected_counts_genes.tsv (from $n samples)"
