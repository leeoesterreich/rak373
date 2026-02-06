#!/bin/bash
# Consolidate per-sample end motif GC-corrected RDS into one summary tibble

#SBATCH --job-name=EndMotif_combine
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=slurm_log/Endmotif/Endmotif_combine_%j.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_combine_%j.err
#SBATCH -t 06:00:00

set -euo pipefail

module purge
module load r/4.5.0

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif"

motifdir="${designated_path}/motif_gc"      # contains *_4bpmotif_gc.rds from previous step #5
outdir="${designated_path}/combined_motif"  # output folder for merged tibble

mkdir -p "$outdir" "slurm_log/Endmotif"

rscript="${designated_path}/combine_motif.R"
if [[ ! -s "$rscript" ]]; then
  echo "ERROR: R script not found/empty: $rscript"
  exit 1
fi

n_files=$(find "$motifdir" -maxdepth 1 -type f -name "*_4bpmotif_gc.rds" | wc -l)
echo "Found ${n_files} *_4bpmotif_gc.rds files in: $motifdir"
if [[ "$n_files" -eq 0 ]]; then
  echo "ERROR: No *_4bpmotif_gc.rds files found in $motifdir"
  exit 1
fi

Rscript "$rscript" \
  --motifdir "$motifdir" \
  --outdir "$outdir"

echo "DONE: combined motif summary written to $outdir"
