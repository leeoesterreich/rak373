#!/bin/bash
# Consolidate per-sample end motif GC-corrected RDS files
# Produces:
#   endmotif_4bp_gc_summary.rds
#   endmotif_4bp_gc_summary.tsv.gz
#   endmotif_4bp_gc_summary_ATGC_only.rds
#   endmotif_4bp_gc_summary_ATGC_only.tsv.gz

#SBATCH --job-name=EndMotif_combine
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=slurm_log/Endmotif/Endmotif_combine_%j.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_combine_%j.err
#SBATCH -t 01:00:00

set -euo pipefail

module purge
module load r/4.6.0

export R_LIBS_USER="/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.6"

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Endmotif"

motifdir="${designated_path}/motif_gc"
outdir="${designated_path}/combined_motif"
logdir="${designated_path}/slurm_log/Endmotif"
rscript="${designated_path}/End_motif_6_combine_motif.R"

mkdir -p "$outdir" "$logdir"

echo "============================================================"
echo "End motif combine job started"
echo "Date: $(date)"
echo "Input motifdir: $motifdir"
echo "Output outdir:  $outdir"
echo "R script:       $rscript"
echo "============================================================"

if [[ ! -d "$motifdir" ]]; then
  echo "ERROR: motifdir does not exist:"
  echo "$motifdir"
  exit 1
fi

if [[ ! -s "$rscript" ]]; then
  echo "ERROR: R script not found or empty:"
  echo "$rscript"
  exit 1
fi

n_files=$(find "$motifdir" -maxdepth 1 -type f -name "*_4bpmotif_gc.rds" | wc -l)

echo "Found ${n_files} *_4bpmotif_gc.rds files in: $motifdir"

if [[ "$n_files" -eq 0 ]]; then
  echo "ERROR: No *_4bpmotif_gc.rds files found in $motifdir"
  exit 1
fi

echo "============================================================"
echo "Running R script now..."
echo "============================================================"

Rscript "$rscript" \
  --motifdir "$motifdir" \
  --outdir "$outdir"

echo "============================================================"
echo "Checking expected output files..."
echo "============================================================"

expected_files=(
  "endmotif_4bp_gc_summary.tsv.gz"
  "endmotif_4bp_gc_summary.rds"
  "endmotif_4bp_gc_summary_ATGC_only.tsv.gz"
  "endmotif_4bp_gc_summary_ATGC_only.rds"
)

for ff in "${expected_files[@]}"; do
  if [[ -s "${outdir}/${ff}" ]]; then
    echo "FOUND: ${outdir}/${ff}"
    ls -lh "${outdir}/${ff}"
  else
    echo "MISSING or EMPTY: ${outdir}/${ff}"
    exit 1
  fi
done

echo "============================================================"
echo "DONE: combined motif summaries written to:"
echo "$outdir"
echo "Job completed at: $(date)"
echo "============================================================"