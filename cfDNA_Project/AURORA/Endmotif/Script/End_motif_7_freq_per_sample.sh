#!/bin/bash
# Compute per-sample 4bp end-motif frequencies from combined motif RDS.
# Runs both:
#   1) full motif summary
#   2) ATGC-only motif summary

#SBATCH --job-name=EndMotif_freq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm_log/Endmotif/Endmotif_freq_%j.out
#SBATCH --error=slurm_log/Endmotif_freq_%j.err
#SBATCH -t 04:00:00

set -euo pipefail

module purge
module load r/4.6.0

export R_LIBS_USER="/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.6"

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Endmotif"

combined_dir="${designated_path}/combined_motif"
outdir="${designated_path}/motif_analysis"
logdir="${designated_path}/slurm_log/Endmotif"

rscript="${designated_path}/End_motif_7_freq_only.R"

mkdir -p "$outdir" "$logdir"

echo "============================================================"
echo "End motif frequency job started"
echo "Date: $(date)"
echo "Input combined dir: $combined_dir"
echo "Output dir:         $outdir"
echo "R script:           $rscript"
echo "============================================================"

if [[ ! -s "$rscript" ]]; then
  echo "ERROR: R script not found or empty:"
  echo "$rscript"
  exit 1
fi

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
in_full="${combined_dir}/endmotif_4bp_gc_summary.rds"
in_atgc="${combined_dir}/endmotif_4bp_gc_summary_ATGC_only.rds"

if [[ ! -s "$in_full" ]]; then
  echo "ERROR: Full motif input RDS missing or empty:"
  echo "$in_full"
  exit 1
fi

if [[ ! -s "$in_atgc" ]]; then
  echo "ERROR: ATGC-only motif input RDS missing or empty:"
  echo "$in_atgc"
  exit 1
fi

# ------------------------------------------------------------
# Run full motif frequency
# ------------------------------------------------------------
echo "============================================================"
echo "Running full motif frequency calculation"
echo "Input: $in_full"
echo "============================================================"

Rscript "$rscript" \
  --in_rds "$in_full" \
  --outdir "$outdir" \
  --prefix "endmotif_4bp_freq_per_sample"

# ------------------------------------------------------------
# Run ATGC-only motif frequency
# ------------------------------------------------------------
echo "============================================================"
echo "Running ATGC-only motif frequency calculation"
echo "Input: $in_atgc"
echo "============================================================"

Rscript "$rscript" \
  --in_rds "$in_atgc" \
  --outdir "$outdir" \
  --prefix "endmotif_4bp_freq_per_sample_ATGC_only"

echo "============================================================"
echo "DONE: motif frequency outputs written to:"
echo "$outdir"
echo "Job completed at: $(date)"
echo "============================================================"