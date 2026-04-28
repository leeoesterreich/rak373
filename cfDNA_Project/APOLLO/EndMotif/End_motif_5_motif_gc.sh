#!/bin/bash
# Array-per-sample GC bias correction for end motifs
# Each array task processes ONE: motif_merge/*_4bp_motif.bed

#SBATCH --job-name=EndMotif_motifGC
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm_log/Endmotif/Endmotif_motifGC_%A_%a.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_motifGC_%A_%a.err
#SBATCH -t 1-00:00
#SBATCH --array=0-261%50   # <-- adjust if file count changes

set -euo pipefail

module purge
module load r/4.5.0

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif"

motifdir="${designated_path}/motif_merge"        # motif_merge
outdir="${designated_path}/motif_gc"             # 20-motif_gc
plotdir="${designated_path}/gcbias_plots"        # 20-gcbias_plots
statdir="${designated_path}/motif_gc_stats"      # 20-motif_gc_stats

mkdir -p "$outdir" "$plotdir" "$statdir" "slurm_log/Endmotif"

rscript="${designated_path}/20-motif_gc.R"
if [[ ! -s "$rscript" ]]; then
  echo "ERROR: R script not found/empty: $rscript"
  exit 1
fi

# Collect motif files
motif_files=( "$motifdir"/*_4bp_motif.bed )
if [[ ${#motif_files[@]} -eq 0 || ! -e "${motif_files[0]}" ]]; then
  echo "ERROR: No *_4bp_motif.bed files found in $motifdir"
  exit 1
fi

task_id=${SLURM_ARRAY_TASK_ID}
if (( task_id < 0 || task_id >= ${#motif_files[@]} )); then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${task_id} out of range (0..$(( ${#motif_files[@]} - 1 )))"
  exit 1
fi

f="${motif_files[$task_id]}"
id=$(basename "$f" "_4bp_motif.bed")

# Create per-sample plot folder (your R plotting path expects this)
mkdir -p "${plotdir}/${id}"

echo "Task ${task_id} / ${#motif_files[@]}"
echo "Sample: $id"
echo "File:   $f"

Rscript "$rscript" \
  --file "$f" \
  --outdir "$outdir" \
  --plotdir "$plotdir" \
  --statdir "$statdir" \
  --nmer 4 \
  --min_frag 90 \
  --max_frag 650

echo "DONE: ${id}"
