#!/bin/bash
#SBATCH --job-name=BreakId
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 1-00:00
#SBATCH --array=0-7
#SBATCH --output=slurm_log/BreakId/BreakId_%A_%a.out
#SBATCH --error=slurm_log/BreakId/BreakId_%A_%a.err

set -euo pipefail

module purge
module load breakid/dev
module load samtools/1.22.1

# ============================================================
# Define paths
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
nib_dir="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/hg19_nib_files"

pairs_csv="${designated_path}/sample_sheet_012425.csv"
bqsr_dir="${designated_path}/3_BQSR_012425"

output_dir="${designated_path}/BreakId/012425"
mkdir -p "${output_dir}" "slurm_log/BreakId"

# ============================================================
# Sanity checks
# ============================================================
[[ -d "$nib_dir" ]] || { echo "[ERROR] nib_dir not found: $nib_dir" >&2; exit 1; }
[[ -f "$nib_dir/chr1.nib" ]] || { echo "[ERROR] chr1.nib not found in: $nib_dir" >&2; exit 1; }
[[ -f "$pairs_csv" ]] || { echo "[ERROR] sample sheet not found: $pairs_csv" >&2; exit 1; }

# ============================================================
# BreakID requires ref_names.txt inside nib_dir
# (create it automatically if missing)
# ============================================================
if [[ ! -s "$nib_dir/ref_names.txt" ]]; then
  echo "[INFO] Creating $nib_dir/ref_names.txt from nib files"
  ls "$nib_dir"/*.nib | xargs -n1 basename | sed 's/\.nib$//' | sort > "$nib_dir/ref_names.txt"
fi
[[ -s "$nib_dir/ref_names.txt" ]] || { echo "[ERROR] ref_names.txt missing/empty in: $nib_dir" >&2; exit 1; }

# ============================================================
# Get sample for this array task (skip header)
# ============================================================
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$pairs_csv" | tr -d '\r')
[[ -n "${line:-}" ]] || { echo "[INFO] No line for task ${SLURM_ARRAY_TASK_ID}; exiting."; exit 0; }

IFS=',' read -r SAMPLE _REST <<< "$line"
SAMPLE=$(echo "$SAMPLE" | xargs)
[[ -n "${SAMPLE:-}" ]] || { echo "[ERROR] Failed to parse SAMPLE from: $line" >&2; exit 1; }

# ============================================================
# Input BAM
# ============================================================
TUMOR_BAM_PATH="${bqsr_dir}/${SAMPLE}/${SAMPLE}.bqsr.bam"
[[ -f "$TUMOR_BAM_PATH" ]] || { echo "[ERROR] BAM not found: $TUMOR_BAM_PATH" >&2; exit 1; }

# Index could be .bai OR .bam.bai depending on how you made it
if [[ ! -f "${TUMOR_BAM_PATH}.bai" && ! -f "${TUMOR_BAM_PATH%.bam}.bai" ]]; then
  echo "[ERROR] BAM index not found for: $TUMOR_BAM_PATH" >&2
  exit 1
fi

# ============================================================
# Output per sample
# ============================================================
sample_outdir="${output_dir}/${SAMPLE}"
mkdir -p "$sample_outdir"
cd "$sample_outdir"

echo "[INFO] Processing sample: $SAMPLE"
echo "[INFO] BAM: $TUMOR_BAM_PATH"
echo "[INFO] NIB_DIR: $nib_dir"
echo "[INFO] OUTDIR: $sample_outdir"

# ============================================================
# Pre-check: ensure every BAM contig has a matching .nib
# ============================================================
missing=$(samtools idxstats "$TUMOR_BAM_PATH" | cut -f1 | grep -v '^*$' | while read -r c; do
  [[ -f "$nib_dir/${c}.nib" ]] || echo "$c"
done | head -n 20)

if [[ -n "$missing" ]]; then
  echo "[ERROR] Missing .nib files for BAM contigs (showing up to 20):"
  echo "$missing"
  echo "[ERROR] Fix nibs or filter BAM before running BreakID."
  exit 1
fi

# ============================================================
# Run BreakID
# ============================================================
out_prefix="${SAMPLE}.breakid"

BreakID \
  -i "$TUMOR_BAM_PATH" \
  -o "$out_prefix" \
  -n "$nib_dir" \
  -q 20 \
  -fast \
  2>&1 | tee "${out_prefix}.breakid.log"

echo "[INFO] Finished sample: $SAMPLE"
echo "[INFO] Log: ${sample_outdir}/${out_prefix}.breakid.log"
