#!/bin/bash
#SBATCH --job-name=FACTERA
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 1-00:00
#SBATCH --array=0-7
#SBATCH --output=slurm_log/FACTERA/FACTERA_%A_%a.out
#SBATCH --error=slurm_log/FACTERA/FACTERA_%A_%a.err

set -euo pipefail

# ============================================================
# Modules
# ============================================================
module load samtools/1.22.1
module load bedtools/2.31.1
module load kentutils/v479
module load blast-plus/2.14.1
module load perl/5.38.0
export PERL5LIB="$HOME/perl5/lib/perl5:${PERL5LIB:-}"
export PATH="$HOME/perl5/bin:${PATH}"

# ============================================================
# Define paths
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
reference="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/hg19.2bit"
exons_bed="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/FACTERA/factera-v1.4.4/exons.bed"
targets_bed="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Bed_file/Lee_2x_probes_Good.bed"
factera="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/FACTERA/factera-v1.4.4/factera.pl"

pairs_csv="${designated_path}/sample_sheet_010526.csv"
dedup_dir="${designated_path}/3_BQSR_010526"

output_dir="${designated_path}/FACTERA/010526"
mkdir -p "${output_dir}" "slurm_log/FACTERA"

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
TUMOR_BAM_PATH="${dedup_dir}/${SAMPLE}/${SAMPLE}.bqsr.bam"
[[ -f "$TUMOR_BAM_PATH" ]] || { echo "[ERROR] bqsr BAM not found: $TUMOR_BAM_PATH" >&2; exit 1; }

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

# ============================================================
# Required UCSC tool for hg19.2bit extraction
# ============================================================
if ! command -v twoBitToFa >/dev/null 2>&1; then
  echo "[ERROR] twoBitToFa not found in PATH. FACTERA requires UCSC twoBitToFa to read .2bit." >&2
  echo "        Try: module spider ucsc / kent / twoBitToFa, or add UCSC tools to PATH." >&2
  exit 1
fi
echo "[INFO] twoBitToFa: $(command -v twoBitToFa)"

# ============================================================
# Work from writable local copy in output directory
# ============================================================
work_bam="${sample_outdir}/${SAMPLE}.bam"
work_bai="${work_bam}.bai"

cp -f "$TUMOR_BAM_PATH" "$work_bam"
# copy whichever index exists
if [[ -f "${TUMOR_BAM_PATH}.bai" ]]; then
  cp -f "${TUMOR_BAM_PATH}.bai" "$work_bai"
elif [[ -f "${TUMOR_BAM_PATH%.bam}.bai" ]]; then
  cp -f "${TUMOR_BAM_PATH%.bam}.bai" "$work_bai"
fi

TUMOR_BAM_PATH="$work_bam"

# ============================================================
# Run FACTERA
# ============================================================
# FACTERA usage:
#   perl factera.pl [options] tumor.bam exons.bed hg19.2bit [optional: targets.bed]


perl "$factera" \
  "$TUMOR_BAM_PATH" \
  "$exons_bed" \
  "$reference" \
  "$targets_bed" \
  > "${SAMPLE}.factera.log" 2>&1

echo "[INFO] FACTERA finished for ${SAMPLE}"
echo "[INFO] Outputs in: $sample_outdir"
