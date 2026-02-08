#!/bin/bash
#SBATCH --job-name=SViCT
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 1-00:00
#SBATCH --array=0-7
#SBATCH --output=slurm_log/SViCT/SViCT_%A_%a.out
#SBATCH --error=slurm_log/SViCT/SViCT_%A_%a.err

set -euo pipefail

module purge
module load python/ondemand-jupyter-python3.11

# To set the default Python interpreter
export PATH=/ihome/alee/rak373/.conda/envs/SViCT/bin:$PATH

# sanity checks (good for logs)
echo "[INFO] python: $(which python)"
echo "[INFO] python version: $(python --version)"
echo "[INFO] svict: $(which svict)"
# ============================================================
# Define paths
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"

pairs_csv="${designated_path}/sample_sheet_010526.csv"
bqsr_dir="${designated_path}/3_BQSR_010526"

output_dir="${designated_path}/SViCT/010526"
mkdir -p "${output_dir}" "slurm_log/SViCT"

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
[[ -f "$TUMOR_BAM_PATH" ]] || { echo "[ERROR] BQSR BAM not found: $TUMOR_BAM_PATH" >&2; exit 1; }

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
echo "[INFO] Tumor BAM (BQSR): $TUMOR_BAM_PATH"
echo "[INFO] Output dir: $sample_outdir"

# ============================================================
# Run SViCT
# NOTE: -o is a PREFIX; SViCT will write: <prefix>.vcf
# ============================================================
out_prefix="${sample_outdir}/${SAMPLE}.svict"

svict \
  -i "$TUMOR_BAM_PATH" \
  -r "$reference" \
  -o "$out_prefix"

vcf="${out_prefix}.vcf"
[[ -s "$vcf" ]] || { echo "[ERROR] Expected VCF not created or empty: $vcf" >&2; exit 1; }

# ============================================================
# Filter: keep header + records with SUPPORT > 3
# Output: <prefix>_filtered.vcf
# ============================================================
filtered_vcf="${sample_outdir}/${SAMPLE}_filtered.vcf"

awk 'BEGIN{FS=OFS="\t"}
     /^#/ {print; next}
     {
       sup=0
       if (match($8, /SUPPORT=([0-9]+)/, a)) sup=a[1]
       if (sup > 3) print
     }' "$vcf" > "$filtered_vcf"

echo "[INFO] Wrote: $vcf"
echo "[INFO] Wrote: $filtered_vcf"
echo "[INFO] Done for $SAMPLE"
