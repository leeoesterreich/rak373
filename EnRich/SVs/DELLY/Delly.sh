#!/bin/bash
#SBATCH --job-name=Delly
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 1-00:00
#SBATCH --array=0-7
#SBATCH --output=slurm_log/Delly/Delly_%A_%a.out
#SBATCH --error=slurm_log/Delly/Delly_%A_%a.err

set -euo pipefail

module load bcftools/1.22
module load delly/1.5.0

# ============================================================
# Define file paths and directories
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"

pairs_csv="$designated_path/sample_sheet.csv"
bqsr_dir="$designated_path/3_BQSR_010526"

output_dir="$designated_path/Delly_010526"
mkdir -p "$output_dir" slurm_log/Delly

# ============================================================
# Get sample corresponding to SLURM_ARRAY_TASK_ID (skip header)
# ============================================================
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$pairs_csv" | tr -d '\r')
[[ -n "${line:-}" ]] || { echo "[INFO] No line for task ${SLURM_ARRAY_TASK_ID}; exiting."; exit 0; }

IFS=',' read -r SAMPLE _BAMPATH <<< "$line"
SAMPLE=$(echo "$SAMPLE" | xargs)

[[ -n "${SAMPLE:-}" ]] || { echo "[ERROR] Failed to parse SAMPLE from: $line" >&2; exit 1; }

# ============================================================
# Use BQSR BAM as input
# ============================================================
TUMOR_BAM_PATH="${bqsr_dir}/${SAMPLE}/${SAMPLE}.bqsr.bam"

[[ -f "$TUMOR_BAM_PATH" ]] || { echo "[ERROR] BQSR BAM not found: $TUMOR_BAM_PATH" >&2; exit 1; }
[[ -f "${TUMOR_BAM_PATH}.bai" ]] || { echo "[ERROR] BQSR BAM index not found: ${TUMOR_BAM_PATH}.bai" >&2; exit 1; }

# ============================================================
# Define output directory for this sample
# ============================================================
sample_outdir="$output_dir/$SAMPLE"
mkdir -p "$sample_outdir"
cd "$sample_outdir"

echo "[INFO] Processing sample: $SAMPLE"
echo "[INFO] Tumor BAM (BQSR): $TUMOR_BAM_PATH"

# ============================================================
# Run DELLY for each SV type
# ============================================================
for svtype in DEL DUP INV BND INS; do
  echo "[INFO] Calling $svtype ..."
  delly call \
    -t "$svtype" \
    -g "$reference" \
    -o "${SAMPLE}_${svtype}.bcf" \
    "$TUMOR_BAM_PATH"

  # Convert BCF to VCF and clean up
  if [[ -f "${SAMPLE}_${svtype}.bcf" ]]; then
    bcftools view "${SAMPLE}_${svtype}.bcf" > "${SAMPLE}_${svtype}.vcf"
    if [[ -s "${SAMPLE}_${svtype}.vcf" ]]; then
      rm -f "${SAMPLE}_${svtype}.bcf"
      echo "[INFO] Converted and removed ${SAMPLE}_${svtype}.bcf"
    else
      echo "[WARNING] VCF empty for ${svtype} (${SAMPLE})"
    fi
  else
    echo "[ERROR] Missing BCF for ${svtype} (${SAMPLE})"
  fi
done

echo "[INFO] Finished all SV types for $SAMPLE"
