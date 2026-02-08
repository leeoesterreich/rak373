#!/bin/bash
#SBATCH --job-name=manta
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --array=0-7
#SBATCH --output=slurm_log/Manta/manta_run_%A_%a.out
#SBATCH --error=slurm_log/Manta/manta_run_%A_%a.err

set -euo pipefail

# ============================================================
# Load required modules
# ============================================================
module purge
module load gcc/8.2.0
module load manta/1.6.0

# path to manta container used by the module
manta_sif="/software/rhel9/manual/install/manta/manta-1.6.0.sif"

# ============================================================
# Define paths and constants (EnRich)
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"

pairs_csv="$designated_path/sample_sheet.csv"
bqsr_dir="$designated_path/3_BQSR_010526"

output_dir="$designated_path/Manta_010526"
log_dir="$designated_path/slurm_log/Manta"
mkdir -p "$output_dir" "$log_dir"

status_log="$output_dir/manta_status.log"

# ============================================================
# Get sample line for this array task (skip header)
# CSV header: Sample,BAM_Path
# ============================================================
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$pairs_csv" | tr -d '\r')
[[ -n "${line:-}" ]] || { echo "[INFO] No line for task ${SLURM_ARRAY_TASK_ID}; exiting."; exit 0; }

IFS=',' read -r SAMPLE _BAMPATH <<< "$line"
SAMPLE=$(echo "$SAMPLE" | xargs)

[[ -n "${SAMPLE:-}" ]] || { echo "[ERROR] Failed to parse SAMPLE from: $line" >&2; exit 1; }

# Use BQSR BAM
TUMOR_BAM_PATH="${bqsr_dir}/${SAMPLE}/${SAMPLE}.bqsr.bam"

echo "============================================================"
echo "[INFO] SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "[INFO] Starting sample: $SAMPLE"
echo "[INFO] Tumor BAM (BQSR): $TUMOR_BAM_PATH"
echo "============================================================"

sample_outdir="$output_dir/$SAMPLE"
run_dir="$sample_outdir/manta_run"
mkdir -p "$sample_outdir" "$run_dir"
cd "$sample_outdir" || exit 1

candidate_sv="${sample_outdir}/${SAMPLE}_manta_candidateSV.vcf.gz"
diploid_sv="${sample_outdir}/${SAMPLE}_manta_diploidSV.vcf.gz"
manta_log="${sample_outdir}/${SAMPLE}_manta.log"

if [[ -f "$candidate_sv" || -f "$diploid_sv" ]]; then
  echo "[SKIP] $SAMPLE already processed (VCF exists)"
  echo "$SAMPLE,SKIPPED (already exists),$(date)" >> "$status_log"
  echo ""
  exit 0
fi

echo "[INFO] Configuring Manta for $SAMPLE..."
start_time=$(date +%s)

# configManta.py provided by module (wrapper).
configManta.py \
  --bam "$TUMOR_BAM_PATH" \
  --referenceFasta "$reference" \
  --runDir "$run_dir" \
  > "$manta_log" 2>&1

echo "[INFO] Running Manta workflow for $SAMPLE..."

# runWorkflow.py inside manta container
singularity exec "$manta_sif" python \
  "$run_dir/runWorkflow.py" \
  -m local \
  -j 32 >> "$manta_log" 2>&1

status=$?
end_time=$(date +%s)
elapsed=$((end_time - start_time))

if [[ $status -ne 0 ]]; then
  echo "[ERROR] Manta failed for $SAMPLE"
  echo "$SAMPLE,FAILED,Elapsed:${elapsed}s,$(date)" >> "$status_log"
  echo ""
  exit 1
fi

manta_results_dir="$run_dir/results/variants"

# Copy results with sample-prefixed names
if [[ -f "${manta_results_dir}/candidateSV.vcf.gz" ]]; then
  cp "${manta_results_dir}/candidateSV.vcf.gz" "$candidate_sv"
  cp "${manta_results_dir}/candidateSV.vcf.gz.tbi" "${candidate_sv}.tbi"
fi

if [[ -f "${manta_results_dir}/diploidSV.vcf.gz" ]]; then
  cp "${manta_results_dir}/diploidSV.vcf.gz" "$diploid_sv"
  cp "${manta_results_dir}/diploidSV.vcf.gz.tbi" "${diploid_sv}.tbi"
fi

echo "[SUCCESS] Manta completed for $SAMPLE"
echo "$SAMPLE,SUCCESS,Elapsed:${elapsed}s,$(date)" >> "$status_log"
