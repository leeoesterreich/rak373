#!/bin/bash
#SBATCH --job-name=manta
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --array=1-4
#SBATCH --output=slurm_log/Manta/manta_run_%A_%a.out
#SBATCH --error=slurm_log/Manta/manta_run_%A_%a.err

set -euo pipefail

# ============================================================
# Load required modules
# ============================================================
module purge
module load gcc/8.2.0
module load manta/1.6.0      # defines configManta.py as a Singularity wrapper

# path to manta container used by the module
manta_sif="/software/rhel9/manual/install/manta/manta-1.6.0.sif"

# ============================================================
# Define paths and constants
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/WGS"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
pairs_csv="$designated_path/Samples.csv"

output_dir="$designated_path/Manta"
log_dir="$designated_path/slurm_log/Manta"
mkdir -p "$output_dir" "$log_dir"

status_log="$output_dir/manta_status.log"

# ============================================================
# Get the sample line for this array task
# Header: TUMOR_BAM,TUMOR_BAM_PATH
# SLURM_ARRAY_TASK_ID = 1..N  ? line = ID+1 (skip header)
# ============================================================
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$pairs_csv")
IFS=',' read -r TUMOR_BAM TUMOR_BAM_PATH <<< "$line"
SAMPLE="$TUMOR_BAM"

echo "============================================================"
echo "[INFO] SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "[INFO] Starting sample: $SAMPLE"
echo "Tumor BAM: $TUMOR_BAM_PATH"
echo "============================================================"

sample_outdir="$output_dir/$SAMPLE"
run_dir="$sample_outdir/manta_run"
mkdir -p "$sample_outdir" "$run_dir"
cd "$sample_outdir" || exit 1

# ---------- validate inputs ----------
missing="false"
for f in "$TUMOR_BAM_PATH" "$reference"; do
    if [ ! -f "$f" ]; then
        echo "[ERROR] Missing file: $f"
        missing="true"
    fi
done

if [ ! -f "${TUMOR_BAM_PATH}.bai" ]; then
    echo "[ERROR] Missing BAM index: ${TUMOR_BAM_PATH}.bai"
    missing="true"
fi

if [ "$missing" = "true" ]; then
    echo "$SAMPLE,FAILED (missing input),$(date)" >> "$status_log"
    echo ""
    exit 1
fi

candidate_sv="${sample_outdir}/${SAMPLE}_manta_candidateSV.vcf.gz"
diploid_sv="${sample_outdir}/${SAMPLE}_manta_diploidSV.vcf.gz"
manta_log="${sample_outdir}/${SAMPLE}_manta.log"

if [ -f "$candidate_sv" ] || [ -f "$diploid_sv" ]; then
    echo "[SKIP] $SAMPLE already processed (VCF exists)"
    echo "$SAMPLE,SKIPPED (already exists),$(date)" >> "$status_log"
    echo ""
    exit 0
fi

echo "[INFO] Configuring Manta for $SAMPLE..."
start_time=$(date +%s)

# configManta.py is a shell function that already wraps singularity exec
configManta.py \
    --bam "$TUMOR_BAM_PATH" \
    --referenceFasta "$reference" \
    --runDir "$run_dir" \
    > "$manta_log" 2>&1

echo "[INFO] Running Manta workflow for $SAMPLE..."

# IMPORTANT: run runWorkflow.py inside the same Singularity container
singularity exec "$manta_sif" python \
    "$run_dir/runWorkflow.py" \
    -m local \
    -j 32 >> "$manta_log" 2>&1

status=$?
end_time=$(date +%s)
elapsed=$((end_time - start_time))

if [ $status -ne 0 ]; then
    echo "[ERROR] Manta failed for $SAMPLE"
    echo "$SAMPLE,FAILED,Elapsed:${elapsed}s,$(date)" >> "$status_log"
    echo ""
    exit 1
fi

manta_results_dir="$run_dir/results/variants"

if [ -f "${manta_results_dir}/candidateSV.vcf.gz" ]; then
    cp "${manta_results_dir}/candidateSV.vcf.gz" "$candidate_sv"
    cp "${manta_results_dir}/candidateSV.vcf.gz.tbi" "${candidate_sv}.tbi"
fi

if [ -f "${manta_results_dir}/diploidSV.vcf.gz" ]; then
    cp "${manta_results_dir}/diploidSV.vcf.gz" "$diploid_sv"
    cp "${manta_results_dir}/diploidSV.vcf.gz.tbi" "${diploid_sv}.tbi"
fi

echo "[SUCCESS] Manta completed for $SAMPLE"
echo "$SAMPLE,SUCCESS,Elapsed:${elapsed}s,$(date)" >> "$status_log"
echo ""
