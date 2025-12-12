#!/bin/bash
#SBATCH --job-name=gridss
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 6-00:00:00
#SBATCH --output=slurm_log/Gridss/gridss_run_%j.out
#SBATCH --error=slurm_log/Gridss/gridss_run_%j.err

set -euo pipefail

# ============================================================
# Load required modules
# ============================================================
module purge
module load gcc/8.2.0
module load gridss/2.13.2

# ============================================================
# Define paths and constants
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/WGS"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
exclude_list="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Blacklist_file/hg38_unified_blacklist.bed"
gridss_jar="/software/rhel9/manual/install/gridss/2.13.2/python3.11/share/gridss-2.13.2-5/gridss.jar"
pairs_csv="$designated_path/Samples.csv"

output_dir="$designated_path/Gridss"
log_dir="$designated_path/slurm_log/Gridss"
mkdir -p "$output_dir" "$log_dir"

status_log="$output_dir/gridss_status.log"

# ============================================================
# Process each sample in CSV (tumor-only)
# Expecting header: TUMOR_BAM,TUMOR_BAM_PATH
# ============================================================
tail -n +2 "$pairs_csv" | while IFS=',' read -r TUMOR_BAM TUMOR_BAM_PATH
do
    SAMPLE="$TUMOR_BAM"

    echo "============================================================"
    echo "[INFO] Starting sample: $SAMPLE"
    echo "Tumor BAM: $TUMOR_BAM_PATH"
    echo "============================================================"

    sample_outdir="$output_dir/$SAMPLE"
    mkdir -p "$sample_outdir"
    cd "$sample_outdir" || exit 1

    # Validate input files
    for f in "$TUMOR_BAM_PATH" "$reference" "$gridss_jar"; do
        if [ ! -f "$f" ]; then
            echo "[ERROR] Missing file: $f"
            echo "$SAMPLE,FAILED (missing $f),$(date)" >> "$status_log"
            continue 2
        fi
    done

    # Check BAM index
    if [ ! -f "${TUMOR_BAM_PATH}.bai" ]; then
        echo "[ERROR] Missing BAM index: ${TUMOR_BAM_PATH}.bai"
        echo "$SAMPLE,FAILED (missing index),$(date)" >> "$status_log"
        continue
    fi

    output_vcf="${SAMPLE}_gridss.vcf"
    log_file="${SAMPLE}_gridss.log"

    if [ -f "$output_vcf" ]; then
        echo "[SKIP] $SAMPLE already processed"
        echo "$SAMPLE,SKIPPED (already exists),$(date)" >> "$status_log"
        continue
    fi

    echo "[INFO] Running GRIDSS..."
    start_time=$(date +%s)

    gridss \
        -r "$reference" \
        -j "$gridss_jar" \
        -o "$output_vcf" \
        -b "$exclude_list" \
        --threads 32 \
        "$TUMOR_BAM_PATH" > "$log_file" 2>&1

    status=$?
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    if [ $status -ne 0 ]; then
        echo "[ERROR] GRIDSS failed for $SAMPLE"
        echo "$SAMPLE,FAILED,Elapsed:${elapsed}s,$(date)" >> "$status_log"
    else
        echo "[SUCCESS] GRIDSS completed for $SAMPLE"
        echo "$SAMPLE,SUCCESS,Elapsed:${elapsed}s,$(date)" >> "$status_log"
    fi

    echo ""
done
