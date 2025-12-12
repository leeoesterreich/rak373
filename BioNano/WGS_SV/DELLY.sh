#!/bin/bash

# This script is used for structural variant discovery using Delly through integrated paired-end and split-read analysis.

#SBATCH --job-name=Delly
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 1-00:00
#SBATCH --array=1-4
#SBATCH --output=slurm_log/Delly/Delly_%A_%a.out
#SBATCH --error=slurm_log/Delly/Delly_%A_%a.err

# Load required modules
module load bcftools/1.22
module load delly/1.5.0

# ============================================================
# Define file paths and directories
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/WGS"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
pairs_csv="$designated_path/Samples.csv"
output_dir="$designated_path/Delly"
mkdir -p "$output_dir"

# ============================================================
# Get line corresponding to SLURM_ARRAY_TASK_ID
# ============================================================
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$pairs_csv")  # line 1 is header
IFS=',' read -r TUMOR_BAM TUMOR_BAM_PATH <<< "$line"

SAMPLE="$TUMOR_BAM"
# ============================================================
# Define output directory for this sample
# ============================================================
sample_outdir="$output_dir/$SAMPLE"
mkdir -p "$sample_outdir"
cd "$sample_outdir"

echo "[INFO] Processing sample: $SAMPLE"
echo "[INFO] Tumor BAM: $TUMOR_BAM_PATH"

# ============================================================
# Run DELLY for each SV type
# ============================================================
for svtype in DEL DUP INV BND INS; do
    echo "[INFO] Calling $svtype ..."
    delly call \
        -t $svtype \
        -g "$reference" \
        -o "${SAMPLE}_${svtype}.bcf" \
        "$TUMOR_BAM_PATH"

    # Convert BCF to VCF
    if [ -f "${SAMPLE}_${svtype}.bcf" ]; then
        bcftools view "${SAMPLE}_${svtype}.bcf" > "${SAMPLE}_${svtype}.vcf"
        if [ -s "${SAMPLE}_${svtype}.vcf" ]; then
            rm -f "${SAMPLE}_${svtype}.bcf"
            echo "[INFO] Converted and removed ${SAMPLE}_${svtype}.bcf"
        else
            echo "[WARNING] VCF file empty or missing for ${svtype}"
        fi
    else
        echo "[ERROR] Missing BCF for ${svtype}"
    fi
done

echo "[INFO] Finished all SV types for $SAMPLE"
