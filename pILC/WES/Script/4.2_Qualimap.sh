#!/bin/bash
#SBATCH --job-name=qualimap_only
#SBATCH --output=slurm_log/qualimap/qualimap_%A_%a.out
#SBATCH --error=slurm_log/qualimap/qualimap_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=0-47

set -euo pipefail

# ==== Load modules ====
module purge
module load gcc/8.2.0 samtools/1.22.1 qualimap/2.3

# ==== File paths ====
dedup_base="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/4.1_Dedup/"
qualimap_base="${dedup_base}/Qualimap"
Feature_BED="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"

mkdir -p "slurm_log/qualimap" "$qualimap_base"

# ==== Get dedup sample folders ====
mapfile -t sample_dirs < <(find "$dedup_base" -mindepth 1 -maxdepth 1 -type d ! -name "Qualimap" | sort)

sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "${sample_dir}" || ! -d "${sample_dir}" ]]; then
    echo "[ERROR] No valid sample directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

sample_name=$(basename "$sample_dir")

echo "[INFO] Running Qualimap for sample: $sample_name"

# ==== Set paths ====
bam_file="${dedup_base}/${sample_name}/${sample_name}.dedup.bam"
out_dir="${qualimap_base}/${sample_name}"

# ==== Check BAM file ====
if [[ ! -f "$bam_file" ]]; then
    echo "[SKIP] BAM file not found for $sample_name. Skipping..."
    exit 0
fi

if [[ ! -f "${bam_file}.bai" ]]; then
    echo "[SKIP] BAM index file not found for $sample_name. Skipping..."
    exit 0
fi

# ==== Prevent X11-related report failures in batch jobs ====
unset DISPLAY

# ==== Run Qualimap ====
qualimap bamqc \
    -bam "$bam_file" \
    -outdir "$out_dir" \
    -outformat HTML \
    -nt "${SLURM_CPUS_PER_TASK}" \
    --feature-file "${Feature_BED}" \
    --java-mem-size=128G

# ==== Confirm output ====
if [[ -f "${out_dir}/qualimapReport.html" && -d "${out_dir}/raw_data_qualimapReport" ]]; then
    echo "[SUCCESS] Full Qualimap HTML report generated at $out_dir"
elif [[ -f "${out_dir}/genome_results.txt" ]]; then
    echo "[WARN] Only genome_results.txt was generated for $sample_name"
    echo "[WARN] Check the SLURM log for report-generation errors"
    exit 1
else
    echo "[ERROR] Qualimap output missing for $sample_name"
    exit 1
fi

echo "[INFO] Qualimap completed for $sample_name"