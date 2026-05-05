#!/bin/bash
#SBATCH --job-name=bwa_mem_align
#SBATCH --output=slurm_log/align/align_%A_%a.out
#SBATCH --error=slurm_log/align/align_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-47

set -euo pipefail

# Load required modules
module purge
module load gcc/8.2.0 samtools/1.22.1 bwa/0.7.17

# ==== File paths ====
trim_base="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/2_TrimGalore"
output_base="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/3_Alignment"
ref="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"

mkdir -p "$output_base"
mkdir -p slurm_log/align

# ==== Get trimmed sample folders ====
mapfile -t sample_dirs < <(find "$trim_base" -mindepth 1 -maxdepth 1 -type d ! -name "Fastqc" | sort)

sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "${sample_dir}" || ! -d "${sample_dir}" ]]; then
    echo "[ERROR] No valid sample directory found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

sample_name=$(basename "$sample_dir")

# ==== Detect trimmed FASTQ files ====
R1_path=$(find "$sample_dir" -maxdepth 1 -type f -name "*_val_1.fq.gz" | sort | head -n 1)
R2_path=$(find "$sample_dir" -maxdepth 1 -type f -name "*_val_2.fq.gz" | sort | head -n 1)

echo "[INFO] Processing sample: $sample_name"
echo "  R1: ${R1_path:-NOT FOUND}"
echo "  R2: ${R2_path:-NOT FOUND}"

# ==== Check input files ====
if [[ -z "${R1_path}" || ! -f "${R1_path}" ]]; then
    echo "[ERROR] Missing R1 trimmed FASTQ for $sample_name"
    exit 1
fi

if [[ -z "${R2_path}" || ! -f "${R2_path}" ]]; then
    echo "[ERROR] Missing R2 trimmed FASTQ for $sample_name"
    exit 1
fi

# ==== Set output paths ====
out_dir="${output_base}/${sample_name}"
bam_file="${out_dir}/${sample_name}.sorted.bam"
hash_file="${out_dir}/${sample_name}.md5"

mkdir -p "$out_dir"

# ==== Compute file hash for caching ====
new_hash=$(md5sum "$R1_path" "$R2_path" | md5sum | cut -d' ' -f1)

# ==== Skip if already aligned and hash matches ====
if [[ -f "$bam_file" && -f "${bam_file}.bai" && -f "$hash_file" ]]; then
    old_hash=$(cat "$hash_file")
    if [[ "$old_hash" == "$new_hash" ]]; then
        echo "[SKIP] Alignment already done and up-to-date for $sample_name"
        exit 0
    else
        echo "[WARN] Hash mismatch. Re-aligning $sample_name"
    fi
else
    echo "[INFO] No previous alignment found. Starting $sample_name"
fi

# ==== BWA-MEM alignment ====
echo "[INFO] Running BWA MEM for $sample_name"

bwa mem -t "$SLURM_CPUS_PER_TASK" "$ref" "$R1_path" "$R2_path" \
    | samtools view -@ "$SLURM_CPUS_PER_TASK" -bS - \
    | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$bam_file"

# ==== Index BAM ====
samtools index -@ "$SLURM_CPUS_PER_TASK" "$bam_file"

# ==== Save hash ====
echo "$new_hash" > "$hash_file"
echo "[INFO] Hash saved for $sample_name"
echo "[SUCCESS] Alignment complete for $sample_name"