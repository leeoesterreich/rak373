#!/bin/bash
#SBATCH --job-name=bwa_mem_align
#SBATCH --output=slurm_log/align/align_%A_%a.out
#SBATCH --error=slurm_log/align/align_%A_%a.err
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-19

# Load required modules
module purge
module load gcc/8.2.0 samtools/1.22.1 bwa/0.7.17 qualimap/2.3

# ==== File paths ====
sample_sheet="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/3_Alignment/Sample_sheet.csv"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/3_Alignment/"
ref="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
mkdir -p "slurm_log/align"

# ==== Get line from sample sheet ====
line=$(tail -n +2 "$sample_sheet" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
IFS=',' read -r sample_name R1_path R2_path <<< "$line"

echo "[INFO] Processing sample: $sample_name"
echo "  R1: $R1_path"
echo "  R2: $R2_path"

# ==== Check input files ====
if [[ ! -f "$R1_path" || ! -f "$R2_path" ]]; then
    echo "[SKIP] Missing FASTQ files for $sample_name. Skipping..."
    exit 0
fi

# ==== Set output paths ====
out_dir="${output_dir}/${sample_name}"
bam_file="${out_dir}/${sample_name}.sorted.bam"
hash_file="${out_dir}/${sample_name}.md5"

mkdir -p "$out_dir"

# ==== Compute file hash for caching ====
new_hash=$(md5sum "$R1_path" "$R2_path" | md5sum | cut -d' ' -f1)

# ==== Skip if already aligned and hash matches ====
if [[ -f "$bam_file" && -f "${bam_file}.bai" && -f "$hash_file" ]]; then
    old_hash=$(cat "$hash_file")
    if [[ "$old_hash" == "$new_hash" ]]; then
        echo "[SKIP] Alignment already done and up-to-date for $sample_name."
        exit 0
    else
        echo "[WARN] Hash mismatch. Re-aligning $sample_name..."
    fi
else
    echo "[INFO] No previous alignment found. Starting $sample_name..."
fi

# ==== BWA-MEM2 Alignment ====
echo "[INFO] Running BWA MEM for $sample_name"
bwa mem -t $SLURM_CPUS_PER_TASK "$ref" "$R1_path" "$R2_path" \
    | samtools view -@ $SLURM_CPUS_PER_TASK -bS - \
    | samtools sort -@ $SLURM_CPUS_PER_TASK -o "$bam_file"

# ==== Index BAM ====
samtools index "$bam_file"

# ==== Save hash ====
echo "$new_hash" > "$hash_file"
echo "[INFO] Hash saved for $sample_name"

# ==== Qualimap QC ====
echo "[INFO] Running Qualimap QC"
qualimap bamqc -bam "$bam_file" \
               -outdir "${out_dir}/qualimap" \
               -outfile "${sample_name}_qualimap_report.pdf" \
               -outformat PDF

echo "Alignment + QC completed for $sample_name"
