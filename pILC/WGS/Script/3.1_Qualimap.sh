#!/bin/bash
#SBATCH --job-name=qualimap_only
#SBATCH --output=slurm_log/qualimap/qualimap_%A_%a.out
#SBATCH --error=slurm_log/qualimap/qualimap_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=0-19

# ==== Load modules ====
module purge
module load gcc/8.2.0 samtools/1.22.1 qualimap/2.3

# ==== File paths ====
sample_sheet="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/3_Alignment/Sample_sheet.csv"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/3_Alignment/"
qualimap_base="${output_dir}/Qualimap"
mkdir -p "slurm_log/qualimap" "$qualimap_base"

# ==== Get sample info ====
line=$(tail -n +2 "$sample_sheet" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
IFS=',' read -r sample_name R1_path R2_path <<< "$line"

echo "[INFO] Running Qualimap for sample: $sample_name"

# ==== Set paths ====
bam_file="${output_dir}/${sample_name}/${sample_name}.sorted.bam"
out_dir="${qualimap_base}/${sample_name}"

# ==== Check BAM file ====
if [[ ! -f "$bam_file" ]]; then
    echo "[SKIP] BAM file not found for $sample_name. Skipping..."
    exit 0
fi

# ==== Run Qualimap ====
mkdir -p "$out_dir"

qualimap bamqc -bam "$bam_file" \
               -outdir "$out_dir" \
               -nt ${SLURM_CPUS_PER_TASK} \
               --java-mem-size=128G

# ==== Confirm output ====
if [[ -f "${out_dir}/genome_results.txt" ]]; then
    echo "[SUCCESS] Qualimap output generated at $out_dir"
else
    echo "[ERROR] Qualimap output missing for $sample_name"
    exit 1
fi

echo "[INFO] Qualimap completed for $sample_name"
