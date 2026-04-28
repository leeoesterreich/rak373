#!/bin/bash
#SBATCH --job-name=qualimap_batch
#SBATCH --output=slurm_log/qualimap/qualimap_%A_%a.out
#SBATCH --error=slurm_log/qualimap/qualimap_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --array=0-261%50  # 262 BAMs, max 50 jobs at once

# ==== Load modules ====
module purge
module load gcc/8.2.0 samtools/1.22.1 qualimap/2.3

# ==== Paths ====
bam_dir="/ix1/alee/LO_LAB/General/Lab_Data/20250827_NST_ctDNA_shWGS_02_Daisong_Rahul/bam_file"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Qualimap"
bam_list="bam_list.txt"

mkdir -p slurm_log/qualimap "$output_dir"

# ==== Generate BAM list (once) ====
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    find "$bam_dir" -type f -name "*.bam" | sort > "$bam_list"
    echo "[INFO] BAM list created with $(wc -l < "$bam_list") files"
fi
sleep 5  # wait for file creation across nodes

# ==== Read BAM file ====
bam_file=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$bam_list")
sample_name=$(basename "$bam_file" .bam)
out_dir="${output_dir}/${sample_name}"

echo "[INFO] Processing $sample_name"
echo "[INFO] Input: $bam_file"
echo "[INFO] Output: $out_dir"

# ==== Check if BAM exists ====
if [[ ! -f "$bam_file" ]]; then
    echo "[ERROR] BAM file not found: $bam_file"
    exit 1
fi

# ==== Run Qualimap ====
mkdir -p "$out_dir"
qualimap bamqc -bam "$bam_file" \
               -outdir "$out_dir" \
               -nt ${SLURM_CPUS_PER_TASK} \
               --java-mem-size=128G

# ==== Check output ====
if [[ -f "${out_dir}/genome_results.txt" ]]; then
    echo "[SUCCESS] Completed $sample_name"
else
    echo "[ERROR] Qualimap failed for $sample_name"
    exit 1
fi
