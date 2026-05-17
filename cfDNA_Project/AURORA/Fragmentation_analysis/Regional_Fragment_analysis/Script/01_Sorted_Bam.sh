#!/bin/bash
#SBATCH --job-name=Sorted_bam
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm_log/Regional_Fragmentation/Sorted_bam_%A_%a.out
#SBATCH --error=slurm_log/Fragmentation/Sorted_bam_%A_%a.err
#SBATCH -t 2-00:00
#SBATCH --array=0-384%50  # 384 BAMs, max 50 jobs at once

module load samtools/1.22.1

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Fragmentation_analysis/Regional_Fragment_analysis"

output_dir="${designated_path}/Sorted_bam"
log_dir="${designated_path}/slurm_log/Regional_Fragmentation"

mkdir -p "$output_dir" "$log_dir"

bam_list="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Fragmentation_analysis/bam_samples_2.txt"

task_id=${SLURM_ARRAY_TASK_ID}
line_number=$((task_id + 1))

# Read sample ID and BAM path from the two-column file
sample=$(awk -v line="$line_number" 'NR==line {print $1}' "$bam_list")
bam=$(awk -v line="$line_number" 'NR==line {print $2}' "$bam_list")

# 1) Filter
# Keep proper pairs: -f 2
# Remove secondary, QC-fail, duplicates, supplementary: -F 3844
# Keep MAPQ >= 30: -q 30
filter_bam="${output_dir}/${sample}.filter.bam"

echo "Filtering -> ${filter_bam}"

samtools view \
  -@ "${SLURM_CPUS_PER_TASK}" \
  -bh \
  -f 2 \
  -F 3844 \
  -q 30 \
  "$bam" > "$filter_bam"

# 2) Flagstat
stat_file="${output_dir}/${sample}.filter.flagstat.txt"

echo "Flagstat -> ${stat_file}"

samtools flagstat \
  -@ "${SLURM_CPUS_PER_TASK}" \
  "$filter_bam" > "$stat_file"

# 3) Coordinate sort
sorted_bam="${output_dir}/${sample}.filter.coordsort.bam"

echo "Coordinate sorting -> ${sorted_bam}"

samtools sort \
  -@ "${SLURM_CPUS_PER_TASK}" \
  -m 2G \
  -o "$sorted_bam" \
  "$filter_bam"

# 4) Index
echo "Indexing -> ${sorted_bam}.bai"

samtools index \
  -@ "${SLURM_CPUS_PER_TASK}" \
  "$sorted_bam"

# Remove intermediate filtered unsorted BAM to save space
rm -f "$filter_bam"

echo "Done: ${sample}"
echo "Final BAM: ${sorted_bam}"
echo "Final index: ${sorted_bam}.bai"