#!/bin/bash
#SBATCH --job-name=Sorted_bam
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm_log/Fragmentation/Sorted_bam_%A_%a.out
#SBATCH --error=slurm_log/Fragmentation/Sorted_bam_%A_%a.err
#SBATCH -t 1-00:00
#SBATCH --array=0-261%50  # 262 BAMs, max 50 jobs at once

module load samtools/1.22.1

designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST"
output_dir="$designated_path/Fragmentation_analysis/Regional_Fragment_analysis/Sorted_bam"
mkdir -p "$output_dir" slurm_log/Fragmentation

bam_list="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/bam_list.txt"
task_id=${SLURM_ARRAY_TASK_ID}
bam=$(sed -n "$((task_id+1))p" "$bam_list")

if [[ -z "$bam" ]]; then
  echo "No BAM file found for SLURM_ARRAY_TASK_ID=${task_id}"
  exit 1
fi

filename=$(basename "$bam")
sample=${filename%.bam}

echo "SLURM_ARRAY_TASK_ID=${task_id}"
echo "BAM file: ${bam}"
echo "Sample: ${sample}"

# 1) Filter (proper pairs, MAPQ>=30, remove secondary/QCfail/dups/supplementary via -F 3844)
filter_bam="${output_dir}/${sample}.filter.bam"
echo "Filtering -> ${filter_bam}"
samtools view -bh -f 2 -F 3844 -q 30 "$bam" > "$filter_bam"

# 2) Flagstat
stat_file="${output_dir}/${sample}.filter.flagstat.txt"
echo "Flagstat -> ${stat_file}"
samtools flagstat "$filter_bam" > "$stat_file"

# 3) Coordinate sort (IMPORTANT: no -n)
sorted_bam="${output_dir}/${sample}.filter.coordsort.bam"
echo "Coordinate sorting -> ${sorted_bam}"
samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 2G -o "$sorted_bam" "$filter_bam"

# 4) Index
echo "Indexing -> ${sorted_bam}.bai"
samtools index -@ ${SLURM_CPUS_PER_TASK} "$sorted_bam"

# Remove intermediate filtered unsorted bam to save space
rm -f "$filter_bam"

echo "Done: ${sample}"
