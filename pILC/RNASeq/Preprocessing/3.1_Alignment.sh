#!/bin/bash
#SBATCH --job-name=Alignment_11samples
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00
#SBATCH --output=slurm_log/align/align_%A_%a.out
#SBATCH --error=slurm_log/align/align_%A_%a.err
#SBATCH --array=0-10

set -euo pipefail

module purge
module load gcc/8.2.0
module load star/2.7.11b

# -----------------------------
# Paths / directories
# -----------------------------
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
output_dir="${designated_path}/2_alignment_file"
mkdir -p "$output_dir" slurm_log/align

# Reference genome
genome_dir="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7"
gtf_file="${genome_dir}/gencode.v47.basic.annotation.gtf"   # kept for completeness (STAR index should already include it)

threads="${SLURM_CPUS_PER_TASK:-16}"

# Sample sheet (CSV)
sample_sheet="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/pILC_RNASeq_11samples.csv"

# -----------------------------
# Read the row for this task
# -----------------------------
if [ ! -f "$sample_sheet" ]; then
  echo "Sample sheet not found: $sample_sheet" >&2
  exit 1
fi

# Extract the (task_id+2)th line (skip header)
line=$(awk -v n=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==n {print; exit}' "$sample_sheet" || true)

if [ -z "${line:-}" ]; then
  echo "No line found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID in $sample_sheet" >&2
  echo "Check --array range and sample sheet rows." >&2
  exit 1
fi

# Parse CSV (assuming no commas inside fields)
IFS=',' read -r sample read1 read2 <<< "$line"

# Trim any accidental whitespace
sample="$(echo "$sample" | xargs)"
read1="$(echo "$read1" | xargs)"
read2="$(echo "$read2" | xargs)"

output_prefix="${output_dir}/${sample}_"
log_file="${output_dir}/${sample}.log"

# -----------------------------
# Input checks
# -----------------------------
if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
  echo "FASTQ files for sample $sample not found." | tee "$log_file"
  echo "R1: $read1" | tee -a "$log_file"
  echo "R2: $read2" | tee -a "$log_file"
  exit 1
fi

echo "Aligning ${sample} with STAR..." | tee "$log_file"
echo "R1: ${read1}" | tee -a "$log_file"
echo "R2: ${read2}" | tee -a "$log_file"
echo "Threads: ${threads}" | tee -a "$log_file"

STAR --genomeDir "$genome_dir" \
     --readFilesIn "$read1" "$read2" \
     --readFilesCommand zcat \
     --runThreadN "$threads" \
     --twopassMode Basic \
     --outFileNamePrefix "$output_prefix" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM \
     --quantMode TranscriptomeSAM GeneCounts \
     > "$log_file" 2>&1

echo "Sample $sample alignment complete." | tee -a "$log_file"