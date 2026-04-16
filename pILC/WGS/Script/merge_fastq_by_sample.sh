#!/bin/bash

# Merge all .fastq.gz fragments per sample

#SBATCH --job-name=MergeFASTQ
#SBATCH --output=slurm_log/merge_fastq/merge_%A_%a.out
#SBATCH --error=slurm_log/merge_fastq/merge_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --array=0-3  # Set based on number of unique samples

# Input/output
INPUT_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20230717_jis117_pILC_split_Jagmohan"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/Raw_fastq_split_sample"
mkdir -p "$OUTPUT_DIR" "slurm_log/merge_fastq"

# Generate list of unique sample IDs (once, indexed by array)
sample_list=($(cd "$INPUT_DIR" && ls *.fastq.gz | sed -E 's/_R[12]_001.*.fastq.gz$//' | sort -u))

# Get sample for this array task
sample="${sample_list[$SLURM_ARRAY_TASK_ID]}"
echo "[INFO] SLURM task ID: $SLURM_ARRAY_TASK_ID ? Processing sample: $sample"

# Merge R1
R1_files=$(ls "$INPUT_DIR/${sample}_R1_001"*.fastq.gz 2>/dev/null)
if [ -n "$R1_files" ]; then
    echo "[INFO] Merging R1 for $sample"
    cat $R1_files > "$OUTPUT_DIR/${sample}_R1.fastq.gz"
else
    echo "[WARN] No R1 files found for $sample"
fi

# Merge R2
R2_files=$(ls "$INPUT_DIR/${sample}_R2_001"*.fastq.gz 2>/dev/null)
if [ -n "$R2_files" ]; then
    echo "[INFO] Merging R2 for $sample"
    cat $R2_files > "$OUTPUT_DIR/${sample}_R2.fastq.gz"
else
    echo "[WARN] No R2 files found for $sample"
fi

echo "[DONE] Merging complete for sample: $sample"
