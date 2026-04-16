#!/bin/bash

# This script is used to annotate CNV by CNVkit.

#SBATCH --job-name=cnvkit
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 1-00:00:00
#SBATCH --array=0-39
#SBATCH --output=slurm_log/cnvkit_%A_%a.out
#SBATCH --error=slurm_log/cnvkit_%A_%a.err

set -euo pipefail

# Load required modules
module purge
module load cnvkit/0.9.12
module load python/ondemand-jupyter-python3.11

# Define paths and directories
input_dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/cna_seg_file_1000kb"
refFlatb38="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/refFlat.txt"
cnvkit_script="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/cnv_annotate.py"
dir="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit"
output_dir="$dir/Annotated_cnr_file"
temp_dir="$dir/temp"
log_dir="slurm_log"

# Create necessary directories
mkdir -p "$temp_dir" "$output_dir" "$log_dir"

# Activate conda environment
source activate /ihome/alee/rak373/.conda/envs/cnvkit_env

# Ensure the cnv_annotate.py script is executable
chmod +x "$cnvkit_script"

# Collect all input files safely
mapfile -t files < <(find "$input_dir" -maxdepth 1 -type f -name "*.cna.seg" | sort)

total_files=${#files[@]}
batch_size=20
start=$((SLURM_ARRAY_TASK_ID * batch_size))
end=$((start + batch_size))

echo "[INFO] Total files found      : $total_files"
echo "[INFO] Array task ID          : ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Processing file index  : $start to $((end - 1))"

# Exit cleanly if this array task has no files left
if (( start >= total_files )); then
    echo "[INFO] No files assigned to this array task. Exiting."
    exit 0
fi

# Process only this batch
for ((i=start; i<end && i<total_files; i++)); do
    cna_seg="${files[$i]}"

    # Extract sample name from the file name
    sample=$(basename "$cna_seg" .cna.seg)

    echo "[INFO] Processing sample: $sample"
    echo "[INFO] Input file: $cna_seg"

    # Temporary files for intermediate steps
    file1="$cna_seg"
    file2="$temp_dir/${sample}_edit1.cnr"
    file3="$temp_dir/${sample}_edit2.cnr"
    file4="$temp_dir/${sample}_edit3.cnr"

    # Step 1: Append 0 to the third column
    awk '{$3=$3" "0; print }' OFS='\t' "$file1" > "$file2"

    # Step 2: Modify header columns as needed
    awk 'NR==1 {gsub("chr","chromosome"); gsub(".*.logR","log2",$7); gsub("0","gene",$4); print} NR>1 {print}' OFS='\t' "$file2" > "$file3"

    # Step 3: Replace "0" with "-" in the 4th column (for non-header rows)
    awk 'NR>1 {gsub("0","-",$4)} {print}' OFS='\t' "$file3" > "$file4"

    # Annotate the CNV file with cnv_annotate.py
    python "$cnvkit_script" "$refFlatb38" "$file4" -o "$output_dir/${sample}.annotated.cnr"

    echo "[INFO] Annotated CNR file generated for $sample"
done

echo "[INFO] Array task ${SLURM_ARRAY_TASK_ID} completed successfully."