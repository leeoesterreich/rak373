#!/bin/bash
# This is a continuation script (III).
# Use GC-filtered fragments (03-filter_frags_gc.sh output) to
# obtain paired-end alignments associated with filtered fragments.

#SBATCH --job-name=End_motif_analysis_step3
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --output=slurm_log/Endmotif/Endmotif_step3_%A_%a.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_step3_%A_%a.err
#SBATCH -t 1-00:00
#SBATCH --array=0-399%50  

# Define input/output directories
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Endmotif"
bedpe_dir="${designated_path}/bamtobed"          # 02-bamtobed
frags_dir="${designated_path}/filter_frags_gc"   # 03-filter_frags_gc output (file1)
output_dir="${designated_path}/filter_bedpe"     # 05-filter_bedpe
mkdir -p "$output_dir"

# Collect all BEDPE files into an array
bedpe_files=( "$bedpe_dir"/*.bedpe )

# Safety: check that we actually found files
if [[ ${#bedpe_files[@]} -eq 0 || ! -e "${bedpe_files[0]}" ]]; then
    echo "No .bedpe files found in $bedpe_dir"
    exit 1
fi

echo "Found ${#bedpe_files[@]} .bedpe files in $bedpe_dir"

task_id=${SLURM_ARRAY_TASK_ID}

if (( task_id < 0 || task_id >= ${#bedpe_files[@]} )); then
    echo "No BEDPE file for SLURM_ARRAY_TASK_ID=${task_id}"
    exit 1
fi

# Select the BEDPE file for this array task
bedpe="${bedpe_files[$task_id]}"

# Extract sample name
sample=$(basename "$bedpe" .bedpe)

echo "SLURM_ARRAY_TASK_ID=${task_id}"
echo "Processing BEDPE: $bedpe"
echo "Sample: ${sample}"

## To obtain paired-end alignments associated with filtered fragments.
## file1 = ${sample}_frags_gc.bed  (from 03-filter_frags_gc.sh)
## file2 = ${sample}.bedpe         (original BEDPE)

awk '
    BEGIN {
        FS = OFS = "\t"
    }
    FNR == NR {
        # file1: *_frags_gc.bed
        # columns: 1=chr, 2=start, 3=end, 4=QNAME, 5=frag_len, 6=GC
        arr[$4] = $6  # key = QNAME, value = GC
        next
    }
    # file2: *.bedpe
    # assume column 7 = QNAME (read name)
    ($7 in arr) {
        # append GC content to BEDPE line
        print $0, arr[$7]
    }
' "${frags_dir}/${sample}_frags_gc.bed" "$bedpe" > "${output_dir}/${sample}_filtered.bedpe"

echo "Sample ${sample} processing complete. Output saved to ${output_dir}/${sample}_filtered.bedpe."

echo "All done for this array task."