#!/bin/bash

# This script is used to filter Delly output vcf file.
# High-confidence filter:
            #  - FILTER must be PASS
            #  - at least 1 split-read (SR>=1)
            #  - PE + SR >= 8 total support
            #  - mapping quality for paired-ends and split-reads >= 30

#SBATCH --job-name=Delly_filter
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 01:00:00
#SBATCH --array=1-4
#SBATCH --output=slurm_log/Delly/Delly_filter_%A_%a.out
#SBATCH --error=slurm_log/Delly/Delly_filter_%A_%a.err

module load bcftools/1.22

base_dir="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/WGS/Delly"

# Make sure log dir exists
mkdir -p slurm_log/Delly

# Get a list of sample dirs
sample_dirs=("${base_dir}"/*)

# Pick the one corresponding to this array task
idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
sample_dir="${sample_dirs[$idx]}"

[ -d "$sample_dir" ] || { echo "[ERROR] No sample dir for index $idx"; exit 1; }

sample_name=$(basename "$sample_dir")
echo "[INFO] Processing sample: $sample_name"

for sv in BND DEL DUP INS INV; do
    vcf="${sample_dir}/${sample_name}_${sv}.vcf"
    out="${sample_dir}/${sample_name}_filtered_${sv}.vcf"

    if [ -f "$vcf" ]; then
        echo "  -> Filtering $sv VCF: $(basename "$vcf")"

        bcftools view \
            -i 'FILTER="PASS" && INFO/SR>=1 && (INFO/PE + INFO/SR) >= 8 && INFO/MAPQ>=30 && INFO/SRMAPQ>=30' \
            -O v -o "$out" "$vcf"

        echo "     Saved: $(basename "$out")"
    else
        echo "  -> VCF not found for $sv: $vcf"
    fi
done

echo "[INFO] All done for $sample_name."
