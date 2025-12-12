#!/bin/bash

# This script is used to filter GRIDSS output VCF files.
# High-confidence filter:
        #  - FILTER must be PASS
        #  - at least 1 split/soft-clipped read in assemblies (ASSR>=1)
        #  - ASRP + ASSR >= 8 total reads in assemblies

#SBATCH --job-name=Gridss_filter
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --output=slurm_log/Gridss/Gridss_filter_%j.out
#SBATCH --error=slurm_log/Gridss/Gridss_filter_%j.err

# Load required modules
module load bcftools/1.22

# Base directory where all GRIDSS sample folders live
base_dir="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/WGS/Gridss/"

# Make sure log dir exists
mkdir -p slurm_log/Gridss

# Loop over each sample directory
for sample_dir in "${base_dir}"/*; do
    # Skip if not a directory
    [ -d "$sample_dir" ] || continue

    sample_name=$(basename "$sample_dir")
    echo "[INFO] Processing sample: $sample_name"

    vcf="${sample_dir}/${sample_name}_gridss.vcf"
    out="${sample_dir}/${sample_name}_filtered_gridss.vcf"

    if [ -f "$vcf" ]; then
        echo "  -> Filtering GRIDSS VCF: $(basename "$vcf")"
        bcftools view \
            -i 'FILTER="PASS" && INFO/ASSR>=1 && (INFO/ASRP + INFO/ASSR) >= 8' \
            -O v -o "$out" "$vcf"

        echo "     Saved: $(basename "$out")"
    else
        echo "  -> VCF not found: $vcf"
    fi

    echo
done

echo "[INFO] All done."
