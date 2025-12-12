#!/bin/bash

# This script is used to filter Manta output VCF files.
# Filtering criteria:
    #  1) FILTER == PASS
    #  2) SVTYPE in {DEL,DUP,INV,BND,INS}
    #  3) Evidence: PR[1] + SR[1] >= 5 alt-support reads
    #  4) Exclude blacklist regions

#SBATCH --job-name=manta_filter
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --array=1-4
#SBATCH --output=slurm_log/Manta/manta_run_%j.out
#SBATCH --error=slurm_log/Manta/manta_run_%j.err

set -euo pipefail

# ============================================================
# User settings
# ============================================================
# Base path for your WGS project
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/WGS"

# Where Manta results live
manta_dir="${designated_path}/Manta"

# Blacklist BED (same as used for GRIDSS)
blacklist="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Blacklist_file/hg38_unified_blacklist.bed"

# Load bcftools
module purge
module load bcftools/1.22

# ============================================================
# Loop over all Manta candidateSV VCFs
# Expected structure: ${manta_dir}/${SAMPLE}/${SAMPLE}_manta_candidateSV.vcf.gz
# ============================================================
for vcf in "${manta_dir}"/*/*_manta_candidateSV.vcf.gz; do
    [ -e "$vcf" ] || { echo "No candidateSV VCFs found in ${manta_dir}"; exit 0; }

    sample_dir=$(dirname "$vcf")
    sample_file=$(basename "$vcf")
    sample="${sample_file%_manta_candidateSV.vcf.gz}"

    out_vcf="${sample_dir}/${sample}_manta_filteredSV.vcf.gz"

    echo "============================================================"
    echo "[INFO] Processing sample: $sample"
    echo "[INFO] Input VCF:  $vcf"
    echo "[INFO] Output VCF: $out_vcf"
    echo "============================================================"

    if [ -f "$out_vcf" ]; then
        echo "[SKIP] Filtered VCF already exists for $sample"
        continue
    fi

    bcftools view \
        -T ^"$blacklist" \
        -i 'FILTER="PASS" && INFO/SVTYPE ~ "DEL|DUP|INV|BND|INS" && (FORMAT/PR[1] + FORMAT/SR[1]) >= 5' \
        -Oz \
        -o "$out_vcf" \
        "$vcf"

    bcftools index -t "$out_vcf"

    echo "[SUCCESS] Finished filtering for $sample"
    echo ""
done

echo "[INFO] All done."
