#!/bin/bash

# This script is used to create PON for somatic mutation calling (template provided by Brent)

#SBATCH --job-name=gatk_pon
#SBATCH --output=slurm_log/mutect2/pon/gatk_pon_%j.out
#SBATCH --error=slurm_log/mutect2/pon/gatk_pon_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32

module load singularity/4.3.2

REF="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
INTERVALS="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/wgs_calling_regions.hg38.interval_list"
WORKSPACE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/pilc_WES_pon_db_042026_rak373"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_PON="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/pon_rak373_${TIMESTAMP}.vcf.gz"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif"
BIND="--bind /ix1"

VCF_LIST="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/5_Mutect2/pon_vcfs.txt"

# Read all lines (with --variant) into an array for direct passing
VCF_ARGS=()
while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    VCF_ARGS+=("$line")
done < "$VCF_LIST"

echo "Starting GenomicsDBImport at $(date)"

singularity exec $BIND $SIF gatk GenomicsDBImport \
  -R "$REF" \
  -L "$INTERVALS" \
  --merge-input-intervals \
  --genomicsdb-shared-posixfs-optimizations true \
  --genomicsdb-workspace-path "$WORKSPACE" \
  ${VCF_ARGS[@]}

echo "Finished GenomicsDBImport at $(date)"
echo "Starting CreateSomaticPanelOfNormals at $(date)"

singularity exec $BIND $SIF gatk CreateSomaticPanelOfNormals \
  -R "$REF" \
  -V "gendb://${WORKSPACE}" \
  -O "$OUTPUT_PON"

echo "Finished CreateSomaticPanelOfNormals at $(date)"
echo "PoN VCF output: $OUTPUT_PON"