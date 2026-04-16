#!/bin/bash

#SBATCH --job-name=reformat_vcf
#SBATCH --output=slurm_log/reformat_vcf/reformat_vcf_%A_%a.out
#SBATCH --error=slurm_log/reformat_vcf/reformat_vcf_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# This script reformats VCF files for PCGR (adds TAF/TDP/NAF/NDP fields) 
# Based on pysam parsing and bcftools normalization #

# ======================================
# Load modules
# ======================================
module load bcftools/1.16
module load singularity/4.3.2

# ======================================
# Define paths
# ======================================
VCF_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/PCGR"
VCF_LIST="$VCF_DIR/filtered_vcf_file_list.txt"
REF_FASTA="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta"
OUT_DIR="$VCF_DIR/vcf_pcgr"
PYTHON_SCRIPT="$VCF_DIR/reformat_vcf.py"
SINGULARITY_IMAGE="/ix1/alee/LO_LAB/Personal/Brent_Schlegel/envs/pysam-xcmds_latest.sif"

# ======================================
# Setup output/log folders
# ======================================
mkdir -p "$OUT_DIR"
mkdir -p "slurm_log/reformat_vcf"

# ======================================
# Read VCF entry for this SLURM task
# ======================================
VCF_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$VCF_LIST" | tr -d '\r')
VCF_BASENAME=$(basename "$VCF_FILE" .vcf.gz)
SAMPLE_NAME=$(echo "$VCF_BASENAME" | cut -d'.' -f1)

INPUT_VCF="$VCF_DIR/$VCF_FILE"
OUTPUT_VCF="$OUT_DIR/${SAMPLE_NAME}.taf.vcf.gz"
WORK_DIR="$OUT_DIR/tmp_${SLURM_ARRAY_TASK_ID}_${SAMPLE_NAME}"

mkdir -pv "$WORK_DIR"
cd "$WORK_DIR"

echo "=========================================================="
echo "Processing sample: $SAMPLE_NAME"
echo "Resolved input VCF: $INPUT_VCF"
echo "Expected output VCF: $OUTPUT_VCF"
echo "Python script: $PYTHON_SCRIPT"
echo "Work directory: $WORK_DIR"
echo "=========================================================="

# ======================================
# Validate input
# ======================================
if [[ ! -f "$INPUT_VCF" ]]; then
  echo "ERROR: Input VCF not found: $INPUT_VCF"
  echo "Contents of $VCF_LIST:"
  cat -A "$VCF_LIST"
  exit 1
fi

# ======================================
# Decompress if necessary
# ======================================
if [[ "$INPUT_VCF" == *.gz ]]; then
  echo "Decompressing VCF..."
  gunzip -c "$INPUT_VCF" > input.vcf
  INPUT_VCF_FOR_SCRIPT="input.vcf"
else
  INPUT_VCF_FOR_SCRIPT="$INPUT_VCF"
fi

# ======================================
# Run Python script inside container
# ======================================
echo "Running reformat_vcf.py via Singularity..."

singularity exec --bind /ix1 "$SINGULARITY_IMAGE" python "$PYTHON_SCRIPT" \
  -vcf_file "$INPUT_VCF_FOR_SCRIPT" \
  -out "${OUTPUT_VCF%.gz}" \
  -reference "$REF_FASTA" \
  -tumor_sample "$SAMPLE_NAME"

# ======================================
# Post-run: check if output was created
# ======================================
echo "Files in work dir:"
ls -lh

if [[ -f "$OUTPUT_VCF" ]]; then
  echo "SUCCESS: Output created: $OUTPUT_VCF"
elif [[ -f "${OUTPUT_VCF%.gz}" ]]; then
  echo "Uncompressed output exists — bgzip may have failed"
else
  echo "ERROR: No output VCF created for $SAMPLE_NAME"
  exit 1
fi

# ======================================
# Cleanup
# ======================================
cd "$OUT_DIR"
rm -rf "$WORK_DIR"
echo "Cleaned up $WORK_DIR"
echo "Job complete for $SAMPLE_NAME"
