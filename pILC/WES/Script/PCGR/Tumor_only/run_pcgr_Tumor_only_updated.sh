#!/bin/bash

# This script is used to run PCGR
#SBATCH --job-name=pcgr
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-31
#SBATCH --output=slurm_log/pcgr/pcgr_%A_%a.out
#SBATCH --error=slurm_log/pcgr/pcgr_%A_%a.err

# Load Singularity
module load singularity/4.3.2

export VEP_VERSION="113"
export GENOME="GRCh38"
export CACHE="homo_sapiens_vep_${VEP_VERSION}_${GENOME}.tar.gz"

# echo Singularity version
singularity --version

# ====== Configurable paths ======
PCGR_DIR="/ix1/alee/LO_LAB/Personal/Brent_Schlegel/tools/PCGR"
VEP_DIR="/ix1/alee/LO_LAB/Personal/Brent_Schlegel/tools/VEP"
SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/pcgr_2.2.3.sif"

# Input data
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Tumor_only"
VCF_DIR="$INPUT_DIR/vcf_pcgr"
TPM_DIR="$INPUT_DIR/tpm_pcgr"
CNA_DIR="$INPUT_DIR/cna_segs"
INDEX_FILE="$INPUT_DIR/With_Tumor_only_index.csv"

# Output
OUT_DIR="$INPUT_DIR/pcgr_outputs"
mkdir -pv "$OUT_DIR"

# ====== Get sample info from CSV ======
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$INDEX_FILE")
IFS=',' read -r vcf_file tpm_file cna_file <<< "$LINE"
sample_id=$(basename "$vcf_file" .taf.vcf.gz)

echo "Running PCGR on Sample: $sample_id"

# ====== Validate files ======
vcf_path="$VCF_DIR/$vcf_file"
[ ! -f "$vcf_path" ] && echo "Missing VCF: $vcf_path" && exit 1

# ====== RNA Expression Flag ======
rna_expr_flag=""
if [ -n "$tpm_file" ] && [ -f "$TPM_DIR/$tpm_file" ]; then
  rna_expr_flag="--input_rna_expression /mnt/tpm/$tpm_file --expression_sim --expression_sim_db tcga,depmap,treehouse"
  echo "[INFO] RNA file exists: $TPM_DIR/$tpm_file"
  echo "[INFO] RNA expression flag: $rna_expr_flag"
else
  echo "[INFO] RNA expression file missing or not provided: $TPM_DIR/$tpm_file"
fi

# ====== CNA Segments Flag ======
cna_flag=""
if [ -n "$cna_file" ] && [ -f "$CNA_DIR/$cna_file" ]; then
  cna_flag="--input_cna /mnt/cna/$cna_file"
  echo "[INFO] CNA file exists: $CNA_DIR/$cna_file"
  echo "[INFO] CNA flag: $cna_flag"
else
  echo "[INFO] CNA file missing or not provided: $CNA_DIR/$cna_file"
fi

# Create sample output dir
sample_out_dir="$OUT_DIR/$sample_id"
mkdir -p "$sample_out_dir"

# ---- Quarto + R for report generation ----
QUARTO_ROOT="/ix1/alee/LO_LAB/Personal/Rahul/tools/quarto/quarto-1.5.57"
QUARTO_BIN="$QUARTO_ROOT/bin/quarto"

if [[ ! -x "$QUARTO_BIN" ]]; then
  echo "ERROR: Quarto not found/executable at: $QUARTO_BIN"
  exit 1
fi
echo "Quarto on host: $($QUARTO_BIN --version)"

# Ensure PCGR can see Quarto + R (pcgrr env) inside container
export SINGULARITYENV_PREPEND_PATH="/mnt/quarto/bin:/opt/mambaforge/envs/pcgrr/bin:/opt/mambaforge/envs/pcgr/bin"

# Quarto writable dirs
export SINGULARITYENV_XDG_CACHE_HOME="/tmp/quarto_cache_home"
export SINGULARITYENV_XDG_CONFIG_HOME="/tmp/quarto_config_home"
export SINGULARITYENV_XDG_DATA_HOME="/tmp/quarto_data_home"
export SINGULARITYENV_QUARTO_CACHE_DIR="/tmp/quarto_cache"
export SINGULARITYENV_TMPDIR="/tmp"
export SINGULARITYENV_LANG="C.UTF-8"
export SINGULARITYENV_LC_ALL="C.UTF-8"

# ====== Run PCGR via Singularity ======
singularity exec \
  --writable-tmpfs \
  --home /tmp \
  --env "XDG_CACHE_HOME=/tmp/quarto_cache_home" \
  --bind "/ix1" \
  -B "$VEP_DIR:/mnt/.vep" \
  -B "$PCGR_DIR/bundle:/mnt/bundle" \
  -B "$QUARTO_ROOT:/mnt/quarto" \
  -B "$VCF_DIR:/mnt/pcgr_inputs" \
  -B "$OUT_DIR:/mnt/pcgr_outputs" \
  -B "$TPM_DIR:/mnt/tpm" \
  -B "$CNA_DIR:/mnt/cna" \
  "$SIF" \
  pcgr \
    --input_vcf "/mnt/pcgr_inputs/$vcf_file" \
    --vep_dir "/mnt/.vep" \
    --refdata_dir "/mnt/bundle" \
    --output_dir "/mnt/pcgr_outputs/$sample_id" \
    --genome_assembly "grch38" \
    --sample_id "$sample_id" \
    --assay "WES" \
    --tumor_only \
    --tumor_site 6 \
    --tumor_dp_tag "TDP" \
    --tumor_af_tag "TAF" \
    --estimate_signatures \
    --estimate_msi \
    --estimate_tmb \
    --vcf2maf \
    --tumor_dp_min 4 \
    --tumor_af_min 0.01 \
    --force_overwrite \
    $rna_expr_flag \
    $cna_flag
