#!/bin/bash
#SBATCH --job-name=VEP_GridssSV
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00
#SBATCH --array=1-8
#SBATCH --output=slurm_logs/VEP/vep_gridss_%A_%a.out
#SBATCH --error=slurm_logs/VEP/vep_gridss_%A_%a.err

set -euo pipefail

# ---------------- Modules ----------------
module purge
module load singularity/4.3.2
module load bcftools/1.22
module load htslib/1.22

# ---------------- Paths ----------------
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
gridss_base="${designated_path}/GRIDSS_012425"
out_base="${gridss_base}/GRIDSS_ANNOT_VEP"
mkdir -p "$out_base" slurm_logs/VEP

VEP_SIF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/VEP/vep.sif"   # container (VEP 115.x)
VEP_DIR="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/VEP"         # cache dir (VEP 113)
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"

# ---------------- VEP cache settings ----------------
cache_species="homo_sapiens"
cache_version="113"
cache_assembly="GRCh37"
cache_mount="/mnt/.vep"

# Allow large SVs
max_sv_size="500000000"

echo "[INFO] VEP_SIF:       $VEP_SIF"
echo "[INFO] VEP_DIR:       $VEP_DIR"
echo "[INFO] FASTA:         $reference"
echo "[INFO] cache_version: $cache_version"
echo "[INFO] max_sv_size:   $max_sv_size"

# ---------------- Resolve samples (FILTERED ONLY) ----------------
mapfile -t vcfs < <(
  find "$gridss_base" \
    -mindepth 2 -maxdepth 2 \
    -type f \
    -name "*.gridss.filtered.vcf.gz" \
    | sort
)

if [[ "${#vcfs[@]}" -eq 0 ]]; then
  echo "[ERROR] No filtered GRIDSS VCFs found under: $gridss_base"
  echo "        Expected pattern: */*.gridss.filtered.vcf.gz"
  exit 1
fi

# ---------------- Pick VCF for this array task ----------------
idx=$((SLURM_ARRAY_TASK_ID - 1))
if [[ "$idx" -lt 0 || "$idx" -ge "${#vcfs[@]}" ]]; then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID out of range (n=${#vcfs[@]})"
  exit 1
fi

in_vcf="${vcfs[$idx]}"
sample_dir="$(dirname "$in_vcf")"
sample="$(basename "$sample_dir")"

sample_outdir="${out_base}/${sample}"
mkdir -p "$sample_outdir"

# Keep base name consistent with your DELLY outputs
base_name="$(basename "$in_vcf" .vcf.gz)"

out_vcfgz="${sample_outdir}/${base_name}.vep${cache_version}.vcf.gz"
out_tsv="${sample_outdir}/${base_name}.vep${cache_version}.tsv"
log="${sample_outdir}/${sample}.gridss.vep.log"
warn="${sample_outdir}/${sample}_gridss.vep_warnings.txt"

echo "============================================================"
echo "[INFO] SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "[INFO] Sample: $sample"
echo "[INFO] Input:  $in_vcf"
echo "[INFO] Output: $out_vcfgz"
echo "============================================================"

# ---------------- Input validation ----------------
if [[ ! -s "$in_vcf" ]]; then
  echo "[ERROR] Missing/empty input VCF: $in_vcf"
  exit 1
fi

if [[ ! -s "${in_vcf}.tbi" ]]; then
  echo "[INFO] Indexing input VCF..."
  tabix -f -p vcf "$in_vcf" || true
fi

if [[ -f "$out_vcfgz" ]]; then
  echo "[SKIP] Already annotated: $out_vcfgz"
  exit 0
fi

# ---------------- Run VEP inside Singularity ----------------
# Notes for GRIDSS:
# - GRIDSS emits BND breakends; VEP handles BND and large SVs when --max_sv_size is set
# - --nearest symbol --distance helps tag nearby genes even for intergenic breakends
singularity exec \
  --writable-tmpfs \
  --no-home \
  --bind /ix1 \
  -B "$VEP_DIR:${cache_mount}" \
  "$VEP_SIF" \
  vep \
    --cache --offline \
    --dir_cache "${cache_mount}" \
    --species "${cache_species}" \
    --cache_version "${cache_version}" \
    --assembly "${cache_assembly}" \
    --fasta "$reference" \
    --max_sv_size "${max_sv_size}" \
    --vcf --compress_output bgzip --force_overwrite \
    --symbol --nearest symbol --distance 100000 \
    --no_stats \
    --fork 16 \
    --input_file "$in_vcf" \
    --output_file "$out_vcfgz" \
    > "$log" 2> "$warn"

# Index output VCF
tabix -f -p vcf "$out_vcfgz"

# Quick QC: confirm CSQ exists
if bcftools view -h "$out_vcfgz" | grep -q '##INFO=<ID=CSQ'; then
  echo "[SUCCESS] CSQ header present"
else
  echo "[WARNING] CSQ header missing – check: $log"
fi

# ---------------- TSV summary (GRIDSS BNDs) ----------------
bcftools +split-vep -s worst \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%Gene\t%SYMBOL\t%Consequence\t%IMPACT\t%DISTANCE\t%NEAREST\n' \
  -o "$out_tsv" -O v "$out_vcfgz" || true

# ---------------- Counts ----------------
echo "[INFO] input_variants:  $(bcftools view -H "$in_vcf" | wc -l)"
echo "[INFO] output_variants: $(bcftools view -H "$out_vcfgz" | wc -l)"
echo "[INFO] annotated_rows:  $(wc -l < "$out_tsv" 2>/dev/null || echo 0)"

echo "[DONE] $sample"
