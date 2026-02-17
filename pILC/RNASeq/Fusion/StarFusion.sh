#!/usr/bin/env bash
#SBATCH --job-name=StarFusion_pILC
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -t 2-00:00
#SBATCH --array=1-31
#SBATCH --partition=htc
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/slurm_log/%x.%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/slurm_log/%x.%A_%a.err

set -euo pipefail

module load star-fusion/1.15.1

# -----------------------------
# Paths
# -----------------------------
SAMPLE_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/pILC_RNASeq_31samples.csv"
CTAT_LIB="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/StarFusion_db/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir"
OUT_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/StarFusion"
LOG_DIR="${OUT_BASE}/slurm_log"

mkdir -p "${OUT_BASE}" "${LOG_DIR}"

# -----------------------------
# Pull sample for this array task
#   - Skip header (NR>1)
#   - Use SLURM_ARRAY_TASK_ID as the (1-based) data row index
# -----------------------------
LINE="$(awk -F',' -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i+1 {print; exit}' "${SAMPLE_SHEET}")"

if [[ -z "${LINE}" ]]; then
  echo "ERROR: No line found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${SAMPLE_SHEET}"
  exit 1
fi

SAMPLE="$(echo "${LINE}" | awk -F',' '{print $1}')"
LEFT="$(echo  "${LINE}" | awk -F',' '{print $2}')"
RIGHT="$(echo "${LINE}" | awk -F',' '{print $3}')"

# Basic cleanup (just in case there are stray CRLFs)
SAMPLE="${SAMPLE//$'\r'/}"
LEFT="${LEFT//$'\r'/}"
RIGHT="${RIGHT//$'\r'/}"

OUTDIR="${OUT_BASE}/${SAMPLE}"
TMPDIR="${OUTDIR}/_tmp"
mkdir -p "${OUTDIR}" "${TMPDIR}"

echo "HOST: $(hostname)"
echo "START: $(date)"
echo "TASK: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "SAMPLE: ${SAMPLE}"
echo "LEFT:   ${LEFT}"
echo "RIGHT:  ${RIGHT}"
echo "OUTDIR:  ${OUTDIR}"
echo "CTAT:    ${CTAT_LIB}"
echo "CPUs:    ${SLURM_CPUS_PER_TASK}"

# -----------------------------
# Sanity checks
# -----------------------------
[[ -s "${LEFT}"  ]] || { echo "ERROR: LEFT fq missing/empty: ${LEFT}"; exit 2; }
[[ -s "${RIGHT}" ]] || { echo "ERROR: RIGHT fq missing/empty: ${RIGHT}"; exit 2; }

[[ -s "${CTAT_LIB}/ref_genome.fa"  ]] || { echo "ERROR: CTAT ref_genome.fa not found"; exit 3; }
[[ -s "${CTAT_LIB}/ref_annot.gtf"  ]] || { echo "ERROR: CTAT ref_annot.gtf not found"; exit 3; }
[[ -d "${CTAT_LIB}/ref_genome.fa.star.idx" || -d "${CTAT_LIB}/star.idx" ]] || {
  echo "ERROR: CTAT STAR index dir not found (expected ref_genome.fa.star.idx or star.idx)"
  exit 3
}

# -----------------------------
# Run STAR-Fusion
# -----------------------------
STAR-Fusion \
  --genome_lib_dir "${CTAT_LIB}" \
  --left_fq "${LEFT}" \
  --right_fq "${RIGHT}" \
  --CPU "${SLURM_CPUS_PER_TASK}" \
  --output_dir "${OUTDIR}" \
  --tmpdir "${TMPDIR}"

echo "DONE: $(date)"

# Quick summary
PRED="${OUTDIR}/star-fusion.fusion_predictions.tsv"
if [[ -s "${PRED}" ]]; then
  echo "Predictions: ${PRED}"
  echo "Num lines (incl header): $(wc -l < "${PRED}")"
else
  echo "WARNING: Predictions TSV not found or empty at ${PRED}"
fi
