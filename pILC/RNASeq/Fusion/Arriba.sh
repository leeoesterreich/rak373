#!/bin/bash
#SBATCH --job-name=ARRIBA
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 1-00:00
#SBATCH --array=0-30
#SBATCH --partition=htc
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/logs/%x.%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/logs/%x.%A_%a.err

module load arriba/2.5.0

SAMPLE_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/pILC_RNASeq_31samples.csv"

ALIGN_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/Alignment_arribaSTAR"
OUT_BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/ARRIBA"

FASTA="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/Human_genome_hg38_p7_chr_only.fna"
GTF="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/gencode.v47.basic.annotation.gtf"

BLACKLIST="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Arriba_db/database/blacklist_hg38_GRCh38_v2.5.0.tsv.gz"
KNOWN_FUSIONS="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Arriba_db/database/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz"
PROTEIN_DOMAINS="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Arriba_db/database/protein_domains_hg38_GRCh38_v2.5.0.gff3"

mkdir -p "${OUT_BASE}"

row=$((SLURM_ARRAY_TASK_ID + 2))
line=$(sed -n "${row}p" "${SAMPLE_SHEET}" || true)
[[ -n "${line}" ]] || { echo "No line for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

IFS=',' read -r SAMPLE R1 R2 <<< "${line}"
SAMPLE=$(echo "$SAMPLE" | tr -d '\r' | xargs)

BAM="${ALIGN_BASE}/${SAMPLE}/${SAMPLE}.Aligned.out.bam"
[[ -s "${BAM}" ]] || { echo "Missing BAM: ${BAM}"; exit 1; }

SAMPLE_OUT="${OUT_BASE}/${SAMPLE}"
mkdir -p "${SAMPLE_OUT}"

LOG="${SAMPLE_OUT}/${SAMPLE}.arriba.log"

echo "[$(date)] Arriba: ${SAMPLE}" | tee "${LOG}"
echo "BAM=${BAM}" | tee -a "${LOG}"

arriba \
  -x "${BAM}" \
  -o "${SAMPLE_OUT}/${SAMPLE}.arriba.fusions.tsv" \
  -O "${SAMPLE_OUT}/${SAMPLE}.arriba.fusions.discarded.tsv" \
  -a "${FASTA}" \
  -g "${GTF}" \
  -b "${BLACKLIST}" \
  -k "${KNOWN_FUSIONS}" \
  -t "${PROTEIN_DOMAINS}" \
  >> "${LOG}" 2>&1

echo "[$(date)] Done Arriba: ${SAMPLE}" | tee -a "${LOG}"
