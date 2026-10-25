#!/bin/bash
#SBATCH --job-name=copy_res
#SBATCH --clusters=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/copy_res_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/copy_res_%j.err
set -euo pipefail

SRC=/vast/alee/Rahul/PDO/RNAFusion_4.1.3/results
DEST=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/20240315_Organoid_RNAseq/results
mkdir -p "${DEST}"

# Safety check: stop unless at least 300 GB is free on /ix1
AVAIL_GB=$(df -BG --output=avail "${DEST}" | tail -n 1 | tr -dc '0-9')
echo "free on ix1: ${AVAIL_GB} GB"
if [[ "${AVAIL_GB}" -lt 300 ]]; then echo "Not enough free space, stopping"; exit 1; fi

# Folders copied completely
FULL="arriba arriba_visualisation fusioncatcher fusioninspector fusionreport starfusion ctatsplicing stringtie salmon agat multiqc fastqc fastqc_for_fastp pipeline_info vcf"
# Folders where only small QC files are kept (no BAM / CRAM / FASTQ files)
QC="fastp star samtools picard"

copy_full() { rsync -a "${SRC}/$1/" "${DEST}/$1/"; }
copy_qc()   { rsync -a --exclude='*.bam' --exclude='*.bai' --exclude='*.cram' \
                 --exclude='*.crai' --exclude='*.fastq.gz' --exclude='*.fq.gz' \
                 "${SRC}/$1/" "${DEST}/$1/"; }
export -f copy_full copy_qc
export SRC DEST

# 8 folders at a time
echo ${FULL} | tr ' ' '\n' | xargs -P 8 -I{} bash -c 'copy_full {}'
echo ${QC}   | tr ' ' '\n' | xargs -P 8 -I{} bash -c 'copy_qc {}'

echo "copy finished"
du -sh "${DEST}"

# ---- Verify: rsync dry run, lists files missing or different in DEST ----
echo "Verifying ..."
TOTAL=0
for f in ${FULL}; do
    N=$(rsync -an --itemize-changes "${SRC}/${f}/" "${DEST}/${f}/" | grep -cE '^(>f|cd|cL|cS)' || true)
    echo "${f}  ${N}"; TOTAL=$((TOTAL + N))
done
for f in ${QC}; do
    N=$(rsync -an --itemize-changes --exclude='*.bam' --exclude='*.bai' --exclude='*.cram' \
          --exclude='*.crai' --exclude='*.fastq.gz' --exclude='*.fq.gz' \
          "${SRC}/${f}/" "${DEST}/${f}/" | grep -cE '^(>f|cd|cL|cS)' || true)
    echo "${f}  ${N}"; TOTAL=$((TOTAL + N))
done
echo "TOTAL differences: ${TOTAL}"
if [[ ${TOTAL} -eq 0 ]]; then echo "VERIFIED: copy is complete"; else echo "NOT VERIFIED: do NOT delete anything on /vast"; exit 1; fi
