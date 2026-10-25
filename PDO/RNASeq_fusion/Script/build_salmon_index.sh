#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --clusters=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/salmon_index_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/salmon_index_%j.err
# Build a Salmon index from a transcriptome (cDNA) FASTA. Needed once per species/reference;
# check_strand_salmon.sh uses it to call library strandedness.
#
# Usage: sbatch Script/build_salmon_index.sh [CDNA_GZ] [INDEX]
#   CDNA_GZ  gzipped cDNA FASTA   (default: human GRCh38, Ensembl 114)
#   INDEX    output index folder  (default: salmon_index_GRCh38_114 in the project folder)
# Another species: give its cDNA FASTA and a NEW index folder name, e.g. salmon_index_mouse
# Rebuild an existing index: FORCE=1 sbatch Script/build_salmon_index.sh

set -u
module load salmon/1.10.3

BASE=${BASE:-/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion}
CDNA_GZ=${1:-/ix1/alee/LO_LAB/General/References/homo_sapiens_sapiens/GRCh38_114/Homo_sapiens.GRCh38.cdna.all.fa.gz}
INDEX=${2:-${BASE}/salmon_index_GRCh38_114}
THREADS=${SLURM_CPUS_PER_TASK:-16}
TMP_INDEX=${INDEX}.building
TMP_FA=${INDEX}.cdna.fa

[[ -f ${CDNA_GZ} ]] || { echo "ERROR: FASTA not found: ${CDNA_GZ}" >&2; exit 1; }
mkdir -p "${BASE}/logs" "$(dirname ${INDEX})"

if [[ -f ${INDEX}/info.json && "${FORCE:-0}" != "1" ]]; then
    echo "Index already exists: ${INDEX} (set FORCE=1 to rebuild). Nothing to do."
    exit 0
fi

echo "Building ${INDEX} from ${CDNA_GZ}"
rm -rf "${TMP_INDEX}" "${TMP_FA}"
zcat ${CDNA_GZ} > ${TMP_FA}
# -k 31 = k-mer size, suited to reads of about 75 bp or longer
salmon index -t ${TMP_FA} -i ${TMP_INDEX} -k 31 -p ${THREADS} \
    && [[ -f ${TMP_INDEX}/info.json ]] \
    || { echo "ERROR: index build failed" >&2; rm -rf "${TMP_INDEX}" "${TMP_FA}"; exit 1; }

rm -rf "${INDEX}"; mv "${TMP_INDEX}" "${INDEX}"     # only a finished index gets the final name
rm -f "${TMP_FA}"
echo "Done: ${INDEX}"; du -sh "${INDEX}"
