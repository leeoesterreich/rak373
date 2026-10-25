#!/bin/bash
#SBATCH --job-name=strand_check
#SBATCH --clusters=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/strand_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/strand_%j.err
# Usage: sbatch Script/check_strand_salmon.sh <FASTQ_DIR> <OUT_NAME> [R1_SUFFIX]
#   FASTQ_DIR  folder with paired FASTQ files
#   OUT_NAME   result folder name inside RNASeq_fusion (one per batch)
#   R1_SUFFIX  optional, default _R1_001.fastq.gz
# Optional: NREADS, CDNA_GZ, INDEX, CLEAN_NAMES

set -u
module load salmon/1.10.3

FASTQ_DIR=${1:-/ix1/alee/LO_LAB/General/Lab_Data/20240315_Organoid_RNAseq/Organoid_mRNA_03_20_2024_Children_Hosp_RAW}
OUT_NAME=${2:-strand_check}
R1_SUFFIX=${3:-_R1_001.fastq.gz}

BASE=${BASE:-/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion}
CDNA_GZ=${CDNA_GZ:-/ix1/alee/LO_LAB/General/References/homo_sapiens_sapiens/GRCh38_114/Homo_sapiens.GRCh38.cdna.all.fa.gz}
INDEX=${INDEX:-${BASE}/salmon_index_GRCh38_114}
NREADS=${NREADS:-1000000}
CLEAN_NAMES=${CLEAN_NAMES:-1}
THREADS=${SLURM_CPUS_PER_TASK:-16}
OUT=${BASE}/${OUT_NAME}

R2_SUFFIX=${R1_SUFFIX/R1/R2}
[[ "${R2_SUFFIX}" == "${R1_SUFFIX}" ]] && R2_SUFFIX=${R1_SUFFIX/_1./_2.}
[[ "${R2_SUFFIX}" == "${R1_SUFFIX}" ]] && { echo "ERROR: cannot derive R2 suffix" >&2; exit 1; }

[[ -d "${FASTQ_DIR}" ]] || { echo "ERROR: folder not found: ${FASTQ_DIR}" >&2; exit 1; }
FILES=$(ls ${FASTQ_DIR}/*${R1_SUFFIX} 2>/dev/null || true)
[[ -n "${FILES}" ]] || { echo "ERROR: no files ending ${R1_SUFFIX} in ${FASTQ_DIR}" >&2; exit 1; }
N_TOTAL=$(echo "${FILES}" | wc -l)
mkdir -p "${OUT}" "${BASE}/logs"

make_name() {   # IPM-BO-287_S1 -> IPM_BO_287
    if [[ "${CLEAN_NAMES}" == "1" ]]; then echo "$1" | sed -E 's/_S[0-9]+$//; s/-/_/g'; else echo "$1"; fi
}
NAMES=$(for R1 in ${FILES}; do make_name "$(basename ${R1} ${R1_SUFFIX})"; done)
DUP=$(echo "${NAMES}" | sort | uniq -d)
[[ -z "${DUP}" ]] || { echo "ERROR: duplicate sample names: ${DUP}" >&2; exit 1; }
echo "Samples found: ${N_TOTAL} | Output: ${OUT} | Reads/sample: ${NREADS}"

# Salmon index: built once, reused by every batch
if [[ ! -f ${INDEX}/info.json ]]; then
    rm -rf "${INDEX}"
    zcat ${CDNA_GZ} > ${INDEX}.cdna.fa
    salmon index -t ${INDEX}.cdna.fa -i ${INDEX} -k 31 -p ${THREADS}
    rm -f ${INDEX}.cdna.fa
fi

I=0; N_FAIL=0
for R1 in ${FILES}; do
    I=$((I+1))
    R2=${R1/%${R1_SUFFIX}/${R2_SUFFIX}}
    S=$(make_name "$(basename ${R1} ${R1_SUFFIX})")
    if [[ -f ${OUT}/${S}/lib_format_counts.json ]]; then echo "[${I}/${N_TOTAL}] ${S}: done, skipping"; continue; fi
    if [[ ! -f ${R2} ]]; then echo "[${I}/${N_TOTAL}] ${S}: R2 missing" >&2; N_FAIL=$((N_FAIL+1)); continue; fi
    echo "[${I}/${N_TOTAL}] ${S}: running"
    salmon quant -i ${INDEX} -l A -p ${THREADS} \
        -1 <(zcat ${R1} | head -n $((NREADS*4))) \
        -2 <(zcat ${R2} | head -n $((NREADS*4))) \
        -o ${OUT}/${S} > ${OUT}/${S}.salmon.log 2>&1 \
        || { echo "[${I}/${N_TOTAL}] ${S}: FAILED, see ${OUT}/${S}.salmon.log" >&2; N_FAIL=$((N_FAIL+1)); }
done

${BASE}/Script/summarise_strand.sh ${OUT}
N_DONE=$(ls ${OUT}/*/lib_format_counts.json 2>/dev/null | wc -l)
echo "Finished: ${N_DONE} of ${N_TOTAL} samples have a result, ${N_FAIL} problems"
[[ "${N_FAIL}" -eq 0 ]] || exit 1
