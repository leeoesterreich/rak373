#!/bin/bash
#
# Run LUMPY Express using BWA-MEM + samblaster,
# then extract/sort splitters + discordants and call SVs with lumpyexpress.

#SBATCH --job-name=LUMPY
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G
#SBATCH --array=0-7
#SBATCH --output=slurm_log/LUMPY/LUMPY_%A_%a.out
#SBATCH --error=slurm_log/LUMPY/LUMPY_%A_%a.err

set -euo pipefail

module purge
module load bwa/0.7.17
module load samtools/1.22.1
module load singularity/4.3.2
module load lumpy-sv/0.3.1

# =========================
# Paths
# =========================
FASTQ_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20250124_Rahul_Jagmohan_EnRich/Illumina_Enrichment"
OUTROOT="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Lumpy/Lumpy_012425"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"
SIF="/software/rhel9/manual/install/lumpy-sv/lumpy-sv-0.3.1.sif"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "${OUTROOT}" slurm_log/LUMPY

# =========================
# Build sample list (TP*_R1 only), map array index -> sample
# =========================
mapfile -t R1_LIST < <(ls -1 "${FASTQ_DIR}"/TP*_R1_001.fastq.gz 2>/dev/null | sort)

if [[ "${#R1_LIST[@]}" -eq 0 ]]; then
  echo "ERROR: No TP*_R1_001.fastq.gz files found in ${FASTQ_DIR}" >&2
  exit 1
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${#R1_LIST[@]}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (N=${#R1_LIST[@]})" >&2
  exit 1
fi

R1="${R1_LIST[${SLURM_ARRAY_TASK_ID}]}"
R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"

# Sample name = prefix before _S{number}_R1_001.fastq.gz
# e.g., TP24_M292_S2_R1_001.fastq.gz -> TP24_M292
SAMPLE="$(basename "${R1}" | sed -E 's/_S[0-9]+_R1_001\.fastq\.gz$//')"

echo "[`date`] Task ${SLURM_ARRAY_TASK_ID}/${#R1_LIST[@]} -> SAMPLE=${SAMPLE}"
echo "[`date`] R1=${R1}"
echo "[`date`] R2=${R2}"

# =========================
# Sanity checks
# =========================
[[ -s "${reference}" ]] || { echo "ERROR: reference missing: ${reference}" >&2; exit 1; }
[[ -s "${R1}" ]] || { echo "ERROR: R1 missing: ${R1}" >&2; exit 1; }
[[ -s "${R2}" ]] || { echo "ERROR: R2 missing: ${R2}" >&2; exit 1; }
[[ -s "${SIF}" ]] || { echo "ERROR: SIF missing: ${SIF}" >&2; exit 1; }

# =========================
# Per-sample output paths
# =========================
OUTDIR="${OUTROOT}/${SAMPLE}"
mkdir -p "${OUTDIR}"

BAM_UNSORT="${OUTDIR}/${SAMPLE}.unsorted.bam"
BAM="${OUTDIR}/${SAMPLE}.bam"

DISC_UNSORT="${OUTDIR}/${SAMPLE}.discordants.unsorted.bam"
SPLIT_UNSORT="${OUTDIR}/${SAMPLE}.splitters.unsorted.bam"

DISC_SORT="${OUTDIR}/${SAMPLE}.discordants.bam"
SPLIT_SORT="${OUTDIR}/${SAMPLE}.splitters.bam"

VCF="${OUTDIR}/${SAMPLE}.lumpy.vcf"

# Optional: skip if already finished
if [[ -s "${VCF}" ]]; then
  echo "[`date`] VCF exists, skipping: ${VCF}"
  exit 0
fi

# =========================
# Step 1: Align -> UNSORTED BAM
# =========================
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:EnRICH\tPL:ILLUMINA"

echo "[`date`] Aligning (bwa mem) | samblaster (SIF) -> ${BAM_UNSORT}"
bwa mem -t "${THREADS}" -R "${RG}" "${reference}" "${R1}" "${R2}" \
  | singularity exec --no-home --writable-tmpfs \
      --bind /ix1,/ihome,/scratch,/software \
      "${SIF}" samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
  | samtools view -@ "${THREADS}" -Sb - \
  > "${BAM_UNSORT}"

# =========================
# Step 1b: Coordinate sort + index
# =========================
echo "[`date`] Sorting BAM -> ${BAM}"
samtools sort -@ "${THREADS}" -o "${BAM}" "${BAM_UNSORT}"
samtools index -@ "${THREADS}" "${BAM}"
rm -f "${BAM_UNSORT}"

# =========================
# Step 2: Extract discordants
# =========================
echo "[`date`] Extracting discordants -> ${DISC_UNSORT}"
samtools view -@ "${THREADS}" -b -F 1294 "${BAM}" > "${DISC_UNSORT}"

# =========================
# Step 3: Extract splitters (extractSplitReads_BwaMem inside SIF)
# =========================
echo "[`date`] Extracting splitters -> ${SPLIT_UNSORT}"
samtools view -@ "${THREADS}" -h "${BAM}" \
  | singularity exec --no-home --writable-tmpfs \
      --bind /ix1,/ihome,/scratch,/software \
      "${SIF}" extractSplitReads_BwaMem -i stdin \
  | samtools view -@ "${THREADS}" -Sb - \
  > "${SPLIT_UNSORT}"

# =========================
# Step 4: Sort + index discordants/splitters
# =========================
echo "[`date`] Sorting discordants -> ${DISC_SORT}"
samtools sort -@ "${THREADS}" -o "${DISC_SORT}" "${DISC_UNSORT}"
samtools index -@ "${THREADS}" "${DISC_SORT}"
rm -f "${DISC_UNSORT}"

echo "[`date`] Sorting splitters -> ${SPLIT_SORT}"
samtools sort -@ "${THREADS}" -o "${SPLIT_SORT}" "${SPLIT_UNSORT}"
samtools index -@ "${THREADS}" "${SPLIT_SORT}"
rm -f "${SPLIT_UNSORT}"

# =========================
# Step 5: LUMPY Express
# =========================
echo "[`date`] Running lumpyexpress -> ${VCF}"
singularity exec --no-home --writable-tmpfs \
  --bind /ix1,/ihome,/scratch,/software \
  "${SIF}" lumpyexpress \
    -B "${BAM}" \
    -S "${SPLIT_SORT}" \
    -D "${DISC_SORT}" \
    -o "${VCF}"

echo "[`date`] Done. VCF: ${VCF}"
