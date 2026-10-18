## Chen, F., Onkar, S., Zou, J. et al.
## Immune infiltration correlates with transcriptomic subtypes in primary estrogen receptor positive invasive lobular breast cancer.
## npj Precision Oncology, 8, 257 (2024) (https://doi.org/10.1038/s41698-024-00746-z)

## Data Source
## The raw gene expression data were downloaded from the NCBI GEO:  
## Accession ID**: [GSE278301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278301)


#!/bin/bash
#SBATCH --job-name=Chelsea_21ILC
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-20%11
#SBATCH --output=slurm_log/Chelsea_%A_%a.out
#SBATCH --error=slurm_log/Chelsea_%A_%a.err

set -euo pipefail

# -------------------- Load modules --------------------
module purge
module load sra-tools/3.0.3

# -------------------- Define paths --------------------
BASE="/ix1/alee/LO_LAB/Personal/Rahul/pILC/Chelsea_21ILC_RNASeq_data"

ACCESSIONS="${BASE}/SRR_accessions.txt"
OUTPUT="${BASE}/Raw_Fastq"
SRA_DIR="${BASE}/SRA_downloads"
TEMP_BASE="${SLURM_TMPDIR:-${BASE}/Temp_SRA}"
TEMP_DIR="${TEMP_BASE}/Chelsea_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$OUTPUT" "$SRA_DIR" "$TEMP_DIR"

# -------------------- Read accession list --------------------
mapfile -t SAMPLES < <(grep -E '^SRR[0-9]+$' "$ACCESSIONS")

if [[ "${#SAMPLES[@]}" -ne 21 ]]; then
    echo "ERROR: Expected 21 SRR samples, found ${#SAMPLES[@]}"
    exit 1
fi

SRR="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "======================================"
echo "Sample: $SRR"
echo "Job: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "======================================"

# -------------------- Skip completed samples --------------------
if [[ -s "${OUTPUT}/${SRR}_1.fastq.gz" && \
      -s "${OUTPUT}/${SRR}_2.fastq.gz" ]] && \
      gzip -t "${OUTPUT}/${SRR}_1.fastq.gz" \
              "${OUTPUT}/${SRR}_2.fastq.gz"; then

    echo "Sample already downloaded and verified. Skipping."
    exit 0
fi

# -------------------- Download SRA --------------------
echo "Downloading SRA: $SRR"

prefetch "$SRR" \
    --output-directory "$SRA_DIR" \
    --max-size u

echo "SRA download complete: $SRR"

# -------------------- Convert SRA to FASTQ --------------------
echo "Starting FASTQ conversion: $SRR"

fasterq-dump "${SRA_DIR}/${SRR}" \
    --split-files \
    --threads "$SLURM_CPUS_PER_TASK" \
    --temp "$TEMP_DIR" \
    --outdir "$TEMP_DIR" \
    --progress

echo "FASTQ conversion completed: $SRR"

# -------------------- Validate and compress --------------------
for READ in 1 2; do

    FASTQ="${TEMP_DIR}/${SRR}_${READ}.fastq"

    if [[ ! -s "$FASTQ" ]]; then
        echo "ERROR: Missing or empty FASTQ: $FASTQ"
        exit 1
    fi

    echo "Compressing $FASTQ"

    gzip "$FASTQ"

    gzip -t "${FASTQ}.gz"

done

# -------------------- Move final files --------------------
mv "${TEMP_DIR}/${SRR}_1.fastq.gz" "$OUTPUT/"
mv "${TEMP_DIR}/${SRR}_2.fastq.gz" "$OUTPUT/"

# -------------------- Cleanup --------------------
rm -rf -- "$TEMP_DIR"

echo "======================================"
echo "Completed: $SRR"
echo "Finished: $(date)"
ls -lh "${OUTPUT}/${SRR}"_*.fastq.gz
echo "======================================"

