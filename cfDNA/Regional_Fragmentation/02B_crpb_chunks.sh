#!/bin/bash
# Run CountReadsPerBin per chunk
#SBATCH --job-name=CRPB_chunk
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH -t 1-00:00
#SBATCH --array=0-10
#SBATCH --output=slurm_log/Fragmentation/crpb_%A_%a.out
#SBATCH --error=slurm_log/Fragmentation/crpb_%A_%a.err

module load deeptools/3.5.6

BLACKLIST="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Blacklist_file/hg38-blacklist.v2.bed"

BINBED="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/hg38.100kb.bed"

CHUNKDIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/bamchunk"
OUTDIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/02_crpb"
PY="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/02B_crpb_chunks.py"

mkdir -p "$OUTDIR" slurm_log/Fragmentation

chunkfile=$(ls ${CHUNKDIR}/bamchunk_* | sort | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
if [[ -z "$chunkfile" ]]; then
  echo "No chunk for task ${SLURM_ARRAY_TASK_ID}"
  exit 0
fi

chunkid=$(basename "$chunkfile")
echo "Task ${SLURM_ARRAY_TASK_ID} processing chunk: $chunkfile"

for mode in short long; do
  outtsv="${OUTDIR}/${mode}.${chunkid}.raw.tsv"
  if [[ -s "$outtsv" ]]; then
    echo "Exists, skipping: $outtsv"
    continue
  fi

  python "$PY" "$chunkfile" "$outtsv" "$mode" 10 "$BLACKLIST" "$BINBED"
done
