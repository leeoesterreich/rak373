#!/bin/bash

# This script is used to Call structural variants from GRIDSS (the Genomic Rearrangement IDentification Software Suite).

#SBATCH --job-name=gridss
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 6-00:00:00
#SBATCH --output=slurm_log/Gridss/gridss_run_%j.out
#SBATCH --error=slurm_log/Gridss/gridss_run_%j.err

set -euo pipefail

# ============================================================
# Load required modules
# ============================================================
module purge
module load gcc/8.2.0
module load gridss/2.13.2

# ============================================================
# Define paths and constants
# ============================================================
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich"
reference="/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg19_010626_rak373/hg19/ucsc.hg19.fasta"
exclude_list="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Blacklist_file/hg19-blacklist.3col.bed"

# Jar (exported so wrapper ALWAYS sees it)
export GRIDSS_JAR="/software/rhel9/manual/install/gridss/2.13.2/python3.11/share/gridss-2.13.2-5/gridss.jar"

pairs_csv="$designated_path/sample_sheet_010526.csv"
bqsr_dir="$designated_path/3_BQSR_010526"

output_dir="$designated_path/GRIDSS_010526"
log_dir="$designated_path/slurm_log/Gridss"
mkdir -p "$output_dir" "$log_dir"

status_log="$output_dir/gridss_status.log"

echo "[INFO] Host: $(hostname)"
echo "[INFO] gridss: $(which gridss)"
echo "[INFO] GRIDSS_JAR=$GRIDSS_JAR"
ls -lh "$GRIDSS_JAR"

# ============================================================
# Process each sample in CSV (Sample,BAM_Path)
# ============================================================
while IFS=',' read -r SAMPLE BAM_PATH; do
  # skip header / empty lines
  [[ "$SAMPLE" == "Sample" ]] && continue
  [[ -z "${SAMPLE// }" ]] && continue

  SAMPLE="$(echo "$SAMPLE" | xargs)"
  BAM_PATH="$(echo "${BAM_PATH:-}" | xargs)"

  # Use BQSR BAM as GRIDSS input
  TUMOR_BAM_PATH="${bqsr_dir}/${SAMPLE}/${SAMPLE}.bqsr.bam"

  echo "============================================================"
  echo "[INFO] Starting sample: $SAMPLE"
  echo "[INFO] CSV BAM_Path (alignment): $BAM_PATH"
  echo "[INFO] Using BQSR BAM: $TUMOR_BAM_PATH"
  echo "============================================================"

  sample_outdir="$output_dir/$SAMPLE"
  mkdir -p "$sample_outdir"

  workdir="$sample_outdir/work"
  tmpdir="$sample_outdir/tmp"
  mkdir -p "$workdir" "$tmpdir"
  export TMPDIR="$tmpdir"

  # Validate required files (BQSR BAM is the input)
  for f in "$TUMOR_BAM_PATH" "${TUMOR_BAM_PATH}.bai" "$reference" "$exclude_list" "$GRIDSS_JAR"; do
    if [[ ! -f "$f" ]]; then
      echo "[ERROR] Missing file: $f"
      echo "$SAMPLE,FAILED (missing $f),$(date)" >> "$status_log"
      continue 1
    fi
  done

  out_vcf="$sample_outdir/${SAMPLE}.gridss.vcf.gz"
  out_assembly="$sample_outdir/${SAMPLE}.gridss.assembly.bam"
  run_log="$sample_outdir/${SAMPLE}.gridss.stdout_stderr.log"

  if [[ -f "$out_vcf" ]]; then
    echo "[SKIP] $SAMPLE already processed: $out_vcf"
    echo "$SAMPLE,SKIPPED (vcf exists),$(date)" >> "$status_log"
    continue
  fi

  echo "[INFO] Running GRIDSS..."
  start_time=$(date +%s)

  # NOTE: GRIDSS CLI uses -t/--threads (not --threads)
  gridss \
    -r "$reference" \
    -o "$out_vcf" \
    -a "$out_assembly" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -j "$GRIDSS_JAR" \
    -w "$workdir" \
    -b "$exclude_list" \
    "$TUMOR_BAM_PATH" \
    > "$run_log" 2>&1

  status=$?
  end_time=$(date +%s)
  elapsed=$((end_time - start_time))

  if [[ $status -ne 0 ]]; then
    echo "[ERROR] GRIDSS failed for $SAMPLE (see $run_log)"
    echo "$SAMPLE,FAILED,Elapsed:${elapsed}s,$(date)" >> "$status_log"
  else
    echo "[SUCCESS] GRIDSS completed for $SAMPLE"
    echo "$SAMPLE,SUCCESS,Elapsed:${elapsed}s,$(date)" >> "$status_log"
  fi

  echo ""
done < <(tail -n +2 "$pairs_csv")
