#!/bin/bash
# End motif processing (protocol steps 06-endmotif_bed, 07-motif_fasta, 08-motif_merge)
# Input: filter_bedpe/*_filtered.bedpe  (Step 3 output)
# Output:
#   06-endmotif_bed/*_{nmer}bp_r1.bed and *_r2.bed
#   07-motif_fasta/*_fa.bed
#   08-motif_merge/*_{nmer}bp_motif.bed

#SBATCH --job-name=EndMotif_step4
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_log/Endmotif/Endmotif_step4_%A_%a.out
#SBATCH --error=slurm_log/Endmotif/Endmotif_step4_%A_%a.err
#SBATCH -t 1-00:00
#SBATCH --array=0-261%50

set -euo pipefail

# -------------------------
# Modules
# -------------------------
module purge
module load bedtools/2.31.1

# -------------------------
# Paths
# -------------------------
designated_path="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif"

bedpe_dir="${designated_path}/filter_bedpe"      # Step 3 output
bed_out="${designated_path}/endmotif_bed"        # protocol 06
fasta_out="${designated_path}/motif_fasta"       # protocol 07
merge_out="${designated_path}/motif_merge"       # protocol 08

mkdir -p "$bed_out" "$fasta_out" "$merge_out"

# hg38 reference fasta (use the same fasta you used in step 2)
reference="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/Human_genome_hg38_p7_chr_only.fna"

# n-mer length
nmer=4

# -------------------------
# Collect inputs
# -------------------------
bedpe_files=( "$bedpe_dir"/*_filtered.bedpe )

if [[ ${#bedpe_files[@]} -eq 0 || ! -e "${bedpe_files[0]}" ]]; then
  echo "No *_filtered.bedpe files found in: $bedpe_dir"
  exit 1
fi

task_id=${SLURM_ARRAY_TASK_ID}
if (( task_id < 0 || task_id >= ${#bedpe_files[@]} )); then
  echo "Invalid SLURM_ARRAY_TASK_ID=${task_id} for ${#bedpe_files[@]} files"
  exit 1
fi

f="${bedpe_files[$task_id]}"
base=$(basename "$f" "_filtered.bedpe")   # sample id
echo "Task: ${task_id}"
echo "Input: $f"
echo "Sample: $base"
echo "nmer: $nmer"

# -------------------------
# 06-endmotif_bed: make r1/r2 BEDs (8 columns)
# BEDPE columns (from your pipeline):
# 1 chr1, 2 start1, 3 end1, 4 chr2, 5 start2, 6 end2, 7 QNAME, 8 score, 9 strand1, 10 strand2, 11 fraglen, 12 GC
# r1 BED: chr start (start+nmer) QNAME score strand1 fraglen GC
# r2 BED: chr (end-nmer) end     QNAME score strand2 fraglen GC
# -------------------------
r1_bed="${bed_out}/${base}_${nmer}bp_r1.bed"
r2_bed="${bed_out}/${base}_${nmer}bp_r2.bed"

awk -v nmer="$nmer" 'BEGIN{OFS="\t"} {print $1, $2, $2+nmer, $7, $8, $9,  $11, $12}' "$f" > "$r1_bed"
awk -v nmer="$nmer" 'BEGIN{OFS="\t"} {print $4, $6-nmer, $6,    $7, $8, $10, $11, $12}' "$f" > "$r2_bed"

# -------------------------
# 07-motif_fasta: extract motif sequence (strand-aware)
# bedtools getfasta -s uses strand from BED col 6.
# Output columns after -bedOut: original BED + sequence as last column
# We keep: chr start end strand GC motif
# -------------------------
r1_fa="${fasta_out}/${base}_${nmer}bp_r1_fa.bed"
r2_fa="${fasta_out}/${base}_${nmer}bp_r2_fa.bed"

bedtools getfasta -fi "$reference" -bed "$r1_bed" -s -bedOut \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$6,$7,$8,toupper($9)}' > "$r1_fa"

bedtools getfasta -fi "$reference" -bed "$r2_bed" -s -bedOut \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$6,$7,$8,toupper($9)}' > "$r2_fa"

# -------------------------
# 08-motif_merge: merge r1 + r2 motif beds
# -------------------------
merged="${merge_out}/${base}_${nmer}bp_motif.bed"
cat "$r1_fa" "$r2_fa" > "$merged"

echo "DONE: ${base}"
echo "  r1_bed: $r1_bed"
echo "  r2_bed: $r2_bed"
echo "  r1_fa : $r1_fa"
echo "  r2_fa : $r2_fa"
echo "  merged: $merged"
