#!/bin/bash
###############################################################################
# HTSeq gene counts for the 31 pILC tumors, one SLURM array task per BAM file.
#
# Library strandedness was measured from the STAR ReadsPerGene files:
# column 4 (reverse) holds ~87% of the counts in all 31 samples -> -s reverse.
#
# BEFORE SUBMITTING:   mkdir -p slurm_log/htseq
# SUBMIT:              sbatch run_htseq.sh
###############################################################################

#SBATCH --job-name=HTSeq
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -t 1-00:00
#SBATCH --array=0-30
#SBATCH --output=slurm_log/htseq/htseq_%A_%a.out
#SBATCH --error=slurm_log/htseq/htseq_%A_%a.err

# strand setting: yes / reverse / no (the job stops if the value is not one of these)
STRAND="reverse"

module purge
module load htseq/0.13.5
export PYTHONNOUSERSITE=1          # ignore numpy 2 in ~/.local; use the numpy that belongs to the HTSeq module

mkdir -p slurm_log/htseq

# ------------------------------- paths ---------------------------------------
input_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/2_alignment_file"
gtf_file="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/gencode.v47.basic.annotation.gtf"
output_dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded"   # new folder, old counts untouched
mkdir -p "$output_dir"

case "$STRAND" in
  yes|reverse|no) ;;
  *) echo "STRAND must be yes, reverse or no."; exit 1 ;;
esac

# ------------------------- pick the BAM of this task -------------------------
bam_files=("$input_dir"/*_Aligned.sortedByCoord.out.bam)
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#bam_files[@]}" ]; then
  echo "Task ${SLURM_ARRAY_TASK_ID}: no BAM file for this index (found ${#bam_files[@]} files)."; exit 1
fi
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "$bam_file" _Aligned.sortedByCoord.out.bam)
name=${sample%%_S[0-9]*}                # TP03-117_S25 -> TP03-117 ; names without _S stay unchanged
output_file="$output_dir/${name}.txt"

echo "Sample: $sample  ->  $output_file   (htseq-count -s $STRAND)"

# write to a temporary file first, so a failed run never leaves a half-finished count file
htseq-count -f bam -r pos -s "$STRAND" -m union --nonunique none -a 10 \
    "$bam_file" "$gtf_file" > "${output_file}.tmp" \
  && mv "${output_file}.tmp" "$output_file" \
  && echo "Done: $output_file"