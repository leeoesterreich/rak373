#!/bin/bash
# This Script is used to calculate HRD Score (HRD-LOH+ TAI+LST)

#SBATCH --job-name=HRDsignatures
#SBATCH --time=30:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# Define variables correctly
calculate_HRD_LST_LOH_TAI_score="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/calculate_HRD_LST_LOH_TAI_score_smap_cnv_3.7.pl"
centromeres="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/hg38_centromeres.txt"
chrom_length="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/hg38_chrom_length.txt"
cnv_masks="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/hg38_cnv_masks.bed"

# Path to the sample list (having only sample name)
sample_list="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/Sample.txt"

# Read each sample name from the file and execute the script
while IFS= read -r sample; do
  echo "Processing sample: $sample"

  # Define file paths for this sample 
  smapFile="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/smap_files/${sample}_variants_combine_filters_inMoleRefine1.smap"
  cnvFile="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/cnvFile/${sample}_cnv_calls_exp.txt"
  aneuploidyFile="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/aneuploidyFile/${sample}_cnv_chrAneuploidy_exp.txt"
  outputFile="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/HRD_scores/${sample}_HRD_scores.txt"
  logFile="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/HRD/HRD_scores/logs/${sample}_run_log.txt"

  # Check if input files exist before running
  if [[ -f "$smapFile" && -f "$cnvFile" && -f "$aneuploidyFile" ]]; then
    # Run Perl script
    perl "$calculate_HRD_LST_LOH_TAI_score" \
      -smapFile "$smapFile" \
      -cnvFile "$cnvFile" \
      -aneuploidyFile "$aneuploidyFile" \
      -outFile "$outputFile" \
      -centromereFile "$centromeres" \
      -chromLengthFile "$chrom_length" \
      -cnvMaskFile "$cnv_masks" &> "$logFile"

    echo "Completed processing for sample: $sample"
  else
    echo "Skipping $sample: Missing input files."
  fi

done < "$sample_list"

echo "All samples processed."
