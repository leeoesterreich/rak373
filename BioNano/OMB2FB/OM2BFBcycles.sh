#!/bin/bash
# This script is used to analyze Breakage Fusion Bridge (BFB) event

#SBATCH --job-name=BFB
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00

# Load necessary modules
module purge
module load gcc/12.2.0 
module load python/ondemand-jupyter-python3.10
module load java/11.0.2-openjdk
module load r/4.4.0

# Activate the Conda environment
source activate /ihome/alee/rak373/.conda/envs/BFBCycle

# Set the default Python interpreter
export PATH=/ihome/alee/rak373/.conda/envs/BFBCycle/bin:$PATH

# Grant execute permissions to the necessary script
chmod +x /bgfs/alee/LO_LAB/Personal/Rahul/BioNano/BFBCycle/runOM2BFB.py

# Set paths
runOM2BFB="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/BFBCycle/runOM2BFB.py"
input_dir="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/BFBCycle/file"
centromere_regions="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/BFBCycle/hg38_centro.txt"
output_dir="/bgfs/alee/LO_LAB/Personal/Rahul/BioNano/BFBCycle/Output"
bfbfinder="/ihome/crc/install/bfbfinder/BFBFinder/build/libs/BFBFinder.jar"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Run BFB analysis
echo "Starting analysis for each sample..."

for rmcap_file in "$input_dir"/*_cnv_rcmap_exp.txt; do
    sample=$(basename "$rmcap_file" _cnv_rcmap_exp.txt)

    smap_file=$(find "$input_dir" -name "${sample}*_variants_combine_filters_inMoleRefine1.smap" -print -quit)
    xmap_file=$(find "$input_dir" -name "${sample}*_EXP_REFINEFINAL1.xmap" -print -quit)
    cmap_file=$(find "$input_dir" -name "${sample}*_EXP_REFINEFINAL1_q.cmap" -print -quit)
    refine_label_dir=$(find "$input_dir" -name "${sample}*_refine1_ExperimentLabel" -print -quit)

    if [[ -z "$smap_file" || -z "$xmap_file" || -z "$cmap_file" || -z "$refine_label_dir" ]]; then
        echo "Missing one or more input files for sample ${sample}. Skipping..."
        continue
    fi

    output_name="${output_dir}/${sample}"

    echo "Running analysis for ${sample}..."

    python3 "$runOM2BFB" \
        -r "$rmcap_file" \
        -c "$centromere_regions" \
        -n "Bionano" \
        -o "$output_name" \
        -s "$smap_file" \
        -x "$xmap_file" \
        -fol "$refine_label_dir" \
        -cmap "$cmap_file" \
        -cov 100 \
        -bfbfinder "$bfbfinder"

    echo "Analysis for ${sample} completed."
done

echo "Analysis completed for all samples."
