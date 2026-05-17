# Create the 100 kb hg38 bins BED
module load bedtools/2.31.1

CHROMSIZES="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/chrNameLength.txt"

OUTBASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Fragmentation_analysis/Regional_Fragment_analysis"
mkdir -p "$OUTBASE"

bedtools makewindows \
  -g "$CHROMSIZES" \
  -w 100000 \
| awk '$1!="chrM"' > "$OUTBASE/hg38.100kb.bed"
