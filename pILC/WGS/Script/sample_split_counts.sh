#!/bin/bash

# To count split files per sample

# Set data directory
DATA_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20230717_jis117_pILC_split_Jagmohan"

# Output file path
OUTPUT_FILE="${SLURM_SUBMIT_DIR:-$HOME}/sample_split_R1_R2_counts.txt"

# Move to data directory
cd "$DATA_DIR" || { echo "Directory not found"; exit 1; }

# Temporary files
TMP_R1=$(mktemp)
TMP_R2=$(mktemp)

# Extract sample names and count R1 files
find . -name "*_R1_001*.fastq.gz" \
| sed -E 's|.*/||' \
| sed -E 's/_R1_001.*.fastq.gz$//' \
| sort \
| uniq -c \
| awk '{ printf "%s\tR1\t%s\n", $2, $1 }' > "$TMP_R1"

# Extract sample names and count R2 files
find . -name "*_R2_001*.fastq.gz" \
| sed -E 's|.*/||' \
| sed -E 's/_R2_001.*.fastq.gz$//' \
| sort \
| uniq -c \
| awk '{ printf "%s\tR2\t%s\n", $2, $1 }' > "$TMP_R2"

# Combine and pivot into one table
awk '
BEGIN { print "Sample_ID\tR1_Count\tR2_Count" }
{
    count[$1][$2] = $3
    samples[$1] = 1
}
END {
    for (s in samples) {
        r1 = (count[s]["R1"] ? count[s]["R1"] : 0)
        r2 = (count[s]["R2"] ? count[s]["R2"] : 0)
        printf "%s\t%s\t%s\n", s, r1, r2
    }
}
' "$TMP_R1" "$TMP_R2" > "$OUTPUT_FILE"

# Cleanup
rm "$TMP_R1" "$TMP_R2"

# Show result
echo "Sample R1/R2 split counts written to: $OUTPUT_FILE"
cat "$OUTPUT_FILE"
