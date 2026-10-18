#!/bin/bash
# Check RNA-seq library strandedness using STAR ReadsPerGene.out.tab files
# # 0% to <20% = Forward, 20% to 80% = Unstranded, >80% to 100% = Reverse-stranded 

# Input directory
dir="/ix1/alee/LO_LAB/Personal/Rahul/pILC/Chelsea_21ILC_RNASeq_data/2_alignment_file/"

# Output file
output="$dir/strand_check.tsv"

# Check for matching files
shopt -s nullglob
files=("$dir"/*_ReadsPerGene.out.tab)

if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: No ReadsPerGene.out.tab files found in $dir"
    exit 1
fi

# Analyze strandedness
for f in "${files[@]}"; do
    awk -v s="$(basename "$f")" '
        NR > 4 {
            y += $3
            r += $4
        }
        END {
            total = y + r

            if (total == 0) {
                strand = "undetermined"
                share = 0
            } else {
                share = r / total

                if (share > 0.8)
                    strand = "reverse"
                else if (share < 0.2)
                    strand = "yes"
                else
                    strand = "no (unstranded)"
            }

            printf "%s\treverse_share=%.2f\tcall=%s\n", s, share, strand
        }
    ' "$f"
done | tee "$output"

echo
echo "Summary of calls:"
cut -f3 "$output" | sort | uniq -c

echo
echo "Results saved to: $output"