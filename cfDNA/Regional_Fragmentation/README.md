## Step 0 — BAM preprocessing (coordinate sort + filtering)

Generate coordinate-sorted BAM files containing high-confidence cfDNA fragments.
Filtering Criteria: proper pairs, MAPQ>=30, remove secondary/QCfail/dups/supplementary via -F 3844
```bash
module load samtools/1.15

samtools view -bh -f 2 -q 30 -F 3844 sample.bam \
  | samtools sort -o sample.filter.coordsort.bam -

samtools index sample.filter.coordsort.bam
