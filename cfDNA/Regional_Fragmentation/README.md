## Step 0 — BAM preprocessing (coordinate sort + filtering)

Generate coordinate-sorted BAM files containing high-confidence cfDNA fragments.

Filtering Criteria: proper pairs, MAPQ>=30, remove secondary/QCfail/dups/supplementary via -F 3844
```bash
module load samtools/1.15

samtools view -bh -f 2 -q 30 -F 3844 sample.bam \
  | samtools sort -o sample.filter.coordsort.bam -

samtools index sample.filter.coordsort.bam

## Step 1 — Generate fixed 100 kb genomic bins (hg38)

Create identical 100 kb genomic bins across the genome for all samples (excluding chrM).

```bash
module load bedtools/2.31.1

CHROMSIZES="/ix1/alee/LO_LAB/References/hg38/hg38.chrom.sizes"

bedtools makewindows -g "$CHROMSIZES" -w 100000 \
  | awk '$1!="chrM"' > hg38.100kb.bed
