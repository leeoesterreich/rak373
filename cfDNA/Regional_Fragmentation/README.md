Step 0 — BAM preprocessing (coordinate sort + filtering)
Goal

Generate coordinate-sorted BAM files containing high-confidence cfDNA fragments for downstream fragmentation analysis.

Input

Aligned BAM files (hg38):

/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/Aligned_bam/

Output

Filtered, coordinate-sorted BAM files:

/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/Sorted_bam/

Filtering criteria

Reads retained:

Proper pairs (-f 2)

MAPQ ≥ 30 (-q 30)

Reads removed:

Secondary

QC-fail

PCR/optical duplicates

Supplementary alignments

Flag used:

-F 3844

Per-sample command
module load samtools/1.15

in_bam="sample.bam"
out_bam="sample.filter.coordsort.bam"

samtools view -bh -f 2 -q 30 -F 3844 "$in_bam" \
  | samtools sort -o "$out_bam" -

samtools index "$out_bam"

Batch processing (all BAMs)
module load samtools/1.15

IN_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/Aligned_bam"
OUT_DIR="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/Sorted_bam"

mkdir -p "$OUT_DIR"

for in_bam in "$IN_DIR"/*.bam; do
  base=$(basename "$in_bam" .bam)
  out_bam="$OUT_DIR/${base}.filter.coordsort.bam"

  echo "Processing: $base"
  samtools view -bh -f 2 -q 30 -F 3844 "$in_bam" \
    | samtools sort -o "$out_bam" -

  samtools index "$out_bam"
done

Sanity checks

Read count before vs after filtering

samtools view -c sample.bam
samtools view -c sample.filter.coordsort.bam


Confirm coordinate sorting

samtools view -H sample.filter.coordsort.bam | grep '^@HD'

Output naming convention
<SAMPLE_ID>.filter.coordsort.bam
<SAMPLE_ID>.filter.coordsort.bam.bai


These BAM files are used as input for Step 1: fixed 100 kb genomic binning.
