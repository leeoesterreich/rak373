Usage: Regional cfDNA Fragmentation Analysis Pipeline

This section describes the end-to-end workflow used to perform regional cfDNA fragmentation analysis based on short/long fragment ratios across fixed genomic bins.

Step 0. BAM preprocessing

All BAM files must be coordinate-sorted and filtered to retain high-confidence cfDNA fragments.

Filtering criteria

Proper pairs only

MAPQ ≥ 30

Remove secondary, QC-fail, duplicate, and supplementary reads

Example:

samtools view -b \
  -f 2 \
  -q 30 \
  -F 3844 \
  input.bam \
| samtools sort -o sample.filter.coordsort.bam

samtools index sample.filter.coordsort.bam

Step 1. Define fixed genomic bins (hg38, 100 kb)

All samples are summarized over identical genomic bins.

CHROMSIZES=/path/to/hg38.chrom.sizes

bedtools makewindows -g "$CHROMSIZES" -w 100000 \
  | awk '$1!="chrM"' > hg38.100kb.bed


Bin size: 100 kb

chrM excluded

Step 2. Count reads per bin (chunked execution)

Counting reads for hundreds of BAMs in a single job can exceed time/memory limits.
To avoid this, BAMs are processed in chunks.

2A. Split BAM list
ls Sorted_bam/*.filter.coordsort.bam > bams.coordsort.list.txt
split -d -l 25 bams.coordsort.list.txt bamchunk_


Each bamchunk_* file contains ~25 BAMs.

2B. Run CountReadsPerBin per chunk

Run the bin-counting step independently per chunk

Ensure all chunks use the same bin BED file

Generate separate raw TSVs for:

short fragments

long fragments

2C. Merge chunk-level TSVs

After all chunks finish, merge them into:

short_rawdata.tsv

long_rawdata.tsv

Each TSV must share identical bin rows across chunks.

Step 2D. Enforce consistent chromosome and bin order

Create a fixed chromosome order:

cat > chr_order.txt <<'EOF'
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
EOF


Sort each chunk TSV while preserving headers:

for f in short.bamchunk_*.raw.tsv long.bamchunk_*.raw.tsv; do
  out="${f%.raw.tsv}.sorted.tsv"
  head -n 1 "$f" > "$out"
  tail -n +2 "$f" \
    | awk -F'\t' 'BEGIN{OFS="\t"} FNR==NR {rank[$1]=NR; next}
      {r=(($1 in rank)?rank[$1]:999); print r,$0}' chr_order.txt - \
    | sort -t $'\t' -k1,1n -k3,3n -k4,4n \
    | cut -f2- >> "$out"
done

Verify bin consistency
ref=short.bamchunk_00.sorted.tsv

for f in short.bamchunk_*.sorted.tsv; do
  diff -q <(cut -f1-3 "$ref" | tail -n +2) <(cut -f1-3 "$f" | tail -n +2) \
    && echo "$f OK" || echo "$f DIFF"
done

diff -q \
  <(cut -f1-3 short.bamchunk_00.sorted.tsv | tail -n +2) \
  <(cut -f1-3 long.bamchunk_00.sorted.tsv  | tail -n +2) \
  && echo "short/long bins match"

Step 2E. Blacklist filtering & GC correction

Remove bins overlapping genomic blacklists

Compute GC content per bin

Apply LOESS GC correction per sample

(This step produces GC-corrected short and long counts.)

Step 3. Short/long fragment ratio

For each bin and sample:

short_long_ratio = short_count / long_count


Optionally transform to:

log2(short/long)

Z-score per bin across samples

Step 4. CNV masking using ichorCNA

To avoid CNV-driven artifacts:

Perform CNV masking per sample

Retain only copy-neutral bins

Use a shared neutral-bin universe across all samples for comparisons

Step 5. Per-site comparison and effect size analysis

Compare bins across metastatic sites (e.g., liver, brain, bone)

Across all bins: no strong global pattern

Compute effect size (η²) per bin

Select top 1000 bins by effect size

Assign site directionality using median values:

median_liver

median_brain

median_bone

Track site_highest / site_lowest per bin

Step 6. Gene annotation of significant bins
6A. Create BED of BH-significant bins

Example: 986 bins passing BH correction

6B. Intersect bins with genes
module load bedtools/2.31.1

bedtools intersect \
  -a universe50.BHsig.bins986.bed \
  -b gencode.v19.genes.bed \
  -wa -wb > universe50.BHsig.bins986.geneHits.tsv

6C. Collapse gene hits

Collapse genes per bin_id

Merge gene annotations back into the bin-level TSV
