## Step 0 — BAM preprocessing (coordinate sort + filtering)

Generate coordinate-sorted BAM files containing high-confidence cfDNA fragments.

Filtering Criteria: proper pairs, MAPQ>=30, remove secondary/QCfail/dups/supplementary via -F 3844
```bash
module load samtools/1.15

samtools view -bh -f 2 -q 30 -F 3844 sample.bam \
  | samtools sort -o sample.filter.coordsort.bam -

samtools index sample.filter.coordsort.bam
```

## Step 1 — Generate fixed 100 kb genomic bins (hg38)

Create identical 100 kb genomic bins across the genome for all samples (excluding chrM).

```bash
module load bedtools/2.31.1

CHROMSIZES="/ix1/alee/LO_LAB/References/hg38/hg38.chrom.sizes"

bedtools makewindows -g "$CHROMSIZES" -w 100000 \
  | awk '$1!="chrM"' > hg38.100kb.bed
```
## Step 2 — Count reads per 100 kb bin

Count reads per bin across many BAMs by splitting the BAM list into chunks (avoids time/memory crashes).

```bash
# Working directory
cd /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/02_crpb

# 2A) Make BAM list (filtered + coordinate-sorted BAMs)
ls /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/Sorted_bam/*.filter.coordsort.bam \
  > bams.coordsort.list.txt

# Split into chunks of 25 BAMs each
split -d -l 25 bams.coordsort.list.txt bamchunk_   # creates bamchunk_00, bamchunk_01, ...

# 2B) Run CountReadsPerBin per chunk (produce short + long raw TSV per chunk)

```
## Step 3 — Enforce consistent bin order and verify short/long bin matching before merging chunk TSVs into final short/long matrices

Sort each chunk TSV in a fixed chromosome order (chr1–chr22 + chrX) and confirm that all chunk files (and short vs long) share identical bin rows.

```bash
cd /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/02_crpb

# Create a chr order file (chr1–chr22 + chrX)
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

# Sort each TSV using chr_order + genomic coordinates (keeps header)
for f in short.bamchunk_*.raw.tsv long.bamchunk_*.raw.tsv; do
  out="${f%.raw.tsv}.sorted.tsv"
  echo "Sorting $f -> $out"

  # header
  head -n 1 "$f" > "$out"

  # body: add chr rank, sort, then drop rank
  tail -n +2 "$f" \
    | awk -F'\t' 'BEGIN{OFS="\t"} FNR==NR {rank[$1]=NR; next}
        {r=(($1 in rank)?rank[$1]:999); print r,$0}
      ' chr_order.txt - \
    | sort -t $'\t' -k1,1n -k3,3n -k4,4n \
    | cut -f2- >> "$out"
done

# Verify bins match AND order matches across short chunks
ref=short.bamchunk_00.sorted.tsv
for f in short.bamchunk_*.sorted.tsv; do
  echo -n "$f  "
  diff -q <(cut -f1-3 "$ref" | tail -n +2) <(cut -f1-3 "$f" | tail -n +2) >/dev/null \
    && echo "OK" || echo "DIFF"
done

# Verify short vs long bins match
diff -q \
  <(cut -f1-3 short.bamchunk_00.sorted.tsv | tail -n +2) \
  <(cut -f1-3 long.bamchunk_00.sorted.tsv  | tail -n +2) \
  && echo "short/long bins match" || echo "short/long bins DIFF"
```
## Step 4 — Blacklist filtering, GC content, and LOESS GC correction

Remove bins overlapping genomic blacklists, compute GC content per bin, and apply LOESS correction per sample to reduce GC bias before calculating fragmentation ratios.

```text
# (A) Drop blacklist-overlapping bins
# SHORT_IN : short_rawdata_countreadsperbin.tsv;
# LONG_IN : long_rawdata_countreadsperbin.tsv
# BLACKLIST: hg38-blacklist.v2.bed
# Output: *.blacklistFiltered.tsv

# (B) Compute GC content per bin
# Input: OUT_BINS + reference FASTA
# Output: hg38.100kb.noblacklist.gc.tsv

# (C) LOESS GC correction (per sample)
# Apply LOESS correction to short and long count matrices using GC per bin.
# Output examples:
short.gc_corrected.tsv
long.gc_corrected.tsv
```

## Step 5 — Compute short/long fragmentation ratio

Compute the fragmentation signal per bin per sample using the short/long fragment ratio transformed to log2 for downstream comparisons.

```text
# Raw CountReadsPerBin outputs
short_rawdata_countreadsperbin.tsv
long_rawdata_countreadsperbin.tsv

# Blacklist-filtered raw counts
short_rawdata_countreadsperbin.blacklistFiltered.tsv
long_rawdata_countreadsperbin.blacklistFiltered.tsv

# GC content + raw counts (before LOESS correction)
short_readcounts_corr_and_GC.csv
long_readcounts_corr_and_GC.csv

# LOESS GC-corrected read counts
short_corr_readcounts.csv
long_corr_readcounts.csv

# Final fragmentation signal
short_long_ratio.csv
```
## Step 6 — CNV masking and consensus copy-neutral bin universe (ichorCNA)

Identify copy-neutral bins using ichorCNA segments and define consensus neutral-bin universes based on the fraction of samples in which each bin is neutral.

```bash
# Input:
#   short_long_ratio.csv
#   ichorCNA CNV segments (per sample)
#
# CNV masking logic:
#   - CNV status is determined per sample
#   - Each bin is marked neutral (1) or non-neutral (0) per sample
#   - A bin × sample mask matrix is generated
#
# Output files:
#   cnv_mask_matrix_bins_x_samples.tsv
#     Binary matrix indicating copy-neutral status (bins × samples)
#
#   neutral_fraction_per_bin.tsv
#     Fraction of samples in which each bin is copy-neutral
#
#   consensus_neutral_bins_<XX>pct.txt
#     Lists of bins that are copy-neutral in at least XX% of samples
#     (e.g., 50%, 60%, 70%)
#
#   consensus_bin_counts_summary.tsv
#     Summary of number and percentage of bins retained at each threshold
```
## Step 7 — Effect size analysis and site directionality (KW + BH)

Quantify site-specific fragmentation differences per bin using effect size, rather than relying only on statistical significance.

```bash
# Input:
#   short_long_ratio.cnv_neutral_universe50.csv
#   sample metadata with metastatic site labels (id,metastatic_sites)
#   
# Analysis logic:
#   - For each bin, compare fragmentation ratios across metastatic sites
#   - Perform Kruskal–Wallis (KW) test per bin
#   - Compute effect size (η²) to quantify strength of site-specific bias
#   - Apply Benjamini–Hochberg (BH) correction across bins
#
# Additional per-bin summaries:
#   - median_liver, median_brain, median_bone
#   - site_highest: site with highest median fragmentation signal
#   - site_lowest: site with lowest median fragmentation signal
#   - delta_high_low: difference between highest and lowest site medians
#
# Output example:
#   effect_size_kw_bh_results.tsv
#   (includes: kw_p, kw_q_BH, eta2, rank_eta2, site_highest, site_lowest)
```
### Step 8 — Visualization (Z-scored heatmaps) 

## Step 9 — Gene annotation of significant fragmentation bins

Annotate BH-significant fragmentation bins with overlapping genes to enable biological interpretation and downstream enrichment analyses.
  - Make a BED of BH-significant bins
  - Intersect those bins with a gene annotation BED to get overlapping gene names

```bash
module load bedtools/2.31.1

BINS="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/LOLA_inputs_top3sites/Effect_size/KW_BH/universe50.BHsig.bins986.bed"

GENES="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/gencode.v19.genes.bed"

OUT="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/LOLA_inputs_top3sites/Effect_size/KW_BH/universe50.BHsig.bins986.geneHits.tsv"

bedtools intersect -a "$BINS" -b "$GENES" -wa -wb > "$OUT"

echo "Wrote: $OUT"
```

# Intersect BH-significant bins with genes (NO CNV masking)

```bash
module load bedtools/2.31.1

BINS="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/LOLA_inputs_top3sites/Effect_size_no_cnv_masking/KW_BH/no_cnv_masking.BHsig.bins.bed"

GENES="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/gencode.v19.genes.bed"

OUT="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/Fragmentation_analysis/Regional_Fragment_analysis/LOLA_inputs_top3sites/Effect_size_no_cnv_masking/KW_BH/no_cnv_masking.BHsig.bins.geneHits.tsv"

bedtools intersect -a "$BINS" -b "$GENES" -wa -wb > "$OUT"

echo "Wrote: $OUT"
```
