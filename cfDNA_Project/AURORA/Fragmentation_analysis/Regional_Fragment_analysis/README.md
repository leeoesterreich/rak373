# Regional Fragmentation Analysis

This folder contains the workflow for regional cfDNA fragmentation analysis using fixed 100 kb genomic bins. The analysis quantifies short and long cfDNA fragments separately, performs bin-level quality control, applies GC-bias correction, and computes the short/long fragmentation ratio per genomic bin.

## Overview

The workflow uses sorted coordinate BAM files as input and generates a regional fragmentation matrix where:

- rows = samples
- columns = genomic bins
- values = GC-corrected short/long fragmentation ratio

Fragment-length definitions used in this workflow:

| Fragment class | Insert size range |
|---|---:|
| Short fragments | 90-150 bp |
| Long fragments | 151-220 bp |

Fragments outside these ranges are not used for the short/long ratio calculation in this workflow.

## Reference files

- 100 kb genomic bin BED file

- hg38 blacklist BED file

- hg38 2bit genome file

## Software requirements

The workflow uses:

- Python 3
- deepTools
- pandas
- numpy
- py2bit
- statsmodels

On the cluster, deepTools was loaded as:

```bash
module load deeptools/3.5.6
```

## Step 1 - Split BAMs into chunks

BAM files are split into smaller chunk lists to allow parallel processing with a Slurm array job.

Each chunk file contains a subset of BAM paths, one BAM per line.

## Step 2 - Count reads per bin for short and long fragments

Script:

```bash
02B_crpb_chunks.sh
02B_crpb_chunks.py
```

The Slurm script runs one array task per BAM chunk. For each chunk, the Python script runs `deeptools.countReadsPerBin` twice:

1. short fragments: 90-150 bp
2. long fragments: 151-220 bp

The same fixed 100 kb BED file is used for all chunks.

Each raw TSV contains:

```text
chrom    start    end    sample_1    sample_2    ...
```

## Step 3 - Sort, verify, and merge chunk TSVs

Before merging chunk files column-wise, each raw TSV must be sorted by genomic coordinates.

This is important because each chunk contains different samples, but all chunks must have identical bin rows in identical order before merging.

Required check before merging:

```text
chrom/start/end rows must match exactly across all chunks
```
Expected merged outputs:

```bash
02_crpb/short_rawdata_countreadsperbin.tsv
02_crpb/long_rawdata_countreadsperbin.tsv
```

## Step 4 - Remove blacklist-overlapping bins

Merged short and long count matrices are filtered to remove bins overlapping hg38 blacklist regions.

Input files:

```bash
02_crpb/short_rawdata_countreadsperbin.tsv
02_crpb/long_rawdata_countreadsperbin.tsv
```

Output files:

```bash
02_crpb/short_rawdata_countreadsperbin.blacklistFiltered.tsv
02_crpb/long_rawdata_countreadsperbin.blacklistFiltered.tsv
```

The blacklist-filtered files are used for GC correction.

## Step 5 - GC-content calculation and LOESS correction

GC content is calculated for each retained genomic bin using the hg38 `.2bit` genome file.

Input files:

```bash
02_crpb/short_rawdata_countreadsperbin.blacklistFiltered.tsv
02_crpb/long_rawdata_countreadsperbin.blacklistFiltered.tsv
```

The workflow performs LOESS correction per sample using:

```python
frac = 0.75
```

For each sample, the correction is calculated as:

```python
corrected = raw_readcount - loess_fit + median(raw_readcount)
```

Expected outputs:

```bash
02_crpb/short_readcounts_corr_and_GC.csv
02_crpb/long_readcounts_corr_and_GC.csv
02_crpb/short_corr_readcounts.csv
02_crpb/long_corr_readcounts.csv
```

### Output descriptions

| File | Description |
|---|---|
| `short_readcounts_corr_and_GC.csv` | Bin-level short-fragment corrected counts with GC content |
| `long_readcounts_corr_and_GC.csv` | Bin-level long-fragment corrected counts with GC content |
| `short_corr_readcounts.csv` | Matrix of GC-corrected short-fragment counts, samples Ã— bins |
| `long_corr_readcounts.csv` | Matrix of GC-corrected long-fragment counts, samples Ã— bins |

## Step 6 - Compute short/long fragmentation ratio

Input files:

```bash
02_crpb/short_corr_readcounts.csv
02_crpb/long_corr_readcounts.csv
```

The short/long ratio is calculated per sample and per genomic bin:

```python
short_long_ratio = short_corrected_readcount / long_corrected_readcount
```

A small epsilon can be added to the denominator to avoid division by zero:

```python
eps = 1e-6
ratio = short / (long + eps)
```

Output:

```bash
02_crpb/short_long_ratio.csv
```

This file is the main regional fragmentation output for downstream analysis.

## Final expected output

The main final matrix is:

```bash
short_long_ratio.csv
```

Matrix structure:

```text
sample    chr1_0_100000    chr1_100000_200000    chr1_200000_300000    ...
```

Rows are samples and columns are genomic bins.

## Important QC checks

### 1. Confirm sorted bins match before merging

Before merging short or long chunk files, confirm that every sorted chunk has identical coordinate rows:

```text
chrom    start    end
```

If bins do not match, do not merge directly.

### 2. Confirm short and long matrices are aligned

After GC correction, short and long matrices must have:

- same samples
- same sample order
- same bins
- same bin order

Python checks:

```python
assert short_mat.index.equals(long_mat.index)
assert short_mat.columns.equals(long_mat.columns)
```

### 3. Check for zero or negative corrected long-fragment values

Because LOESS correction may produce zero or slightly negative values, inspect the long-fragment corrected matrix before calculating ratios.

```python
print("Long values <= 0:", int((long_df <= 0).sum().sum()))
```

If many values are zero or negative, consider masking those bins before downstream analysis.
