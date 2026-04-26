# FACETS Copy Number Analysis Pipeline for WES Tumor–Normal Pairs

This folder contains the scripts used to run FACETS-based allele-specific copy number analysis for paired tumor–normal whole-exome sequencing (WES) samples. The workflow starts from matched tumor and normal BQSR BAM files, generates SNP pileup count files, runs FACETS/facetsSuite, and extracts major/minor copy number calls from FACETS RDS outputs for PCGR.
Ref: https://github.com/mskcc/facets-suite, https://github.com/mskcc/facets
---

## Overview

The pipeline has three main steps:

1. **Generate SNP pileup files** from matched tumor–normal BAMs using `snp-pileup`.
2. **Run FACETS/facetsSuite** on each `.snp_pileup.gz` file.
3. **Extract CNA segments (for PCGR input)** from FACETS RDS output into a tab-delimited format containing chromosome, start, end, major copy number, and minor copy number.

---

## Repository structure

```text
FACETS/
├── Paired_samples.csv
├── snp-pileup/
│   └── *.snp_pileup.gz
├── Result/
│   └── <SAMPLE_ID>/
├── slurm_log/
├── Script/
│   ├── snp-pileup-wrapper.R
│   ├── run-facets-wrapper.R
│   └── extract_cna_from_facets_rds.R
├── snp_aray_pileup.sh
└── facets.sh
```

---

## Input files

### 1. Paired sample sheet

`Paired_samples.csv` contains one row per matched tumor–normal pair.

Required columns:

```text
SAMPLE,TUMOR_BAM,NORMAL_BAM,TUMOR_BAM_PATH,NORMAL_BAM_PATH
```

Example:

```text
TP03-117,TP03-117_T,TP03-117_N,/path/to/TP03-117_T.bqsr.bam,/path/to/TP03-117_N.bqsr.bam
```

### 2. Tumor and normal BAM files

The pipeline expects BQSR-processed BAM files for each tumor and matched normal sample. The absolute BAM paths are read from `Paired_samples.csv`.

### 3. Common SNP VCF

The SNP pileup step uses a polymorphic common SNP VCF file from NCBI dbSNP build 151 for GRCh38.p7:

```text
https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
```

In this workflow, the downloaded VCF is stored at:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/00-common_all.vcf.gz
```

Example download command:

```bash
mkdir -p /ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human

wget -O /ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/00-common_all.vcf.gz \
  https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
```
---

## Software requirements

The scripts were configured for the HPC module environment:

```bash
module load gcc/12.2.0 r/4.4.0
```

Required R packages include:

```r
argparse
facets
facetsSuite
pctGCdata
dplyr
readr
```

For the clean R 4.4.0 environment, the R package library used was:

```bash
/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4
```

---

## Step 1: Generate SNP pileup files

Submit the SNP pileup SLURM array job:

```bash
sbatch snp_aray_pileup.sh
```

This script reads paired tumor–normal BAM paths from:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Paired_samples.csv
```

and writes output files to:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/snp-pileup
```

Each output file has the format:

```text
<SAMPLE>.snp_pileup.gz
```

---

## Step 2: Run FACETS/facetsSuite

After the `.snp_pileup.gz` files are generated, submit the FACETS SLURM array job:

```bash
sbatch facets.sh
```

The FACETS script reads input files from:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/snp-pileup
```

For each `.snp_pileup.gz` file, it creates a sample-specific output directory:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/FACETS/Result/<SAMPLE_ID>
```

The current FACETS run uses:

```bash
--purity-cval 150
--cval 150
--everything
```

and passes the R package library path using:

```bash
--facets-lib-path /ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.4
```

---

## Extract CNA segments from FACETS RDS output (For PCGR)

The script `extract_cna_from_facets_rds.R` converts FACETS RDS output into a PCGR-compatible TSV file with the following columns:

```text
Chromosome
Start
End
nMajor
nMinor
```

Usage:

```bash
Rscript Script/extract_cna_from_facets_rds.R <input_facets.rds> <output_cna.tsv>
```
The script calculates:

```r
minor_cn <- lcn.em
major_cn <- tcn.em - lcn.em
```

Missing `lcn.em` values are set to 0 before calculating major/minor copy number.

---

## Workflow summary

```bash
# 1. Generate SNP pileup files
sbatch snp_aray_pileup.sh

# 2. Run FACETS
sbatch facets.sh

# 3. Optional: extract CNA table from FACETS RDS (for PCGR).
Rscript Script/extract_cna_from_facets_rds.R <input_facets.rds> <output_cna.tsv>
```
