# PDO Raw Counts (BO-046, BO-114, BO-189, BO-268)

## Overview
This folder contains the **raw RNA-seq gene counts** of four patient-derived organoids (PDOs): BO-046, BO-114, BO-189 and BO-268.

The counts were **extracted from the organoid RNA-seq count matrix** (Michael Gorry, `Org_RNAseq_Data`). They are used for VST and ssGSEA (Hallmark).

## Source files
Location: `/ix1/alee/LO_LAB/Personal/Michael_Gorry/Org_RNAseq_Data/`

| File | Use |
|------|-----|
| `Org_RAWcounts__ENSG_63140.csv` | Raw counts with Ensembl IDs (used for extraction) |
| `colData_organoid_metadata.xlsx` | Sample metadata (checked only; printed to the console) |

## Files in this folder

| File | Description |
|------|-------------|
| `PDO_raw_counts.csv` | Raw integer counts, genes (Ensembl IDs with version) x 4 samples |
| `PDO_sample_to_bam_map.csv` | BAM file from which each sample's counts came |
| `PDO_library_sizes.csv` | Total counts and number of genes with counts per sample |
| `Extract_PDO_raw_counts.R` | Script that creates the three files above |
| `PDO_ssGSEA_Hallmark.R` | Script for PDO-only Hallmark ssGSEA (reads `PDO_raw_counts.csv`) |

## How the counts were extracted

1. Read the raw count file as text.
2. Match each sample name (BO-046, BO-114, BO-189, BO-268) to exactly one column. Names are compared without case, dots, hyphens or underscores.
3. Remove the second header row (`Geneid` + BAM paths, featureCounts style). Without this step, R reads all columns as text. The BAM paths are saved in `PDO_sample_to_bam_map.csv`.
4. Convert to integers and check that no value was lost, all values are whole numbers and not negative, and gene IDs are unique.
5. Save the matrix, the BAM map and the library sizes.

## BAM files

| Sample | BAM name in the count file |
|--------|----------------------------|
| BO-046 | `BO046_S11_Aligned.sortedByCoord.out.bam` |
| BO-114 | `BO114_S9_Aligned.sortedByCoord.out.bam` |
| BO-189 | `BO189_S13_Aligned.sortedByCoord.out.bam` |
| BO-268 | `BO268_S36_BAligned.sortedByCoord.out.bam` |

## Library sizes

| Sample | Total counts | Genes with counts |
|--------|--------------|-------------------|
| BO-046 | 33,990,926 | 28,206 |
| BO-114 | 38,179,906 | 29,705 |
| BO-189 | 38,853,523 | 29,971 |
| BO-268 | 24,265,130 | 31,984 |

## Notes
- **BO-268:** the BAM name contains an extra `B` (`BO268_S36_B...`). It has the lowest depth but the most genes with counts. Check the metadata to see whether this is a batch, a second library or a re-sequencing run.
- Gene IDs keep their Ensembl version numbers (e.g. `ENSG00000290825.1`). Remove them before mapping to gene symbols.
- The matrix contains all annotated genes. Filtering (protein-coding genes, low counts, variance) is done in `PDO_ssGSEA_Hallmark.R`.
- ssGSEA scores from the PDO-only run are relative to these four samples. They cannot be compared directly with scores from the pILC tumor run.