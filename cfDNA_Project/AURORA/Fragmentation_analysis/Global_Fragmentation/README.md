# Global Fragmentation Analysis

This folder contains scripts for global cfDNA fragmentation analysis using paired-end BAM files. The current workflow uses Picard `CollectInsertSizeMetrics` to derive insert-size distributions and read-orientation metrics for each sample.

## Purpose

The goal of this analysis is to summarize global fragment-size patterns across cfDNA samples. For each BAM file, the workflow generates:

- Insert-size metrics table
- Insert-size histogram PDF

These outputs can be used to compare global cfDNA fragmentation profiles across samples or cohorts.

## Main script

```bash
CollectInsertSizeMetrics.sh
```

This is a Slurm array job script. Each array task processes one BAM file from the sample list.

## Software requirements

The script was written for an HPC environment using environment modules.

## For Downstream analysis, 
- follow *Global_fragmentation_analysis.ipynb*
