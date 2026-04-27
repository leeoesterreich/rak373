# pILC WES Analysis Pipeline

This repository contains the SLURM-based workflow used for processing paired-end whole-exome sequencing (WES) data for the pILC project, starting from raw FASTQ files through somatic variant calling with GATK Mutect2 using a Panel of Normals (PoN).

## Overview

The pipeline includes the following major steps:

1. Raw FASTQ quality control
   - FastQC
   - MultiQC
2. Adapter trimming
   - Trim Galore
3. Post-trimming quality control
   - FastQC
   - MultiQC
4. Alignment to hg38
   - BWA-MEM
   - samtools sort/index
5. Alignment quality assessment
   - Qualimap
6. Preprocessing for variant calling
   - AddOrReplaceReadGroups
   - MarkDuplicates
   - BaseRecalibrator
   - ApplyBQSR
7. Panel of Normals (PoN) generation
   - Mutect2 on normal samples only
   - GenomicsDBImport
   - CreateSomaticPanelOfNormals
8. Tumor-only somatic variant calling using PoN
   - Mutect2
   - LearnReadOrientationModel
   - GetPileupSummaries
   - CalculateContamination
   - FilterMutectCalls

---

## Input data structure

Raw FASTQ files are organized into one folder per sample. Each sample may contain multiple lane-level FASTQ files.

Example:

```bash
TP08-356_T_OES8805A14_S20_L004_R1_001.fastq.gz
TP08-356_T_OES8805A14_S20_L003_R1_001.fastq.gz
TP08-356_T_OES8805A14_S20_L004_R2_001.fastq.gz
TP08-356_T_OES8805A14_S20_L003_R2_001.fastq.gz
```

Interpretation:

- `TP08-356_T` = sample ID
- `L003`, `L004` = sequencing lanes
- `R1`, `R2` = paired-end reads

Thus, **4 FASTQ files usually represent 1 biological sample**, not 4 samples.

---

## Sample naming convention

Sample folders and downstream files follow this convention:

- `*_T` = tumor
- `*_N` = normal

Examples:

```bash
TP13-1150_T
TP10-S1190_N
```

This naming is used throughout the workflow to automatically distinguish tumor and normal samples.

---

## Directory structure

### Raw data

```bash
/ix1/alee/LO_LAB/General/Lab_Data/20250611_HTGC_WES_pILC_Jagmohan/
```

### Project working directory

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/
```

### Main output structure

```bash
WES/
|-- 1_fastqc_result/
|-- 2_TrimGalore/
|-- 3_Alignment/
|-- 4_BQSR/
|-- 5_Mutect2/
    |-- PON_normals/
    |-- PON/
|-- 6_Mutect2_with_PON/
```

---

## Reference files

### Reference genome

```bash
/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Homo_sapiens_assembly38.fasta
```

### Known sites for BQSR

```bash
/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
```

### Germline resource for Mutect2 / contamination estimation

```bash
/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/af-only-gnomad.hg38.vcf.gz
```

### WES target intervals

```bash
/ix1/alee/LO_LAB/General/References/GATK_resource_bundle_hg38_060325_bts76/v0/Illumina_Exome_TargetedRegions_v1.2.hg38.bed
```

> Note: For WES PoN generation and Mutect2 calling, the exome target intervals should be used consistently. Do not use whole-genome calling intervals for WES PoN construction.

---

## Software and modules

Examples of modules used in this workflow:

```bash
fastqc/0.11.7
multiqc/1.30
trimgalore/0.6.10
bwa/0.7.17
samtools/1.22.1
qualimap/2.3
singularity/4.3.2
```

GATK is run through a Singularity container:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/env/gatk_latest.sif
```

---

# Pipeline steps

## 1. Raw FASTQ QC

FastQC is run on all raw lane-level FASTQ files in each sample folder.

Output:

```bash
WES/1_fastqc_result/<sample>/
```

MultiQC is then run to summarize all raw FastQC reports.

## 2. Adapter trimming

Lane-level FASTQ files for each sample are merged into one R1 and one R2 file, followed by trimming with Trim Galore.

Output:

```bash
WES/2_TrimGalore/<sample>/
```

Typical trimmed files:

```bash
<sample>_R1.merged_val_1.fq.gz
<sample>_R2.merged_val_2.fq.gz
```

## 3. Post-trim QC

FastQC is run again on trimmed reads, followed by MultiQC.

Output:

```bash
WES/2_TrimGalore/Fastqc/<sample>/
```

## 4. Alignment

Trimmed paired reads are aligned to hg38 using BWA-MEM.

Output:

```bash
WES/3_Alignment/<sample>/<sample>.sorted.bam
WES/3_Alignment/<sample>/<sample>.sorted.bam.bai
```

## 5. Alignment QC

Qualimap is used for BAM quality assessment.

Output:

```bash
WES/3_Alignment/Qualimap/<sample>/
```

Typical files:

- `qualimapReport.html`
- `genome_results.txt`
- `raw_data_qualimapReport/`
- `images_qualimapReport/`

## 6. BQSR preprocessing

Each aligned BAM is processed through:

1. `AddOrReplaceReadGroups`
2. `MarkDuplicates`
3. `BaseRecalibrator`
4. `ApplyBQSR`

Output:

```bash
WES/4_BQSR/<sample>/
```

Typical files:

```bash
<sample>.dedup.bam
<sample>.dedup.bai
<sample>.dedup.metrics.txt
<sample>.recal.table
<sample>.bqsr.bam
<sample>.bqsr.bai
```

### Index naming convention

This workflow uses **GATK-style BAM index files**:

```bash
<sample>.bqsr.bai
```

instead of:

```bash
<sample>.bqsr.bam.bai
```

Intermediate read-group files are removed after successful BQSR:

```bash
<sample>.rg.bam
<sample>.rg.bai
```

## 7. Normal-only Mutect2 calls for PoN generation

Normal samples (`*_N`) are processed using Mutect2 in tumor-only mode with:

```bash
--max-mnp-distance 0
```

This is important for compatibility with `GenomicsDBImport` during PoN generation.

Output:

```bash
WES/5_Mutect2/PON_normals/<sample>/<sample>.vcf.gz
```

## 8. Panel of Normals generation

PoN is generated in two steps:

1. `GenomicsDBImport`
2. `CreateSomaticPanelOfNormals`

Input:

- normal-only Mutect2 VCFs from `PON_normals/`
- optionally additional PoN-compatible normal VCFs from other technically similar WES datasets

Output:

```bash
WES/5_Mutect2/PON/<pon_name>.vcf.gz
```

### Important notes for PoN generation

- Input VCFs must be generated from **normal samples only**
- Input VCFs must be compatible with PoN generation
- VCFs containing MNPs can cause `GenomicsDBImport` failure
- Use `--max-mnp-distance 0` during Mutect2 normal calling
- Use **WES target intervals**, not WGS calling intervals

## 9. Tumor-only Mutect2 using PoN

Tumor samples (`*_T`) are called with:

- `--panel-of-normals`
- `--germline-resource`
- WES target intervals

Additional filtering workflow includes:

1. `LearnReadOrientationModel`
2. `GetPileupSummaries`
3. `CalculateContamination`
4. `FilterMutectCalls`

Output:

```bash
WES/5_Mutect2/Mutect2_with_PON/<sample>/
```

Typical files:

```bash
<sample>.unfiltered.vcf.gz
<sample>.unfiltered.vcf.gz.stats
<sample>.f1r2.tar.gz
<sample>.orientation-priors.tar.gz
<sample>.pileups.table
<sample>.contamination.table
<sample>.segments.table
<sample>.vcf.gz
```

Here:

- `<sample>.unfiltered.vcf.gz` = raw Mutect2 calls
- `<sample>.vcf.gz` = filtered final somatic callset

---

## SLURM usage

Example submission:

```bash
sbatch script_name.sh
```

Example array usage:

- 48 total WES samples:
  ```bash
  #SBATCH --array=0-47
  ```

- 17 normal samples:
  ```bash
  #SBATCH --array=0-16
  ```

- 31 tumor samples:
  ```bash
  #SBATCH --array=0-30
  ```

Optional concurrency limit:

```bash
#SBATCH --array=0-30%10
```

---

## Notes and troubleshooting

### 1. Mixed tumor and normal folders in the same parent directory

This is expected. Tumor and normal samples are separated automatically by suffix:

- `*_T` for tumor workflows
- `*_N` for normal workflows

### 2. PoN generation failure with `GenomicsDBImport`

A common cause is an MNP-containing VCF in the normal set. Ensure all normal Mutect2 VCFs used for PoN were generated with:

```bash
--max-mnp-distance 0
```

### 3. WES vs WGS intervals

For exome data, always use the exome capture intervals consistently across:

- normal Mutect2 for PoN
- GenomicsDBImport
- tumor Mutect2
- GetPileupSummaries
- FilterMutectCalls

### 4. BAM index format

This workflow keeps GATK-style `.bai` indexes rather than `.bam.bai`.

---

## Summary

This repository implements an end-to-end pILC WES analysis workflow using:

- FastQC / MultiQC
- Trim Galore
- BWA
- samtools
- Qualimap
- GATK BQSR
- Mutect2 PoN generation
- tumor-only Mutect2 with PoN-based somatic calling

The workflow is designed for SLURM-based HPC execution and supports automatic separation of tumor and normal samples based on folder naming conventions.

---