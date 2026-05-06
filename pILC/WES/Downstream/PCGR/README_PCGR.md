# PCGR WES Annotation and Tumor-only vs With_PON Tier Comparison

This repository documents the WES downstream workflow used for **PCGR annotation** and subsequent **variant confidence tiering** by comparing two Mutect2 calling modes:

1. **Tumor-only mode**
2. **With_PON mode** — Mutect2 calling with a Panel of Normals

The workflow has two major parts:

```text
Part 1: Tumor-only / With_PON" VCF preparation and PCGR annotation
Part 2: Tumor-only vs With_PON PCGR PASS VCF comparison
```

---

# Part 1: Tumor-only / With_PON VCF preparation and PCGR annotation

The analysis begins with filtered Mutect2 tumor-only / With_PON VCFs. These VCFs are reformatted for PCGR compatibility by adding tumor and normal depth/allele-frequency INFO tags. The reformatted VCFs are then used as PCGR input. After PCGR annotation is complete for both Tumor-only and With_PON workflows, the PCGR PASS VCFs are compared to classify variants into confidence tiers.

```text
Filtered Tumor-only / With_PON VCFs
        │
        ▼
Reformat VCFs for PCGR
        │
        ▼
PCGR Tumor-only annotation
        │
        ▼
PCGR PASS VCFs
        │
        ├───────────────┐
        ▼               ▼
Tumor-only PASS     With_PON PASS
        │               │
        └───────┬───────┘
                ▼
        bcftools isec comparison
                │
        ┌───────┼────────┐
        ▼       ▼        ▼
    Tier 1   Tier 2    Tier 3
```
## Tumor-only or With_PON working directory

Expected directory structure:

```text
Tumor_only/
├── filtered_vcf/
├── filtered_vcf_file_list.txt
├── reformat_vcf.py
├── reformat_vcf.sh
├── vcf_pcgr/
├── vcf_pcgr_file_list.txt
├── tpm_pcgr/
├── cna_segs/
├── With_Tumor_only_index.csv
├── run_pcgr_Tumor_only_updated.sh
├── pcgr_outputs/
└── slurm_log/
```

---

## Step 1: Input filtered VCF list

The workflow starts from filtered tumor-only / With_PON VCFs listed in:

```bash
filtered_vcf_file_list.txt
```
---

## Step 2: Reformat tumor-only / With_PON VCFs for PCGR

### Script

```bash
reformat_vcf.sh
```

This is a Slurm array job:

```bash
#SBATCH --array=1-31
```

Each array task reads one VCF from:

```bash
filtered_vcf_file_list.txt
```

The script then runs:

```bash
reformat_vcf.py
```

inside a Singularity container.

### Output

Reformatted VCFs are written to:

```bash
vcf_pcgr/
```

## Step 3: Prepare PCGR input index

PCGR is run using a CSV index file:

```bash
With_Tumor_only_index.csv
With_PON_index.csv

```

This file is expected to contain one row per sample with:

```text
vcf_file,tpm_file,cna_file
```

The PCGR script reads the sample-specific row using the Slurm array task ID:

```bash
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$INDEX_FILE")
IFS=',' read -r vcf_file tpm_file cna_file <<< "$LINE"
sample_id=$(basename "$vcf_file" .taf.vcf.gz)
```

The `+1` means the script assumes the CSV file has a header row.

---

## Step 4: Run PCGR

### Script

```bash
run_pcgr_Tumor_only_updated.sh
run_pcgr_with_PON_updated.sh

```
---

## Optional RNA expression input

If a TPM file is provided in the index CSV and exists in:

```bash
tpm_pcgr/
```
the script adds RNA expression to the PCGR run

---

## Optional CNA segment input

If a CNA file is provided in the index CSV and exists in:

```bash
cna_segs/
```
the script adds CNA input to the PCGR run 

---

## PCGR command summary

PCGR is run inside a Singularity container with the following key settings:

```bash
pcgr \
  --input_vcf "/mnt/pcgr_inputs/$vcf_file" \
  --vep_dir "/mnt/.vep" \
  --refdata_dir "/mnt/bundle" \
  --output_dir "/mnt/pcgr_outputs/$sample_id" \
  --genome_assembly "grch38" \
  --sample_id "$sample_id" \
  --assay "WES" \
  --tumor_only \
  --tumor_site 6 \
  --tumor_dp_tag "TDP" \
  --tumor_af_tag "TAF" \
  --estimate_signatures \
  --estimate_msi \
  --estimate_tmb \
  --vcf2maf \
  --tumor_dp_min 4 \
  --tumor_af_min 0.01 \
  --force_overwrite
```

Additional RNA and CNA flags are appended when files are available.

---

## Run PCGR tumor-only annotation

```bash
cd /ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Tumor_only
sbatch run_pcgr_Tumor_only_updated.sh
```

---

## PCGR tumor-only outputs

Outputs are written to:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Tumor_only/pcgr_outputs
```

---

# Part 2: PCGR Tumor-only vs With_PON Tier Comparison

## Goal

After PCGR annotation is completed for both Tumor-only and With_PON workflows, the PASS VCFs are compared to classify variants by confidence.

The comparison uses:

```bash
bcftools isec
```

with exact matching of:

```text
CHROM + POS + REF + ALT
```

---

## Input directories for comparison

### Tumor-only PCGR outputs

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Tumor_only/pcgr_outputs
```

### With_PON PCGR outputs

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/With_PON/pcgr_outputs
```

---

## Output directory for comparison

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Compare_TumorOnly_vs_WithPON
```

Expected structure:

```text
Compare_TumorOnly_vs_WithPON/
├── sample_list.txt
├── Tier_1_common/
├── Tier_2_With_PON_only/
├── Tier_3_Tumor_only_only/
├── summaries/
├── logs/
├── temp/
└── summary_all_samples_tier_counts.csv
```

---

## Tier definitions

### Tier 1: common variants

```text
Tier 1 = Tumor-only ∩ With_PON
```

Variants present in both Tumor-only and With_PON PCGR PASS VCFs.

This is the highest-confidence group.

Output:

```bash
Tier_1_common/<SAMPLE>.tier1.common.vcf.gz
```

For common variants, the comparison script writes the record from the With_PON VCF so that the PON-informed PCGR annotation is retained.

---

### Tier 2: With_PON-only variants

```text
Tier 2 = With_PON − Tumor-only
```

Variants present only in the With_PON PCGR PASS VCF.

Output:

```bash
Tier_2_With_PON_only/<SAMPLE>.tier2.with_pon_only.vcf.gz
```

---

### Tier 3: Tumor-only-only variants

```text
Tier 3 = Tumor-only − With_PON
```

Variants present only in the Tumor-only PCGR PASS VCF.

This group should be interpreted cautiously because these variants were not retained in the PON-informed workflow.

Output:

```bash
Tier_3_Tumor_only_only/<SAMPLE>.tier3.tumor_only_only.vcf.gz
```

---

## Create sample list for comparison

Create:

```bash
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Compare_TumorOnly_vs_WithPON/sample_list.txt
```

## Run comparison

```bash
sbatch /ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/compare_tumor_only_vs_pon_tiers.sh
```
---

## Per-sample comparison outputs

For each sample, the script generates:

```bash
Tier_1_common/<SAMPLE>.tier1.common.vcf.gz
Tier_2_With_PON_only/<SAMPLE>.tier2.with_pon_only.vcf.gz
Tier_3_Tumor_only_only/<SAMPLE>.tier3.tumor_only_only.vcf.gz
```

Each VCF is indexed:

```bash
*.vcf.gz.tbi
```

## Combine all summaries

After all Slurm jobs finish:

```bash
cd /ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/PCGR/Compare_TumorOnly_vs_WithPON

awk -F',' '
BEGIN {
  OFS=",";
  print "Sample","Tumor-only PASS variants","With_PON PASS variants","Tier 1 common","Tier 2 With_PON only","Tier 3 Tumor-only only"
}
FNR==1 {next}
{
  print $1,$3,$4,$5,$6,$7
}
' summaries/*.tier_summary.csv > PCGR_TumorOnly_vs_WithPON_variant_counts.csv
```

---

# Recommended interpretation

## Tier 1

Tier 1 variants are present in both Tumor-only and With_PON PCGR PASS VCFs. These are the most robust variants for downstream interpretation.

## Tier 2

Tier 2 variants are present only in the With_PON PCGR PASS VCF. These variants are retained by the PON-informed workflow but were not present in the Tumor-only workflow.

## Tier 3

Tier 3 variants are present only in the Tumor-only PCGR PASS VCF. These variants should be interpreted cautiously because they are not retained when the PON-informed workflow is used.

---