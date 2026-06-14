# All-case CNV extraction from annotated CNVkit CNR files

## Purpose

This folder contains all-case longitudinal gene-level CNV matrices generated from per-sample annotated CNVkit `.cnr` files. These files were created to extend the previous 20 selected-case clinical/CNV extraction workflow to all available AURORA cfDNA samples.

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/Annotated_cnr_file
```

This folder is expected to contain annotated CNR files for approximately 784 longitudinal plasma samples from approximately 270 AURORA cases.

## Output directory

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction/
```

## Generated files

### 1. Gene-level CNV state matrix

```text
All_representative_patients_longitudinal_gene_CNV_state_matrix.tsv
```

Despite the name containing `representative`, this file is intended to contain all samples processed from the full annotated CNR folder.

Structure:

```text
rows    = genes
columns = samples
values  = gene-level CNV state
```

The first column is:

```text
gene
```

All remaining columns are sample IDs, for example:

```text
AL1Y-NB1
AL1Y-NB2
ARWN-NB4
```

CNV states are standardized as:

```text
NEUT   = neutral / no called CNV change
GAIN   = copy-number gain
AMP    = amplification
HLAMP  = high-level amplification
HETD   = heterozygous deletion
HOMD   = homozygous deletion, if present
```

If a gene is represented by multiple genomic bins in a sample, the script keeps the strongest altered state using this priority:

```text
NEUT < HETD < HOMD < GAIN < AMP < HLAMP
```

This creates one CNV state per gene per sample.

---

### 2. CNV sample annotation file

```text
All_representative_patients_longitudinal_CNV_sample_annotation.tsv
```

This file contains one row per sample and is used by the downstream clinical/CNV extraction workflow.

Important columns include:

```text
sample_id
Broad_entity_sample_id
patient_id_from_metasheet
Patient_Identifier
Patient_Identifier_from_metasheet
Patient_Identifier_metasheet_match
draw_num
Tumor_Fraction
Tumor_Fraction_percent
Ploidy
Mean_coverage
Histological_subtype
Histology_group
ER_Status
PR_Status
HER2_Status
Molecular_subtype
Molecular_subtype_clean
Subclone Fraction
Frac. Genome Subclonal
Frac. CNA Subclonal
Frac_Subnucleosomal_90_150
Frac_Mononucleosomal_151_250
Frac_Dinucleosomal_251_450
Frac_Trinucleosomal_451_650
CNV_burden_non_NEUT_genes
n_amp_or_hlamp_genes
n_hlamp_genes
CNR_file
```

Tumor fraction is merged from:

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Clinical_data/Metasheet_Analysis.with_Histology_Receptor_status.xlsx
```

The metasheet stores tumor fraction as a fraction, for example:

```text
0.09397
```

The script keeps the original value and also creates a percent-scaled value:

```text
Tumor_Fraction         = 0.09397
Tumor_Fraction_percent = 9.397
```

Sample ID matching is done by normalizing the metasheet collaborator sample ID:

```text
AUP-AL1Y-NB1 -> AL1Y-NB1
```

This is then matched to the CNR-derived sample ID:

```text
AL1Y-NB1
```

---

### 3. CNR manifest file

```text
All_representative_patients_longitudinal_CNR_manifest.tsv
```

This file records the source CNR file and basic QC metrics for each sample.

Important columns include:

```text
sample_id
Patient_Identifier
draw_num
cnr_file
call_column_used
n_cnr_rows
n_gene_rows_after_explode
n_unique_genes
n_altered_genes_non_NEUT
n_amp_or_hlamp_genes
n_hlamp_genes
```

The `call_column_used` column records which CNV call column was used from each `.annotated.cnr` file. The preferred column is:

```text
<SAMPLE>.Corrected_Call
```

If that is not available, the script uses:

```text
<SAMPLE>.event
```

## Input files

### Annotated CNR files

Each sample has one annotated CNVkit `.cnr` file, for example:

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/Annotated_cnr_file/ARWN-NB4.annotated.cnr
```

Expected structure:

```text
chromosome
start
end
gene
log2
<SAMPLE>.Corrected_Call
<SAMPLE>.Corrected_Copy_Number
<SAMPLE>.copy.number
<SAMPLE>.event
<SAMPLE>.logR_Copy_Number
<SAMPLE>.subclone.status
```

The script uses the `gene` column and the corrected CNV call column.

### Metadata / tumor-fraction file

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Clinical_data/Metasheet_Analysis.with_Histology_Receptor_status.xlsx
```

Used columns include:

```text
Broad_entity:sample_id
Collaborator_sample_id
Patient_id
Patient_Identifier
Histological_subtype
ER_Status
PR_Status
HER2_Status
Molecular_subtype
Mean_coverage
Tumor_Fraction
Ploidy
Subclone Fraction
Frac. Genome Subclonal
Frac. CNA Subclonal
Frac_Subnucleosomal_90_150
Frac_Mononucleosomal_151_250
Frac_Dinucleosomal_251_450
Frac_Trinucleosomal_451_650
```

## Workflow summary

```text
1. Read all *.annotated.cnr files from Annotated_cnr_file/
2. Extract sample_id from filename
      ARWN-NB4.annotated.cnr -> ARWN-NB4
3. Extract Patient_Identifier and draw number
      ARWN-NB4 -> Patient_Identifier = ARWN, draw_num = 4
4. Read gene-level CNV calls from each CNR file
5. Split comma-separated gene names into individual genes
6. Collapse multiple bins per gene to one strongest CNV state
7. Combine all samples into one gene x sample CNV state matrix
8. Read metasheet and normalize sample IDs
      AUP-AL1Y-NB1 -> AL1Y-NB1
9. Merge tumor fraction, histology, receptor status, and molecular subtype
10. Save all-case CNV matrix, sample annotation, and CNR manifest
```

## How to use these files in the all-case clinical/CNV extraction script

After generating these files, update the clinical/CNV extraction script to use:

```python
cnv_matrix_file = BASE_DIR / "CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_gene_CNV_state_matrix.tsv"

cnv_annotation_file = BASE_DIR / "CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_CNV_sample_annotation.tsv"

cnr_manifest_file = BASE_DIR / "CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_CNR_manifest.tsv"
```

Keep:

```python
candidate_cases = None
```

This will make the downstream script use all cases present in the new all-case CNV matrix.

## Important notes

1. The earlier 20-case workflow was limited because it used files from:

```text
CNVKit/Selected_patient_CNV_heatmap/
```

Those files contained only selected representative cases.

2. `candidate_cases = None` means:

```text
Use all cases present in the input CNV matrix.
```

It does not mean all AURORA cases unless the input CNV matrix itself contains all cases.

3. The all-case matrix should now be used for full-cohort clinical/CNV extraction.

4. Tumor fraction is derived from the metasheet, not directly from the CNR files.

5. `Patient_Identifier_metasheet_match` should be checked as a QC column. It confirms whether the patient identifier parsed from the CNR sample name agrees with the patient identifier in the metasheet.

## Recommended QC checks

After generating the files, run:

```bash
cd /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction

wc -l All_representative_patients_longitudinal_CNV_sample_annotation.tsv
wc -l All_representative_patients_longitudinal_CNR_manifest.tsv

head -5 All_representative_patients_longitudinal_CNV_sample_annotation.tsv
head -5 All_representative_patients_longitudinal_CNR_manifest.tsv
```

In Python:

```python
import pandas as pd
from pathlib import Path

out_dir = Path("/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction")

mat = pd.read_csv(
    out_dir / "All_representative_patients_longitudinal_gene_CNV_state_matrix.tsv",
    sep="\t",
    nrows=5
)

ann = pd.read_csv(
    out_dir / "All_representative_patients_longitudinal_CNV_sample_annotation.tsv",
    sep="\t"
)

manifest = pd.read_csv(
    out_dir / "All_representative_patients_longitudinal_CNR_manifest.tsv",
    sep="\t"
)

print("CNV matrix preview shape:", mat.shape)
print("Annotation samples:", ann.shape[0])
print("Manifest samples:", manifest.shape[0])
print("Unique cases:", manifest["Patient_Identifier"].nunique())
print("Samples with TFx:", ann["Tumor_Fraction_percent"].notna().sum())
print("Patient ID match summary:")
print(ann["Patient_Identifier_metasheet_match"].value_counts(dropna=False))
```

## Downstream use

These files are intended as inputs for:

```text
All_cases_clinical_CNV_extraction.xlsx
```

and related TSV outputs that integrate:

```text
clinical data
blood draw timing
treatment timing
active treatment context
tumor fraction
gene-level CNV burden
longitudinal CNV change
```

The same strict treatment-context logic from the 20-case workflow should be used:

```text
Active_treatment_at_draw = treatment_start_day <= blood_draw_day <= treatment_end_day
```

Both treatment start and treatment end days must be present. Missing treatment end day is not treated as ongoing.
