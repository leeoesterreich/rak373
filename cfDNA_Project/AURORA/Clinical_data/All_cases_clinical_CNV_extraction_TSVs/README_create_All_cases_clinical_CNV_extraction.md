# README: `create_All_cases_clinical_CNV_extraction.ipynb`

## Purpose

This notebook generates an all-case clinical + longitudinal CNV/TFx extraction workbook and matching TSV files for AURORA cfDNA samples.

It extends the earlier 20 representative-case workflow to the full all-case CNV matrix generated from annotated CNVkit CNR files.

The notebook integrates:

```text
clinical data
blood draw timing
treatment timing
active treatment context
tumor fraction / TFx
histology and receptor status
gene-level CNV burden
selected key-gene CNV states
longitudinal CNV/TFx change
```

## Key design principles

1. The workflow starts from the CNV matrix, so CNV-only or clinically unmapped samples are retained.
2. `Patient_Identifier` is extracted from CNV sample IDs such as `AL1Y-NB1`.
3. Clinical mapping is performed using clinical `Barcode` values.
4. `Active_treatment_at_draw` uses a strict treatment interval rule:

```text
treatment_start_day <= blood_draw_day <= treatment_end_day
```

Both treatment start and treatment end days must be present.

5. Missing treatment end day is **not** assumed to mean ongoing treatment.
6. Additional treatment-context columns are added because active treatment status alone is not enough for interpretation.
7. For all-case extraction, the notebook uses a memory-safe wide-matrix approach instead of creating a very large `cnv_long` table.

## Input files

### 1. Clinical file

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Clinical_data/AUP_ClinicalData_20260609.NO.DATES.bloodID.UseSet.txt
```

This file provides clinical rows including:

```text
patient_id
Site
Category
Days_From_Primary_Diagnosis
BCR_Status
Clinical_ID
Clinical_ID_2
Barcode
Tissue_Site
Histologic.type
ER_Status
PR_Status
HER2_Status
Met_Therapy_Line
Regimen
EndTrt_Days_From_Primary_Diagnosis
Treatment_Response
Progression
Reason_Tx_Ended
```

The file is read using:

```python
encoding="latin1"
```

This avoids Unicode decoding errors caused by Windows-style characters.

---

### 2. All-case CNV gene-state matrix

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_gene_CNV_state_matrix.tsv
```

Structure:

```text
rows    = genes
columns = longitudinal plasma samples
values  = gene-level CNV state
```

Example CNV states:

```text
NEUT
GAIN
AMP
HLAMP
HETD
HOMD
```

Although the filename contains `representative`, this file is intended to contain all available samples generated from the all-case annotated CNR folder.

---

### 3. All-case CNV sample annotation

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_CNV_sample_annotation.tsv
```

Important columns:

```text
sample_id
Patient_Identifier
draw_num
Tumor_Fraction
Tumor_Fraction_percent
Histology_group
Molecular_subtype_clean
CNV_burden_non_NEUT_genes
n_amp_or_hlamp_genes
n_hlamp_genes
```

`Tumor_Fraction_percent` is used in downstream summaries and blood draw tables.

---

### 4. All-case CNR manifest

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_CNR_manifest.tsv
```

This file records the source `.annotated.cnr` file and QC information for each sample.

## Output directory

The notebook writes outputs to:

```text
/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Clinical_data/clinical_CNV_extraction_TSVs_All_cases/All_cases_clinical_CNV_extraction_TSVs
```

## Main outputs

### Excel workbook

```text
All_cases_clinical_CNV_extraction.xlsx
```

Recommended Excel sheets:

```text
Case_summary
Blood_draws
Treatment_timeline
Diagnosis_receptor
CNV_sample_summary
CNV_gene_summary_key_genes
All_clinical_rows
CNV_rank_features_all
```

For all-case extraction, the full CNV gene summary should **not** be written to Excel because it can contain millions of rows.

---

### TSV files

Recommended TSV outputs:

```text
all_cases_case_summary.tsv
all_cases_blood_draws_with_treatment_CNV_TFx.tsv
all_cases_treatment_timeline.tsv
all_cases_diagnosis_receptor_tissue.tsv
all_cases_CNV_sample_summary.tsv
all_cases_CNV_gene_summary_key_genes.tsv
all_cases_all_clinical_rows.tsv
all_cases_CNV_rank_features.tsv
all_cases_unmapped_or_flagged_CNV_cases.tsv
```

The ZIP output is optional. If not needed, comment out the ZIP section.

## Important mapping logic

### Patient identifier extraction

The notebook extracts a 4-character `Patient_Identifier` from sample or barcode strings.

Examples:

```text
AUP-AQ0Z-TTP2-B -> AQ0Z
AUP-AQ0Z-TTM1-A -> AQ0Z
AUP-AMCZ-TTP1-B -> AMCZ
AMCZ-NB1        -> AMCZ
AQ30-NB1        -> AQ30
AO01-NB1        -> AO01
AM01-NB1        -> AM01
```

### Mapping to clinical patient_id

The notebook builds a clinical crosswalk using:

```text
clinical Barcode -> Patient_Identifier -> patient_id
```

Mapping status categories:

```text
Mapped to clinical patient_id
No clinical Barcode match
Conflicting clinical patient_id mapping
```

After mapping, the notebook prints:

```text
Mapping summary
Unmapped / flagged cases
Unmapped / flagged Patient_Identifier list
```

For the full all-case run, a typical result may look like:

```text
Mapped to clinical patient_id    262
No clinical Barcode match          8
```

This means the CNV matrix contains 270 cases, of which 262 map to clinical `patient_id`.

## Memory-safe CNV handling

The earlier 20-case script created a long table:

```python
cnv_long = cnv_matrix.melt(...)
```

For the all-case matrix, this can create approximately:

```text
24,881 genes × 784 samples = 19,506,704 gene-sample rows
```

Therefore, the all-case notebook avoids `cnv_long` and instead keeps the CNV data in wide format:

```python
state_mat = cnv_matrix.set_index("gene")[cnv_sample_cols].copy()
```

This is more memory-efficient and avoids unnecessary expansion.

## CNV sample summary

The notebook calculates sample-level CNV burden from the wide matrix.

Important columns:

```text
CNV_burden_non_NEUT_genes
Amplified_genes_HLAMP_or_AMP
Altered_genes_all_non_NEUT
```

Definitions:

```text
CNV_burden_non_NEUT_genes = number of genes with CNV state != NEUT
Amplified_genes_HLAMP_or_AMP = genes with AMP or HLAMP
```

For Excel, it is recommended to drop `Altered_genes_all_non_NEUT` because this column can contain very long strings for all-case analysis.

## Selected key CNV genes

The notebook uses a focused key-gene panel instead of generating a full all-gene summary for Excel.

Selected CNV genes:

```text
ERBB2
FGFR1
MYC
CCND1
MCL1
PIK3CA
ESR1
CDH1
TP53
AKT1
MAP3K1
FOXA1
```

These genes are used to create:

```text
all_cases_CNV_gene_summary_key_genes.tsv
```

and the Excel sheet:

```text
CNV_gene_summary_key_genes
```

This avoids creating a very large full gene summary table for all cases.

## CNV rank score

The notebook calculates a heuristic prioritization score for each case.

Features used:

```text
n_draws
max_TFx
range_TFx
mean_CNV_burden
range_CNV_burden
n_variable_genes
total_gene_state_changes
key_amp_calls
```

Formula:

```text
CNV_rank_score =
1.2 × normalized(n_draws)
+ 1.2 × normalized(max_TFx)
+ 1.0 × normalized(range_TFx)
+ 1.2 × normalized(mean_CNV_burden)
+ 1.5 × normalized(range_CNV_burden)
+ 1.5 × normalized(n_variable_genes)
+ 1.0 × normalized(total_gene_state_changes)
+ 1.3 × normalized(key_amp_calls)
```

This is not a clinical response score. It is a prioritization score for identifying cases with informative longitudinal CNV/TFx dynamics.

## Treatment timeline

The notebook extracts metastatic treatment rows using these clinical categories:

```text
Met_Chemotherapy
Met_Hormonal_Therapy
Met_Targeted_Therapy
Met_Immunotherapy
Met_Other_Therapy
```

For each treatment row, it records:

```text
treatment_category
treatment_start_day
treatment_end_day
treatment_line
regimen
treatment_response
progression
reason_tx_ended
treatment_label
treatment_interval
```

## Blood draw treatment context

For each CNV sample / blood draw, the notebook adds:

```text
Active_treatment_at_draw
Treatment_context_at_draw
Most_recent_treatment_before_draw
Most_recent_treatment_start_day
Most_recent_treatment_end_day
Most_recent_treatment_response_recorded
Most_recent_treatment_progression_recorded
Most_recent_treatment_reason_ended_recorded
Most_recent_completed_treatment_before_draw
Most_recent_completed_treatment_end_day
Days_since_most_recent_completed_treatment_end
Most_recent_completed_treatment_response_recorded
Most_recent_completed_treatment_progression_recorded
Most_recent_completed_treatment_reason_ended_recorded
Next_treatment_after_draw
Next_treatment_start_day
Days_until_next_treatment_start
Previous_blood_draw_day
Days_since_previous_blood_draw
Previous_TFx_percent
TFx_change_from_previous_draw
Previous_CNV_burden_non_NEUT_genes
CNV_burden_change_from_previous_draw
```

Treatment context labels:

```text
On-treatment
Post-treatment
Pre-treatment
Between-treatment-lines
No matched treatment context
No clinical blood draw date or patient_id mapping
```
## Treatment-context definition

Treatment context was assigned at the **blood-draw level** by comparing each cfDNA blood draw day with the patient’s metastatic treatment timeline.

Treatment start and end days were derived as:

```text
treatment_start_day = Days_From_Primary_Diagnosis
treatment_end_day   = EndTrt_Days_From_Primary_Diagnosis
```

Only metastatic treatment categories were used to define the treatment timeline, including chemotherapy, hormonal therapy, targeted therapy, immunotherapy, and other metastatic therapies.

---

## Blood-draw treatment-context categories

Each blood draw was assigned to one of the following treatment-context categories:

### On-treatment

A blood draw was classified as `On-treatment` when it occurred within a treatment interval:

```text
treatment_start_day <= blood_draw_day <= treatment_end_day
```

This required both treatment start and treatment end days to be available.

### Between-treatment-lines

A blood draw was classified as `Between-treatment-lines` when:

```text
a completed treatment occurred before the blood draw
AND
a subsequent treatment started after the blood draw
AND
the blood draw did not fall within an active treatment interval
```

Therefore, `Between-treatment-lines` represents an off-treatment interval after completion of a prior treatment line and before initiation of the next treatment line.

### Post-treatment

A blood draw was classified as `Post-treatment` when:

```text
a completed treatment occurred before the blood draw
AND
no subsequent treatment was recorded after the blood draw
AND
the blood draw did not fall within an active treatment interval
```

This means the blood draw occurred after a completed treatment in the available clinical timeline. It does not necessarily imply that no future treatment ever occurred outside the available records.

### Pre-treatment

A blood draw was classified as `Pre-treatment` when:

```text
no completed treatment occurred before the blood draw
AND
a subsequent treatment was recorded after the blood draw
AND
the blood draw did not fall within an active treatment interval
```

This represents a blood draw collected before the first matched treatment line in the available metastatic treatment timeline.

### No matched treatment context

A blood draw was classified as `No matched treatment context` when it could not be assigned to pre-treatment, on-treatment, between-treatment-lines, or post-treatment based on the available treatment timeline.

This can occur when:

```text
no active treatment overlapped the draw
AND
no completed prior treatment was identified
AND
no subsequent treatment was identified
```

or when the blood draw day or patient mapping was missing.

---

## Final treatment-context grouping for analysis

For downstream treatment-context association analysis, the original treatment-context labels were collapsed into three analysis groups:

```text
Pre-treatment
→ Pre-treatment

On-treatment
→ On-treatment

Post-treatment
Between-treatment-lines
→ Post-treatment/off-treatment

No matched treatment context
No clinical blood draw date or patient_id mapping
→ Excluded
```

`Post-treatment` and `Between-treatment-lines` were combined because both represent blood draws collected outside an active treatment interval after a prior treatment exposure. This combined group was labeled `Post-treatment/off-treatment`.

## Strict active-treatment definition

A treatment is considered active only when:

```text
treatment_start_day <= blood_draw_day <= treatment_end_day
```

Both treatment start and treatment end days must be present.

If the treatment end day is missing, that treatment is **not** assumed to be ongoing.

A blank `Active_treatment_at_draw` means:

```text
No treatment row with a clearly defined start/end interval overlapped this blood draw.
```

It does **not** necessarily mean the patient was untreated.

## Recommended Section 17 for Excel writing

For all-case analysis, Section 17 should avoid using the old full-gene object:

```python
candidate_cases_CNV_gene_summary
```

That object is not created in the memory-safe all-case version.

Use only the key-gene summary:

```python
candidate_cases_CNV_sample_summary_excel = candidate_cases_CNV_sample_summary.drop(
    columns=["Altered_genes_all_non_NEUT"],
    errors="ignore"
)

sheet_order = {
    "Case_summary": candidate_cases_case_summary,
    "Blood_draws": candidate_cases_blood_draws_with_treatment_CNV_TFx,
    "Treatment_timeline": candidate_cases_treatment_timeline,
    "Diagnosis_receptor": candidate_cases_diagnosis_receptor_tissue,
    "CNV_sample_summary": candidate_cases_CNV_sample_summary_excel,
    "CNV_gene_summary_key_genes": candidate_cases_CNV_gene_summary_key_genes,
    "All_clinical_rows": candidate_cases_all_clinical_rows_export,
    "CNV_rank_features_all": rank_features_all,
}
```

Do not write a full `CNV_gene_summary_non_NEUT` sheet for all cases unless it has been explicitly filtered to stay below Excel row limits.

## Important path note

Use one of these styles, not both together.

Correct relative-to-BASE_DIR style:

```python
cnv_matrix_file = BASE_DIR / "CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_gene_CNV_state_matrix.tsv"
```

Correct absolute path style:

```python
cnv_matrix_file = Path("/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/CNVKit/All_case_CNV_extraction/All_representative_patients_longitudinal_gene_CNV_state_matrix.tsv")
```

Avoid this mixed style:

```python
cnv_matrix_file = BASE_DIR / "/ix1/alee/LO_LAB/..."
```

Although Python will usually treat the second path as absolute, it is confusing and should be avoided.

## Recommended QC checks

After the notebook finishes, check:

```python
print("Cases present in CNV matrix:", len(cnv_patients_all))
print("CNV samples:", len(cnv_sample_cols))
print("Mapped cases:", int((candidate_patient_map["mapping_status"] == "Mapped to clinical patient_id").sum()))
print("Unmapped/flagged cases:")
print(candidate_patient_map[candidate_patient_map["mapping_status"] != "Mapped to clinical patient_id"])
```

Also check the output files:

```bash
cd /ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Clinical_data/clinical_CNV_extraction_TSVs_All_cases/All_cases_clinical_CNV_extraction_TSVs

ls -lh

wc -l all_cases_case_summary.tsv
wc -l all_cases_blood_draws_with_treatment_CNV_TFx.tsv
wc -l all_cases_CNV_sample_summary.tsv
wc -l all_cases_CNV_gene_summary_key_genes.tsv
```

## Expected interpretation of all-case run

If the notebook prints:

```text
Cases present in CNV matrix: 270
CNV samples: 784
Mapped cases: 262
```

then the all-case CNV matrix is being used correctly.

If it prints only around 20 cases, then the script is still pointing to the selected-case CNV matrix.

## Downstream use

The all-case workbook and TSVs can be used to prioritize cases for:

```text
longitudinal treatment-response interpretation
TFx trajectory analysis
CNV burden dynamics
key-gene amplification/deletion tracking
GRIFFIN/nucleosome profiling candidate selection
```

For GRIFFIN candidate selection, use this all-case clinical/CNV output together with available BAM/GC-bias YAML files.
