# Case Summary and CNV Rank Score Documentation

## Purpose

This document explains how the `Case_summary` sheet was created and how `CNV_rank_score` was calculated for prioritizing candidate cases for downstream treatment-history and GRIFFIN/nucleosome profiling.

The main idea is to start from longitudinal ichorCNA tumor fraction (`TFx`) trajectory categories, select representative patients, build a longitudinal gene-level CNV heatmap, and then rank those representative cases using CNV and TFx dynamics.

---

## Overall workflow

```text
Step 1: Use ichorCNA TFx trajectories
        ↓
Classify patients into trajectory categories
        ↓
Select representative patients from each category
        ↓
Make longitudinal CNV heatmap for those selected patients
        ↓
Calculate CNV_rank_score within this selected representative patient set
        ↓
Shortlist cases for treatment-history + GRIFFIN/nucleosome profiling
```

---

## Important distinction

The trajectory category and `CNV_rank_score` are related but not the same.

| Step | Purpose | Main information used |
|---|---|---|
| TFx trajectory category | Classify patients based on longitudinal tumor fraction behavior | baseline TFx, final TFx, min TFx, max TFx, TFx range, direction changes |
| CNV heatmap | Visualize gene-level CNV states across longitudinal samples | `.annotated.cnr` corrected calls |
| CNV_rank_score | Prioritize representative cases for downstream analysis | longitudinal TFx + gene-level CNV dynamics |
| Final shortlist | Choose cases for treatment-history and GRIFFIN profiling | CNV_rank_score + treatment interpretability + clinical context + GRIFFIN availability |

---

## Case_summary sheet

The `Case_summary` sheet has one row per candidate patient. It summarizes clinical context, longitudinal blood draws, tumor fraction, CNV burden, notable altered genes, treatment history, metastatic sites, and recommended downstream use.

### Case_summary columns

| Column | Meaning | Extracted from / calculated from |
|---|---|---|
| `Patient_Identifier` | Short patient/case ID, such as `ANGM`, `ANTU`, or `APZO`. This is the ID used in CNV sample names such as `ANGM-NB1`. | Extracted from CNV sample barcode/name and matched to clinical `Barcode`. |
| `patient_id` | Clinical patient ID, such as `588-21`. This is the key used in the clinical file to collect all rows for that patient. | Clinical file column `patient_id`. |
| `Site` | Institution/site, such as Mayo, Vanderbilt, UAB, etc. | Clinical file column `Site`. |
| `n_all_clinical_rows` | Total number of clinical rows found for that `patient_id`. This includes diagnosis, blood draw, treatment, accession, progression, survival, and other rows. | Count of all rows in the clinical file for the matched `patient_id`. |
| `n_blood_draw_rows_in_clinical` | Number of rows where `Category == "Blood_Draw"` for that patient. This tells how many blood draws are recorded clinically. | Clinical file rows where `Category == "Blood_Draw"`. |
| `n_CNV_samples_in_matrix` | Number of samples from that patient present in the gene-level CNV matrix. For example, `ANGM` has 6 CNV samples: `ANGM-NB1` to `ANGM-NB6`. | CNV state matrix and CNV sample annotation. |
| `blood_draw_days` | Lists blood draw barcode and day from primary diagnosis, for example `ANGM-NB1:4084; ANGM-NB2:4195`. | Clinical file rows where `Category == "Blood_Draw"`; uses `Barcode` and `Days_From_Primary_Diagnosis`. |
| `TFx_percent_min` | Minimum tumor fraction percentage among CNV samples for this patient. | CNV sample annotation column `Tumor_Fraction_percent`. |
| `TFx_percent_max` | Maximum tumor fraction percentage among CNV samples for this patient. | CNV sample annotation column `Tumor_Fraction_percent`. |
| `TFx_percent_mean` | Mean tumor fraction percentage across that patient’s CNV samples. | Calculated from CNV sample annotation. |
| `histology_from_ranked_or_annotation` | Case-level histology used for interpretation, such as Ductal or Lobular. | Primarily from CNV sample annotation / ranked candidate table; originally derived from clinical histology fields. |
| `subtype_from_ranked_or_annotation` | Molecular subtype used for interpretation, such as `HR+/HER2-`, `HR-/HER2-`, or `HR+/HER2+`. | Primarily from CNV sample annotation / ranked candidate table; originally derived from ER/PR/HER2 clinical annotation. |
| `CNV_rank_score` | Candidate ranking score created to prioritize cases. Higher score means more interesting based on longitudinal CNV/TFx features. | Derived score from ranked candidate table; not a raw clinical value. |
| `mean_CNV_burden` | Average number of altered genes per sample across that patient’s longitudinal samples. | Calculated from gene-level CNV matrix. |
| `max_CNV_burden` | Maximum number of altered genes seen in any one sample from that patient. | Calculated from gene-level CNV matrix. |
| `notable_genes_states` | Important genes and CNV states observed across timepoints, for example `MCL1(GAIN,AMP,HLAMP)`. | Calculated from gene-level CNV matrix. |
| `receptor_summary_from_clinical` | Text summary of receptor/histology information across primary and metastatic rows. Includes tissue site, day, ER, PR, HER2, and histology when available. | Clinical file columns: `Category`, `Days_From_Primary_Diagnosis`, `Tissue_Site`, `Histologic.type`, `ER_Status`, `PR_Status`, `HER2_Status`, etc. |
| `metastatic_sites_summary` | Summary of metastatic diagnosis sites, such as liver, bone, lymph node, brain, pleura, etc. | Clinical file rows where `Category == "Met_Diagnosis"`; uses `Tissue_Site`, `bx`, and `bx.details`. |
| `treatment_summary` | Compact summary of metastatic treatment lines, regimen, start/end days, response, and progression. | Clinical treatment rows such as `Met_Chemotherapy`, `Met_Hormonal_Therapy`, and `Met_Targeted_Therapy`; uses `Met_Therapy_Line`, `Regimen`, `Days_From_Primary_Diagnosis`, `EndTrt_Days_From_Primary_Diagnosis`, `Treatment_Response`, and `Progression`. |
| `survival_summary` | Survival status and survival day, for example `Dead at 1439 days` or `Alive with Disease at 2196 days`. | Clinical rows `Overall_Survival_days` and/or `date_of_death`. |
| `recommended_use` | Interpretation of why this case may be useful for downstream analysis. | Manually derived from TFx pattern, CNV burden, treatment interpretability, and clinical context. |

---

## Data sources used to create Case_summary

| Information type | Source |
|---|---|
| Clinical timeline, blood draws, treatment history, metastatic sites, receptor status, survival | Clinical file |
| TFx percentage, histology annotation, subtype annotation | CNV sample annotation file |
| Gene-level CNV state | Gene-level CNV state matrix |
| CNV burden and notable genes | Calculated from gene-level CNV state matrix |
| Candidate ranking | Calculated from TFx and CNV-derived features |

---

## CNV_rank_score

`CNV_rank_score` is a heuristic score designed to prioritize cases for downstream treatment-history and GRIFFIN/nucleosome profiling.

It was calculated for the representative patients selected from the longitudinal TFx trajectory categories.

The score answers this question:

```text
Which cases look most interesting from longitudinal TFx + gene-level CNV dynamics?
```

It does not directly answer:

```text
Which case has the best treatment-history story?
```

Treatment interpretability is evaluated after the CNV/TFx-based ranking.

---

## Features used in CNV_rank_score

A good candidate case should have:

- more longitudinal samples,
- higher tumor fraction,
- larger TFx change over time,
- higher CNV burden,
- more CNV change over time,
- more variable genes,
- more `AMP` or `HLAMP` events in key genes.

Each feature was normalized across the ranked representative candidate cases.

---

## Min-max normalization

For each feature:

```text
normalized value = (value − minimum value) / (maximum value − minimum value)
```

The minimum and maximum should be calculated across the cases included in the ranking.

For the main `CNV_rank_score`, this means the selected representative patients with CNV data, not only the final top six shortlisted cases.

### Example

If a patient has:

```text
n_draws = 6
```

and across all ranked representative cases:

```text
minimum n_draws = 3
maximum n_draws = 8
```

then:

```text
normalized(n_draws) = (6 − 3) / (8 − 3)
                    = 3 / 5
                    = 0.60
```

If we normalize only within a final 6-case shortlist and the values are:

```text
minimum n_draws = 5
maximum n_draws = 8
```

then:

```text
normalized(n_draws) = (6 − 5) / (8 − 5)
                    = 1 / 3
                    = 0.333
```

Therefore, normalized values depend on the case set used for normalization.

---

## CNV_rank_score formula

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

---

## CNV_rank_score feature definitions

| Feature | Meaning | Why it matters |
|---|---|---|
| `n_draws` | Number of longitudinal CNV samples for that patient. | More timepoints gives a better treatment/GRIFFIN story. |
| `max_TFx` | Highest tumor fraction across samples. | Higher TFx means stronger tumor-derived signal. |
| `range_TFx` | Maximum TFx minus minimum TFx. | Captures dynamic change over time. |
| `mean_CNV_burden` | Average number of non-NEUT genes per sample. | Captures overall CNV alteration burden. |
| `range_CNV_burden` | Maximum CNV burden minus minimum CNV burden. | Captures CNV evolution over time. |
| `n_variable_genes` | Number of genes whose CNV state changes across longitudinal samples. | Captures gene-level evolution. |
| `total_gene_state_changes` | Total number of CNV state changes across genes and timepoints. | Rewards cases with dynamic CNV trajectories. |
| `key_amp_calls` | Number of `AMP` or `HLAMP` calls in key genes. | Prioritizes strong focal/high-level events. |

---

## Interpretation of CNV_rank_score

Higher `CNV_rank_score` suggests that a case may be useful for downstream analysis because it has stronger longitudinal molecular dynamics.

However, the score should be treated as a screening tool only.

Final case selection should also consider:

- availability of treatment history,
- clarity of treatment line timing,
- availability of GRIFFIN/nucleosome profiling data,
- clinical interpretability,
- receptor and histology context,
- whether blood draws align with treatment changes or progression events.

---

## Suggested future use

Use `CNV_rank_score` to shortlist candidate cases, then inspect the sample-level blood draw table to determine whether the case has a clear clinical story.

A useful downstream case-study figure may include:

```text
Panel A: Treatment timeline
Panel B: TFx trajectory
Panel C: CNV burden trajectory
Panel D: Gene-level CNV states for selected genes
Panel E: GRIFFIN/nucleosome profile over the same timepoints
```

---

## Key caution

`CNV_rank_score` is not a formal biological score and should not be used alone to claim clinical or biological significance.

It is intended to prioritize cases for manual inspection and downstream integrative analysis.
