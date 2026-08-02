# Griffin Breast-Lineage ER Analysis

This repository contains two Jupyter notebooks for evaluating whether cfDNA nucleosome features at the breast-lineage transcription-factor binding sites **ESR1, FOXA1, and GATA3** are associated with estrogen-receptor status and improve classification of metastatic ER phenotype in the AURORA breast-cancer cohort.

## Notebooks

### 1. `Griffin_Breast_Lineage_ER_1.ipynb`

Performs association analyses between Griffin breast-lineage features and ER status.

The notebook:

- selects one first cfDNA blood draw per patient using the lowest available `NB` number;
- uses Griffin measurements for `ESR1`, `FOXA1`, and `GATA3`;
- analyzes:
  - `mean_coverage`;
  - `central_coverage`;
  - `amplitude`;
- evaluates primary-tumor ER status and metastatic ER status in separate models;
- constructs standardized breast-lineage module scores;
- fits adjusted regression models with HC3-robust standard errors;
- includes TF-specific and module-level analyses;
- applies Benjamini-Hochberg false-discovery-rate correction;
- generates tabular results and forest plots.

Primary and metastatic ER status are analyzed separately and are not entered together in the same association model.

### 2. `Griffin_Breast_Lineage_ER_2.ipynb`

Performs internally cross-validated classification of metastatic ER status.

The primary question is whether coordinated `ESR1`/`FOXA1`/`GATA3` Griffin signals add information beyond historical primary ER status and plasma tumor fraction.

The notebook compares four models:

1. TFx only
2. Primary ER + TFx
3. Primary ER + TFx + three breast-lineage module scores
4. Primary ER + TFx + nine individual Griffin values

The three module scores correspond to:

- mean coverage;
- central coverage;
- amplitude.

For each module, the `ESR1`, `FOXA1`, and `GATA3` measurements are median-imputed and standardized within each training fold, averaged, and standardized again.

Validation uses repeated stratified cross-validation:

- 5 outer folds;
- 10 repeats;
- elastic-net logistic regression;
- fold-specific imputation and scaling;
- averaged out-of-fold predictions for each patient;
- patient-level bootstrap confidence intervals;
- paired bootstrap comparisons between models.

This is an internal-validation analysis and should not be interpreted as external clinical validation.

## Clinical outcome definitions

### Primary ER status

Primary ER status is obtained from:

```text
Metasheet_Analysis.with_Histology_Receptor_status.First_Primary.xlsx
```

### Metastatic ER status

Metastatic ER status is derived from clinical records with:

```text
Category == "Met_Diagnosis"
```

Metastatic records are searched chronologically. Blank or `Unknown` ER entries are skipped, and the earliest informative `Positive` or `Negative` result is selected.

## Required input files

The notebooks expect the following AURORA files:

```text
Clinical_data/
├── AUP_ClinicalData_20260609.xlsx
├── Metasheet_Analysis.with_Histology_Receptor_status.First_Primary.xlsx
└── Metasheet_Analysis.with_Histology_Receptor_status.Earliest_Met_Diagnosis.xlsx
```

They also require the Griffin summary workbook containing sample-level TFBS measurements:

```text
Griffin_22TF_GC_corrected_summary_by_module.xlsx
```

The classification notebook reads the `All_TFBS_Long` sheet from this workbook.

Before running the notebooks, confirm that all input paths match the local project structure.

## Main Griffin variables

For each of the three breast-lineage TFBSs, the analyses use:

```text
ESR1_mean_coverage
ESR1_central_coverage
ESR1_amplitude

FOXA1_mean_coverage
FOXA1_central_coverage
FOXA1_amplitude

GATA3_mean_coverage
GATA3_central_coverage
GATA3_amplitude
```

## Association-analysis outputs

The first notebook writes results to:

```text
Griffin/Downstream/Griffin_features_Primary_vs_Metastatic_ER/
```

Major outputs include:

```text
Griffin_BreastLineage_Primary_vs_Metastatic_ER_Module_Analysis.xlsx
Griffin_BreastLineage_Module_Model_Results.csv
Griffin_BreastLineage_Module_Analysis_Cohort.csv
Griffin_BreastLineage_Primary_vs_Metastatic_ER_Module_Forest.png
Griffin_BreastLineage_Primary_vs_Metastatic_ER_Module_Forest.pdf
```

The notebook may also create receptor-source-specific Excel workbooks named:

```text
Griffin_BreastLineage_<RECEPTOR_SOURCE>_ER_Analysis.xlsx
```

where `<RECEPTOR_SOURCE>` is either:

```text
First_Primary
Earliest_Met_Diagnosis
```

## Metastatic ER classification outputs

The second notebook writes outputs under its configured ROC-analysis directory.

Major outputs include:

```text
AURORA_Griffin_Metastatic_ER_ROC_Cohort.csv
AURORA_Griffin_Metastatic_ER_ROC_Metrics.csv
AURORA_Griffin_Metastatic_ER_ROC_Paired_Deltas.csv
AURORA_Griffin_Metastatic_ER_ROC_OOF_Predictions.csv
AURORA_Griffin_Metastatic_ER_ROC_Hyperparameters.csv
AURORA_Griffin_Metastatic_ER_ROC_Subgroups.csv
AURORA_Griffin_ER_Concordance.csv
AURORA_Griffin_Metastatic_ER_ROC_Methods.csv
AURORA_Griffin_Metastatic_ER_ROC_Run_Info.json
AURORA_Griffin_Metastatic_ER_ROC_Module.png
AURORA_Griffin_Metastatic_ER_ROC_Module.pdf
AURORA_Griffin_Metastatic_ER_ROC_9Feature_Sensitivity.png
AURORA_Griffin_Metastatic_ER_ROC_9Feature_Sensitivity.pdf
AURORA_Griffin_Metastatic_ER_PR_Module.png
AURORA_Griffin_Metastatic_ER_PR_Module.pdf
```

## Software requirements

The notebooks use Python 3 and the following packages:

```text
pandas
numpy
scipy
statsmodels
scikit-learn
matplotlib
openpyxl
```

Install the required packages with:

```bash
pip install pandas numpy scipy statsmodels scikit-learn matplotlib openpyxl
```

## Recommended execution order

Run the notebooks in this order:

```text
1. Griffin_Breast_Lineage_ER_1.ipynb
2. Griffin_Breast_Lineage_ER_2.ipynb
```

The first notebook evaluates statistical associations between breast-lineage Griffin features and ER status. The second notebook evaluates their incremental classification performance for metastatic ER phenotype.

## Suggested repository structure

```text
Griffin_Breast_Lineage_ER/
├── README.md
├── Griffin_Breast_Lineage_ER_1.ipynb
├── Griffin_Breast_Lineage_ER_2.ipynb
└── .gitignore
```

Remove `-checkpoint` from notebook filenames before committing them to GitHub.

## Data privacy

Clinical spreadsheets, patient-level Griffin results, analysis cohorts, and out-of-fold prediction files should not be committed to a public repository.

A suitable `.gitignore` may include:

```gitignore
.ipynb_checkpoints/
*.xlsx
*.csv
*.json
*.png
*.pdf
__pycache__/
```

## Interpretation notes

- The association notebook uses one first cfDNA draw per patient and is not a longitudinal model.
- Primary and metastatic ER associations are evaluated separately.
- Cross-validation reduces optimism but does not replace validation in an independent cohort.
- AUROC improvement does not by itself establish prediction of receptor conversion, particularly when primary and metastatic ER statuses are highly concordant.
