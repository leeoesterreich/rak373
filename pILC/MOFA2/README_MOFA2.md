# MOFA2 multi-omics integration of pleomorphic invasive lobular carcinoma (pILC)

## Overview

This workflow uses **Multi-Omics Factor Analysis (MOFA2)** to integrate four molecular data views from a pILC cohort: RNA-seq expression, WES-derived copy-number alterations (CNAs), WES somatic mutations, and MSK-IMPACT somatic mutations. It identifies latent factors describing variation across available molecular modalities and subsequently evaluates relationships with clinicopathologic variables and biological pathways.

The workflow is based on the following **three R Markdown scripts**, which must be executed in order:

1. `01_Prepare_MOFA2_input_GENCODEv47.Rmd` — normalize identifiers, filter features, construct four feature-by-sample matrices, and audit sample overlap.
2. `02_Train_MOFA2_model_GENCODEv47.Rmd` — construct and train the MOFA2 model, and export model/quality-control summaries.
3. `03_MOFA2_Downstream_GENCODEv47.Rmd` — attach clinical data; analyze variance explained, factor scores and weights; test clinical associations; and perform RNA-weight gene-set enrichment.

**Scope:** This README documents the *code and its configured behavior*, not validated result counts. The final sample count, feature count, number of retained factors, and significant findings must be read from the generated output tables after a successful run.

## 1. Molecular inputs and preprocessing

| View (exact MOFA2 name) | Input | Encoding / processing |
| --- | --- | --- |
| `RNA_VST_top5000var` | Stranded HTSeq RNA-seq, GENCODE v47; an existing VST-transformed matrix of the top 5,000 variable genes (`.rds`) | Continuous; gene-by-sample matrix. Clean sample IDs, average duplicated feature/sample names, exclude entirely nonfinite and zero-variance features. If >5,000 rows, retain the 5,000 most variable. **VST is performed upstream**, not by these scripts. |
| `WES_CNV_numeric` | Gene-level WES CNV calls from PCGR (`.tsv`) | OncoKB oncogene/TSG genes only. Homozygous deletion = −2; heterozygous deletion/loss = −1; neutral = 0; gain = +1; amplification = +2. For duplicate gene–sample calls retain the largest absolute value; discard zero-variance features. Unrepresented gene–sample entries in the constructed wide matrix are filled with 0. |
| `WES_mutation_recurrent` | Combined WES MAF (`.maf`) | OncoKB oncogene/TSG genes and selected mutation classes only. If `FILTER` exists, keep PASS, missing, or blank values. Binary gene-by-sample encoding (0/1). Retain genes mutated in **at least two distinct samples**; discard zero-variance features. |
| `IMPACT_mutation_OncoKB` | Filtered MSK-IMPACT mutation table (`.tsv`) | OncoKB oncogene/TSG genes and selected mutation classes; binary 0/1. **No recurrence threshold** because this is a targeted panel. Preserve samples present in the input, including those with no retained mutations. |

### Required input files (paths currently hard-coded in script 01)

```text
RNA_VST_FILE
/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_VST_matrix_top5000_variable_genes.rds

WES_MAF_FILE
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Tier_1_common_all_samples.combined.maf

WES_CNV_FILE
/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/CNV_PCGR/Plot/Figures/WES_PCGR_CNV_onecall_per_gene_sample.tsv

IMPACT_FILE
/ix1/alee/LO_LAB/Personal/Rahul/pILC/MSK_IMPACT/pILC_IMPACT_filtered.tsv

ONCOKB_GENE_FILE
/ix1/alee/LO_LAB/Personal/Rahul/pILC/MOFA2/cancerGeneList.tsv
```

The OncoKB gene table is expected to contain gene identifiers (for example, `Hugo Symbol`) and a gene-type field (for example, `Gene Type`); rows whose gene type contains `ONCOGENE` or `TSG` are selected. This is a **gene selection list**, not a separate fifth MOFA2 view.

**Mutation classes retained:** `Frame_Shift_Del`, `Frame_Shift_Ins`, `In_Frame_Del`, `In_Frame_Ins`, `Missense_Mutation`, `Nonsense_Mutation`, `Nonstop_Mutation`, `Splice_Site`, `Translation_Start_Site`.

### Identifier matching

The helper `clean_sample_id()` uppercases and trims names, strips trailing sequencing suffixes such as `_S25`, changes dots to hyphens (`TP03.117` → `TP03-117`), removes selected tumor/normal endings, and removes leading zeros in the numeric TP suffix. Gene names in variant/CNV tables are uppercased and trimmed. **Check the exported raw-to-clean mapping** because collapsing distinct original IDs into a single ID can affect sample membership. The script uses the union of samples across views; it does not restrict training to samples measured in all four modalities.

The first script checks that **at least one sample has two or more views**. The second script explicitly requires all four nonempty views, even though the first script would omit empty ones from its saved list.

## 2. Dependencies and environment

- R with `rmarkdown`, `knitr`, `data.table`, `ggplot2`, `MOFA2`, `fgsea`, and `msigdbr`.
- `AnnotationDbi` and `org.Hs.eg.db` **when RNA feature identifiers are Ensembl IDs** and pathway analysis needs symbol mapping.
- A functional MOFA2 Python/backend environment available to `run_mofa(..., use_basilisk = TRUE)`; ensure the MOFA2/basilisk setup works on the execution host.
- Appropriate R graphics devices for PDF and 300-dpi PNG outputs, plus write access to the output directory.

The scripts prepend this site-specific package library:

```r
USER_LIB <- "/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(USER_LIB, .libPaths()))
```

Change `USER_LIB` to an accessible path or remove it when using another R installation. Package versions are **not pinned** in the uploaded scripts; save `sessionInfo()` from the rendered reports for provenance.

## 3. Run order

Edit the paths in the **`paths` / `paths-and-settings` chunks** before running. All scripts use:

```text
OUT_DIR = /ix1/alee/LO_LAB/Personal/Rahul/pILC/MOFA2
```

Render the R Markdown files from R (or run their chunks sequentially in RStudio):

```r
rmarkdown::render("01_Prepare_MOFA2_input_GENCODEv47.Rmd")
rmarkdown::render("02_Train_MOFA2_model_GENCODEv47.Rmd")
rmarkdown::render("03_MOFA2_Downstream_GENCODEv47.Rmd")
```

Run the next script **only if** the previous one finishes without errors. Script 02 reads the prepared `.rds` file from script 01. Script 03 reads the HDF5 model and sample-view-presence summary produced by script 02.

### Model configuration (script 02)

| Option | Value in script |
| --- | --- |
| Starting number of factors | `MAX_FACTORS = 10` |
| Factor dropping | `DROP_THRESHOLD = 0.02` (configured 2% threshold) |
| Random seed | `MODEL_SEED = 42` |
| Iteration cap | `maxiter = 1000` |
| Convergence mode | `medium` |
| Scale views | `FALSE` |
| RNA VST likelihood | Gaussian |
| WES CNV likelihood | Gaussian |
| WES mutation likelihood | Bernoulli |
| MSK-IMPACT mutation likelihood | Bernoulli |
| Backend call | `run_mofa(..., use_basilisk = TRUE)` |
| Model tag | `K10_drop02_seed42_GENCODEv47` |

MOFA2 can **retain fewer than 10 factors** following training; inspect the actual retained factor names in `MOFA2_training_summary_*.tsv`. A Bernoulli likelihood is used for mutation occurrence, whereas integer-valued copy-number categories are treated as Gaussian in this implementation.

## 4. Downstream analysis (script 03)

The script requires the trained `.hdf5` file, the sample-view-presence table, and:

```text
/ix1/alee/LO_LAB/Personal/Rahul/pILC/MOFA2/ClinicalCorrelation.csv
```

### Clinical metadata

Clinical sample IDs are standardized using the same helper. The code reads age at diagnosis, overall tumor cellularity, pleomorphic invasive-component cellularity, tumor-infiltrating lymphocytes (TILs), grade, ER, PR, and HER2. ER/PR/HER2 are converted to positive/negative categories; hormone-receptor status is defined as **HR+ when either ER or PR is positive**, HR− when both are negative. Combined HR/HER2 categories are `HR+/HER2-`, `HR+/HER2+`, `Triple-negative`, and `HR-/HER2+`. An entry of `<N` in a numeric clinical field is encoded as **N/2**. Only the first entry for duplicated cleaned clinical sample IDs is kept.

### Analyses produced

1. **Variance explained** — total per view; contribution of each factor to each view; factor summaries and dominant contributing view. Automatically focus on factors with at least **5% variance explained in any view** (`FOCUS_R2_CUTOFF = 5`); `FOCUS_FACTORS_OVERRIDE` can manually select factor numbers.
2. **Latent factor visualization** — long/wide factor-score tables and plots colored by available view-completeness and clinical annotations, including pairwise factor comparisons.
3. **Feature weights** — all feature weights; top 25 by absolute weight per view/factor; positive/negative weight tables; interpretation tables and factor-specific weight plots. `HIGHLIGHT_GENES` optionally marks chosen RNA genes in bar plots.
4. **Covariation patterns** — data heatmaps and scatterplots for selected factors/views, including clinical coloring where available.
5. **Clinical association tests** — factor scores versus age, cellularity, pleomorphism and TILs by **Spearman correlation**; two-category clinical fields by **Wilcoxon rank-sum**; and fields with >2 categories by **Kruskal–Wallis**. Benjamini–Hochberg (BH) FDR is calculated across all factor × clinical-feature tests; results and p-value/FDR heatmaps are exported.
6. **Ranked gene-set enrichment** — for each focus factor, rank RNA genes by **signed MOFA2 feature weight** (not differential-expression statistics); run `fgseaMultilevel` against MSigDB **Hallmark**, **Reactome**, and **GO Biological Process**. If needed, map version-stripped Ensembl IDs to human gene symbols using `org.Hs.eg.db`; when multiple RNA features map to the same symbol, keep the largest absolute weight. `minSize = 10`, `maxSize = 500`, significant `padj < 0.05`. Plot up to 40 significant pathways per focus factor (ranked by absolute NES).

**Interpretation:** Factor-weight GSEA assesses pathways associated with latent molecular axes; it does **not** directly test pathways differentially expressed between clinical groups. The sign of a latent factor is arbitrary, so interpret positive/negative weights relative to the factor orientation in this saved model.

## 5. Generated output structure

```text
MOFA2/
├── MOFA2_RDS/
│   ├── MOFA2_prepared_input_view_matrices_GENCODEv47.rds
│   ├── MOFA2_training_input_long_format_<MODEL_TAG>.rds
│   ├── MOFA2_pILC_model_<MODEL_TAG>_trained.rds
│   └── MOFA2_pILC_model_<MODEL_TAG>_with_clinical_metadata.rds
├── MOFA2_Tables/
│   ├── RNA_sample_ID_mapping_GENCODEv47.tsv
│   ├── MOFA2_prepared_view_summary_GENCODEv47.tsv
│   ├── MOFA2_prepared_sample_view_presence_GENCODEv47.tsv
│   ├── MOFA2_training_input_view_summary_<MODEL_TAG>.tsv
│   ├── MOFA2_training_sample_view_presence_<MODEL_TAG>.tsv
│   ├── MOFA2_training_input_long_format_<MODEL_TAG>.tsv
│   └── MOFA2_training_summary_<MODEL_TAG>.tsv
├── MOFA2_Model/
│   └── MOFA2_pILC_model_<MODEL_TAG>.hdf5
├── MOFA2_Plots/
│   ├── MOFA2_training_data_overview_<MODEL_TAG>.pdf
│   └── MOFA2_trained_model_data_overview_<MODEL_TAG>.pdf
├── MOFA2_Downstream/
│   ├── ClinicalCorrelation_clean_for_MOFA2.csv
│   ├── MOFA2_sample_metadata_with_clinical_annotations.tsv
│   ├── MOFA2_metadata_summary.tsv
│   ├── MOFA2_total_variance_explained_by_view.tsv
│   ├── MOFA2_variance_explained_per_factor_per_view.tsv
│   ├── MOFA2_dominant_view_per_factor.tsv
│   ├── MOFA2_factor_variance_summary.tsv
│   ├── MOFA2_factor_scores_wide_with_metadata.tsv
│   ├── MOFA2_feature_weights_long.tsv
│   ├── MOFA2_top25_weighted_features_by_view_factor.tsv
│   ├── MOFA2_factor_clinical_association_tests.tsv
│   ├── MOFA2_factor_clinical_pvalue_matrix.tsv
│   ├── MOFA2_factor_clinical_FDR_matrix.tsv
│   ├── MOFA2_downstream_run_summary.tsv
│   └── [other factor/weight tables and PDF/PNG plots]
└── MOFA2_Pathway_GSEA/
    ├── Focus_factors_all_RNA_weights_mapped.tsv
    ├── Focus_factors_ranked_GSEA_all_results.tsv
    ├── Focus_factors_ranked_GSEA_FDR_significant.tsv
    ├── Focus_factors_Hallmark_GSEA_FDR_significant.tsv
    └── Factor*_significant_GSEA.pdf / .png (when significant)
```

Here `<MODEL_TAG>` is `K10_drop02_seed42_GENCODEv47` for the uploaded version. Some plots are conditional on available clinical variables, views, focus factors, or significant pathways; therefore not every optional plot is guaranteed.

## 6. Notes and limitations

- This is an **unsupervised** latent-factor integration workflow; clinical outcomes are assessed **after** MOFA2 training rather than supplied as training labels.
- Molecular views have different sample coverage, and model training uses all sample IDs represented by the input views. Missing *whole views* are supported conceptually by the representation; the long-format table only contains observed finite values for represented sample–view combinations.
- `TOP_N_RNA_GENES = 5000` limits the RNA feature count only when more than 5,000 genes are loaded; it does not guarantee exactly 5,000 pass subsequent checks.
- WES mutation recurrence and targeted-panel mutation filtering are deliberately different.
- The WES CNV view encodes discrete copy-number categories as numeric data and models them with a Gaussian likelihood; consider this assumption when interpreting factor loadings.
- The scripts do not implement train/test splitting or external factor validation, and do not themselves establish causal links between factors and clinical variables.
- The uploaded R Markdown code includes absolute project-specific paths; adapt them before transferring the workflow to another environment. No job-submission/SLURM script was provided with these three files.
