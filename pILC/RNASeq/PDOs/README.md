# PDO RNA-seq analysis (BO-046, BO-114, BO-189, BO-268)

## Overview
This folder holds the RNA-seq analysis of four patient-derived organoids (PDOs): BO-046, BO-114, BO-189 and BO-268. The analyses repeat what was done for the 31 pILC tumors and place the PDOs within the tumors (RNA clusters k = 3, 50-run NMF, GenePattern numbering).

Root folder: `/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/PDO/`

Important limits for every result:
- **n = 4, one sample per PDO, no replicates.** All PDO-vs-tumor comparisons are descriptive or exploratory. Where a test is reported (Wilcoxon, BH), it is labelled exploratory.
- PDOs differ from tumors in **culture conditions** and have **no stromal or immune cells**. Differences in immune, stromal or fat-related programs are therefore expected and are not tumor-intrinsic findings.
- PDOs have no RNA cluster and no clinical data, so no cluster or clinical tests are run.

## Source files
Location: `/ix1/alee/LO_LAB/Personal/Michael_Gorry/Org_RNAseq_Data/`

| File | Use |
|------|-----|
| `Org_RAWcounts__ENSG_63140.csv` | Raw counts with Ensembl IDs (used for extraction) |
| `colData_organoid_metadata.xlsx` | Sample metadata (checked only; printed to the console) |

Other inputs used from the tumor analysis (read only, not changed):
- Tumor raw counts with gene names: `4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/` (`htseq_raw_counts_with_gene_names.csv`, `htseq_VST_matrix_all_expression_filtered_genes.rds`)
- Tumor TPM: `4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt`
- Tumor clusters: `NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv` (GenePattern numbering, columns `sample`, `cluster`)
- NMF weights: `NMF/R_NMFConsensus_50runs/NMF_best_run_k3_W.csv`
- Reference: GENCODE v47 basic annotation (`gencode.v47.basic.annotation.gtf`), MSigDB Hallmark (`Pathway_analysis/Dataset/h.all.v2025.1.Hs.symbols.gmt`), CollecTRI (`Decoupler/collectri_human.tsv`), PBCA InstaPrism reference (`Immune_deconvolution_BC_Atlas/PBCA_InstaPrism_reference.Rds`)

## Run order

| Step | Script | Needs |
|------|--------|-------|
| 1 | `Extract_PDO_raw_counts.Rmd` | Source count file |
| 2 | `TPM_PDO.R` | Step 1 |
| 3 | `PDO_counts_to_VST_GENCODEv47.R` | Step 1 |
| 4 | `PCA_PDO.R` | Step 3 |
| 5 | `immunedeconv_ESTIMATE_PDO.R` | Steps 2 and 3 (tumor ESTIMATE scores optional, for the comparison figure) |
| 6 | `InstaPrism_PBCA_PDO.R` | Step 2 |
| 7 | `InstaPrism_PBCA_PDO_downstream.R` | Steps 5 (optional) and 6, and the saved tumor InstaPrism results |
| 8 | `InstaPrism_PBCA_PDO_cancer_epithelial_genes.R` | Step 2 and the tumor cancer-epithelial gene ranking (step 6 fractions optional, for the annotation strip) |
| 9 | `PDO_decoupler_TF_activity.R` | Step 3, and the tumor decoupleR script (it creates the CollecTRI file) |
| 10 | `PDO_decoupler_volcano_PDO_vs_tumors.R` | Step 9 |
| 11 | `PDO_decoupler_TF_boxplot.R` | Step 9 |
| 12 | `PDO_PAM50.R` | Step 3 (tumor PAM50 calls optional, for the concordance check) |
| 13 | `PDO_ssGSEA_Hallmark.R` | Step 3, including its PDO-only VST matrix (saved tumor ssGSEA scores optional, for a check) |
| 14 | `PDO_NMF_projection.R` | Step 3 and the NMF W file |
| 15 | `PDO_cGAS_STING.R` | Step 1 and the NMF cluster file (saved tumor VST and metadata optional, for checks) |
| 16 | `PDO_cGAS_STING_cancer_cell.R` | Steps 2 and 6 and the saved tumor InstaPrism result (tumor cancer-epithelial ranking file optional) |

Steps 4-16 can be run in any order after their inputs exist. All scripts use seed 42 and save PNG and PDF figures with a matching CSV table.

## Files in this folder (raw counts, step 1)

| File | Description |
|------|-------------|
| `PDO_raw_counts.csv` | Raw integer counts, genes (Ensembl IDs with version) x 4 samples |
| `PDO_sample_to_bam_map.csv` | BAM file from which each sample's counts came |
| `PDO_library_sizes.csv` | Total counts and number of genes with counts per sample |
| `Extract_PDO_raw_counts.Rmd` | Script that creates the three files above |

### How the counts were extracted
1. Read the raw count file as text.
2. Match each sample name (BO-046, BO-114, BO-189, BO-268) to exactly one column. Names are compared without case, dots, hyphens or underscores.
3. Remove the second header row (`Geneid` + BAM paths, featureCounts style). Without this step, R reads all columns as text. The BAM paths are saved in `PDO_sample_to_bam_map.csv`.
4. Convert to integers and check that no value was lost, all values are whole numbers and not negative, and gene IDs are unique.
5. Save the matrix, the BAM map and the library sizes.

### BAM files

| Sample | BAM name in the count file |
|--------|----------------------------|
| BO-046 | `BO046_S11_Aligned.sortedByCoord.out.bam` |
| BO-114 | `BO114_S9_Aligned.sortedByCoord.out.bam` |
| BO-189 | `BO189_S13_Aligned.sortedByCoord.out.bam` |
| BO-268 | `BO268_S36_BAligned.sortedByCoord.out.bam` |

### Library sizes

| Sample | Total counts | Genes with counts |
|--------|--------------|-------------------|
| BO-046 | 33,990,926 | 28,206 |
| BO-114 | 38,179,906 | 29,705 |
| BO-189 | 38,853,523 | 29,971 |
| BO-268 | 24,265,130 | 31,984 |

## Analyses and outputs

### TPM (`TPM_PDO.R`) → `Max_TPM/`
Gene length is the union of exons (GENCODE v47 basic). Gene symbols come from the GTF, and one Ensembl ID is kept per symbol (highest max TPM).
- `tpm_matrix_PDO.txt`, `tpm_matrix_PDO_log2TPM.txt` (log2(TPM + 1)), `HGNC_to_Ensembl_mapping.txt`, and one `<sample>_TPM.tsv` per PDO.

### VST (`PDO_counts_to_VST_GENCODEv47.R`) → `VST_data_GENCODEv47/`
Counts are labelled with gene names (the duplicated symbols are written as `SYMBOL|ENSG`). Low-expression filter: count >= 5 in at least 10% of samples (1 of 4 samples here). VST with `blind = TRUE` (design `~ 1`). Top 5000 variable genes are also saved.
- Labelled counts (**used by the later scripts**): `htseq_raw_counts_with_gene_names.csv` / `.rds`; filtered counts `htseq_expression_filtered_counts_gene_names.csv` / `.rds`
- VST: `htseq_VST_gene_names_all_expression_filtered_genes.csv`, `htseq_VST_matrix_all_expression_filtered_genes.rds`, `DESeq2_VST_object.rds`, `htseq_VST_gene_names_top5000_variable_genes.csv`, `htseq_VST_matrix_top5000_variable_genes.rds`, `htseq_VST_gene_variance_table.csv`, `VST_sample_names.csv`, `vst_transformed_data.gct`
- Other: `gene_annotation_from_GTF.tsv`, `PDO_assigned_reads_per_sample.csv`

### PCA of the PDOs (`PCA_PDO.R`) → `Dimensionality_plot/`
PDOs only: VST, top 2000 variable genes, filter >= 10 counts in at least 2 samples. With 4 samples there are only 3 informative PCs and no batch variable, so no batch statistics.
- `PCA_PDO.png` / `.pdf`, `PCA_PDO_coordinates.csv`, `PCA_PDO_variance.csv`

### ESTIMATE (`immunedeconv_ESTIMATE_PDO.R`) → `immunedeconv/TPM/` and `immunedeconv/VST/`
ESTIMATE (immunedeconv) on the PDO TPM matrix (primary) and the PDO VST matrix (sensitivity check). No group statistics (n = 4).
- `ESTIMATE_scores_PDO.csv`, `ESTIMATE_scores_PDO.png` / `.pdf`, `ESTIMATE_scores_PDO_vs_tumors.png` / `.pdf` (tumors by cluster from the saved tumor ESTIMATE scores)
- `immunedeconv/VST/ESTIMATE_TPM_vs_VST_correlation_PDO.csv` (TPM versus VST agreement)

### InstaPrism with the PBCA reference (`InstaPrism_PBCA_PDO.R`) → `InstaPrism_PBCA_reference/`
Same PBCA reference as the tumors, run on the PDO TPM matrix. No cluster strip. The script warns if the Cancer Epithelial fraction is below 0.8.
- `PDO_InstaPrism_PBCA_cell_fractions.csv`, `PDO_InstaPrism_PBCA_fraction_with_summary.csv`, `PDO_InstaPrism_PBCA_deconv_result.rds`, `PBCA_InstaPrism_input_summary_PDO.csv`, `R_sessionInfo_PDO.txt`
- Figures: `PDO_InstaPrism_PBCA_stacked_barplot`, `PDO_InstaPrism_PBCA_fraction_heatmap` (PDF), `PDO_InstaPrism_PBCA_total_immune_fraction`

### InstaPrism downstream (`InstaPrism_PBCA_PDO_downstream.R`) → `InstaPrism_PBCA_reference/Downstream/`
Descriptive: PDO key fractions against the tumors (by cluster), cell states, compartments.
- `PDO_InstaPrism_PBCA_key_fractions_per_PDO.csv`, `PDO_InstaPrism_PBCA_key_fractions_vs_tumors.csv` (includes the PDO percentile within the tumors), `PDO_InstaPrism_PBCA_cell_states.csv`, `PDO_InstaPrism_PBCA_compartment_fractions.csv`
- Figures: `PDO_InstaPrism_PBCA_key_fractions_PDO_vs_tumors`, `PDO_InstaPrism_PBCA_cell_states_top`

### Cancer-epithelial genes (`InstaPrism_PBCA_PDO_cancer_epithelial_genes.R`) → `InstaPrism_PBCA_reference/Downstream/Cancer_epithelial_genes/`
Uses the tumor ranking of cancer-epithelial genes (`pILC_InstaPrism_PBCA_cancer_epithelial_gene_ranking.csv`). Genes are z-scored with the **tumor** mean and SD, and the heatmap shows Cluster 1-3 and the PDOs.
- `PDO_cancer_epithelial_signature_scores.csv`, `PDO_cancer_epithelial_genes_zscores.csv`, `PDO_cancer_epithelial_top_genes_heatmap` (PNG + PDF)

### decoupleR TF activity → `Decoupler/`
CollecTRI network, `ulm` method. Genes are z-scored with the **tumor** mean and SD (not across the 4 PDOs), so PDO values are "SDs from the tumor average". Tumor TF scores are recomputed on the same genes, so they can differ slightly from the main tumor run. The only test is PDOs versus all tumors (Wilcoxon, BH), exploratory.
- `PDO_decoupler_TF_activity.R`: `TF_activity_per_sample_PDO.csv`, `TF_activity_per_sample_tumors_common_genes.csv`, `TF_activity_PDO_ulm_pvalues.csv`, `TF_PDO_vs_tumors_comparison.csv`, `TF_PDO_vs_tumors_stats.csv`, `TF_activity_heatmap_PDO_vs_tumors`, `TF_activity_PDO_top_TFs_barplot`
- `PDO_decoupler_volcano_PDO_vs_tumors.R`: `TF_activity_volcano_PDO_vs_tumors` (PNG + PDF) with `_counts.csv`, `_labeled_TFs.csv` and `_significant_TFs.csv` tables
- `PDO_decoupler_TF_boxplot.R`: `TF_activity_selected_TFs_boxplot_PDO_vs_tumors` (PNG + PDF) and `_stats.csv`. The TFs to plot are set at the top of the script (AR, ESR1, FOXA1, HDAC1, CIITA, BARX1 by default).
- If `select()` fails with a "tbl_df" error, restart R. The scripts already contain `select <- dplyr::select` to avoid the masking by AnnotationDbi.

### PAM50 (`PDO_PAM50.R`) → `PAM50/`
genefu PAM50 (Parker 2009 centroids). One **joint** VST is built from the raw counts of the 31 tumors and the 4 PDOs. The primary PDO call comes from genefu run on all 35 samples together. Because calls can depend on the cohort, the script also runs a tumors-only and a PDO-only version and compares them with the original saved tumor calls.
- `PDO_PAM50_subtype_calls.csv`, `PAM50_calls_tumors_and_PDO_joint_run.csv`, `PAM50_tumor_call_concordance.csv`, `PAM50_genes_used.csv`, `joint_tumor_PDO_VST_matrix.rds`
- Figures: `PDO_PAM50_subtype_by_group`, `PDO_PAM50_subtype_probability_by_sample`, `PDO_PAM50_centroid_correlation_by_sample`

### Hallmark ssGSEA (`PDO_ssGSEA_Hallmark.R`) → `Pathway_analysis/ssGSEA/`
**One joint run on the 31 tumors and the 4 PDOs** (same call as the tumor script: `GSVA::ssgseaParam(..., normalize = TRUE)` and `gsva()`; GSVA >= 1.50). The joint VST is built from the raw counts of both sets, so PDO and tumor scores are on one scale. This replaces the earlier PDO-only run.
- Scores: `ssGSEA_Hallmark_scores_PDO.csv`, `ssGSEA_Hallmark_scores_tumors_and_PDO_joint_run.csv`, `joint_tumor_PDO_VST_matrix.rds`
- Check: `ssGSEA_tumor_scores_saved_vs_joint_run.csv` (Spearman rho per pathway between the saved tumor scores and the joint-run tumor scores; shows whether adding the PDOs or refitting the VST changes the tumor scores)
- Comparison: `ssGSEA_Hallmark_PDO_vs_tumors.csv` (PDO scores, tumor median and quartiles, per-cluster medians, PDO mean in tumor SDs, PDO percentile within the tumors, exploratory Wilcoxon p and FDR)
- Figures: `ssGSEA_Hallmark_heatmap_PDO_vs_tumors` (SDs from the tumor mean, clipped at +/-3; tumors by cluster, then PDOs), `ssGSEA_Hallmark_boxplots_PDO_vs_tumors_significant` (significant pathways only, FDR < 0.05; extra `_pageN` files if more than 12), `ssGSEA_Hallmark_PDO_top_pathways_barplot`
- **PDO-only plots** (second ssGSEA run on the 4 PDOs alone, from the PDO-only VST matrix; scores are relative to the PDOs and are **not on the tumor scale**, so do not place them next to joint-run values): `ssGSEA_Hallmark_PDOonly_top_bottom_pathways` (top 5 and bottom 5 by mean score across the 4 PDOs) and `ssGSEA_Hallmark_PDOonly_largest_range_pathways` (top 10 by range across the 4 PDOs). The PDO-only scores are saved in `ssGSEA_Hallmark_scores_PDO_only_run.csv`. Setting `PDO_PLOT_SCORES <- "joint"` in the script makes these two plots use the joint-run scores instead (file names then contain `joint` instead of `PDOonly`).
- Tables for these figures: `ssGSEA_Hallmark_PDO_vs_tumors_significant_pathways.csv`, `ssGSEA_Hallmark_PDOonly_top_bottom_pathways.csv`, `ssGSEA_Hallmark_PDOonly_largest_range_pathways.csv`, `ssGSEA_Hallmark_PDOonly_range_all_pathways.csv`

### NMF projection (`PDO_NMF_projection.R`) → `NMF_projection/`
Places the PDOs in the k = 3 NMF using the saved W matrix (best run). Uses a joint VST of tumors and PDOs (same construction as the ssGSEA script).
- **Projection:** with W fixed, each sample's loadings H are solved by non-negative least squares (NNLS) and compared on the tumor scale. Each PDO is assigned to the nearest tumor cluster in this loading space.
- **Validation:** the 31 tumors are projected the same way (leave-one-out) and the script reports how many return to their own cluster. **If this is below about 80%, the expression input does not match the NMF input, and the PDO placement should not be interpreted.** Metagene numbers are not assumed to equal cluster numbers; the map is learned from the tumors.
- **Top-weighted genes:** top 50 genes per metagene (`TOP_N`) for the heatmap, signature scores and centroid correlation. Each gene belongs to the metagene where its weight is highest.
- **PCA:** top 2000 genes overall by highest weight (`PCA_N`), with axes from the **31 tumors only**; the PDOs are projected into that space and do not influence the axes.
- Tables: `NMF_projection_PDO.csv`, `NMF_projection_tumors.csv`, `NMF_projection_validation_tumors_confusion.csv`, `NMF_metagene_to_cluster_map.csv`, `NMF_top_weighted_genes.csv`, `NMF_top_gene_signature_scores.csv`, `NMF_top_genes_centroid_correlation_PDO.csv`, `NMF_PCA_weighted_genes.csv`, `NMF_PCA_weighted_genes_coordinates.csv`, `NMF_PCA_weighted_genes_variance.csv`
- Figures: `NMF_projection_loadings_heatmap`, `NMF_top_genes_heatmap` (rows ordered Cluster 1, 2, 3, then the PDOs), `NMF_top_gene_signature_scores_boxplot`, `NMF_PCA_weighted_genes_tumors_with_PDO_projected`
- Interpretation: the top genes of two metagenes are mainly microenvironment genes (fat and stroma; immune and matrix), which the PDOs lack. The PDOs therefore score low on those two metagenes by composition. The metagene of luminal tumor-cell genes (ESR1, PGR, GATA3, FOXA1, ...) is the informative one for the PDOs.

### cGAS-STING signature (`PDO_cGAS_STING.R`) → `Pathway_analysis/cGAS_STING/`
Repeats the tumor cGAS-STING score for the PDOs. The 7-gene panel (Chen et al., npj Breast Cancer 2024: CGAS, STING1, TBK1, IRF3, STAT1, CCL5, CXCL10) is scored as the mean of per-gene z-scores, with a sensing core (CGAS, STING1, TBK1, IRF3) and a downstream output (CCL5, CXCL10, STAT1).
- **Scale:** joint VST of tumors and PDOs; z-scores use the mean and SD of the 31 tumors, so a PDO value is "SDs from the tumor mean".
- **Statistics:** only the percentile of each PDO among the 31 tumors and an exploratory Wilcoxon test (4 PDOs vs 31 tumors, BH). The cluster tests, TILs-adjusted models, grade test and clinical correlations of the tumor script cannot be repeated (the PDOs have no cluster, TILs or grade).
- **Checks printed to the console:** the cluster labels against `Sample_metadata.tsv` of the tumor script, and the tumor score from the joint VST against the tumor script's own VST.
- Tables: `PDO_cGAS_STING_gene_map.tsv`, `..._scores_and_gene_expression_per_sample.tsv`, `..._PDO_vs_tumors_per_PDO.tsv`, `..._PDO_vs_tumors_Wilcoxon_exploratory.tsv`
- Figures: `..._boxplots_tumors_by_cluster_with_PDOs`, `..._scores_tumors_by_cluster_with_PDOs`, `..._heatmap_tumors_by_cluster_then_PDOs`
- Interpretation: PDOs have no immune cells, so the output genes (interferon and T-cell genes) are expected to be low. The sensing core is the fairer comparison.

### cGAS-STING in the cancer-cell part (`PDO_cGAS_STING_cancer_cell.R`) → `InstaPrism_PBCA_reference/Downstream/cGAS_STING_cancer_cell/`
Same 7-gene panel and score as `PDO_cGAS_STING.R`, but also using the InstaPrism (PBCA) cancer-epithelial expression estimate (Z) of each tumor, read from the saved deconvolution results (InstaPrism is not rerun). The PDOs are almost pure cancer cells, so this removes the part of the PDO-versus-tumor difference that comes from the tumors' immune and stromal cells.
- **Two views side by side:** whole-sample log2(TPM + 1) and cancer-cell log2(Z + 1), both z-scored with the tumor mean and SD; PDO percentile = % of the 31 tumors at or below the PDO.
- **Caution:** Z is a model estimate, noisier for low-expressed genes and tumors with a low cancer fraction (the fraction is saved and shown as a strip).
- Tables: `..._gene_map.tsv`, `..._values_per_sample.tsv`, `..._PDO_vs_tumors_per_PDO.tsv` (both views), `..._PDO_vs_tumors_Wilcoxon_exploratory.tsv`, `..._cancer_fraction_per_sample.tsv`, `..._rank_in_cancer_epithelial_gene_ranking.tsv` (only if the ranking file is found)
- Figures: `..._PDO_percentile_bulk_vs_cancer_cells`, `..._boxplots_tumors_by_cluster_with_PDOs`, `..._heatmap_tumors_by_cluster_then_PDOs`

## How the PDOs are put on the tumor scale

| Method | How PDO and tumor values are made comparable |
|--------|-----------------------------------------------|
| PAM50, ssGSEA, NMF projection, cGAS-STING | **Joint** VST from the raw counts of tumors and PDOs (blind, intercept only), then one run on all 35 samples |
| decoupleR, cancer-epithelial genes, ssGSEA and NMF heatmaps | Genes or pathways are z-scored with the **tumor** mean and SD, so values are "SDs from the tumor average" (PDOs do not enter the mean or SD) |
| PCA of the PDOs, TPM, ESTIMATE, InstaPrism | PDOs only; compared with the saved tumor results only in comparison figures |

Heatmaps clip at +/-3 SDs. Many PDO cells reach the limit because the PDOs sit outside the tumor range, and clipped colors cannot show how far outside. Use the CSV values (unclipped z-scores and percentiles) for statements.

## Notes
- **BO-268:** the BAM name contains an extra `B` (`BO268_S36_B...`). It has the lowest depth but the most genes with counts. Check the metadata to see whether this is a batch, a second library or a re-sequencing run.
- Gene IDs in `PDO_raw_counts.csv` keep their Ensembl version numbers (e.g. `ENSG00000290825.1`). The scripts remove them before mapping to gene symbols.
- `PDO_raw_counts.csv` contains all annotated genes. Low-expression filtering is done in the VST script and again, with the same rule, inside the joint-VST scripts (PAM50, ssGSEA, NMF projection).
- ssGSEA scores come from the joint run with the tumors, so they can be compared with the tumor scores. The tumor scores from this joint run can differ slightly from the saved tumor-only scores; the check file reports how much.
- Ranking pathways by mean or range of the absolute ssGSEA scores across the 4 PDOs (top/bottom and range figures) describes the PDOs only. It says nothing about how they compare with the tumors; use `ssGSEA_Hallmark_PDO_vs_tumors.csv` for that. Absolute scores also depend on the gene set itself (size, how highly its genes are expressed), so low or high absolute scores can reflect the gene set and not PDO biology.
