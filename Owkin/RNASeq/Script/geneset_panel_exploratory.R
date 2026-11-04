# ==============================================================================
# Script 5: exploratory panel of all 11 lab gene sets + ESR1 percentile
#           8 Primary vs 8 Recurrence (paired), same style as the ESR1 activity figure
#   A. Paired tests per feature (exact Wilcoxon, n_up) - exploratory, no selection
#   B. Panel figure: boxplot + paired lines + points (colour = patient, shape = recurrence type)
#   C. Aids to choose the sets for an ER-output index, WITHOUT using the recurrence result:
#        - Spearman correlation of every feature with ESR1 percentile (16 samples)
#        - correlation matrix between features
#        - gene overlap between sets
# Input : script 2 output (geneset_composite_scores_per_sample.csv) and the lab workbook
# Output: .../Okwin/RNASeq/ER_signaling
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(readxl)
  library(ggplot2)
})

SEED <- 42
set.seed(SEED)

# ---- 0. Settings (CHANGE THESE) --------------------------------------------------
BASE         <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/RNASeq"
OUT_DIR      <- file.path(BASE, "ER_signaling")
SCORES_FILE  <- file.path(OUT_DIR, "geneset_composite_scores_per_sample.csv")   # from script 2
GENESET_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/ER_scoring_gene_sets.xlsx"

# The 11 lab gene sets, in the order of the lab workbook
GENE_SETS <- c("er_upreg", "er_pos_corr", "er_neg_corr",
               "estrogene_v2_esr1_mutant_up", "estrogene_v2_esr1_mutant_down",
               "Y537S_unique_genes", "D538G_unique_genes",
               "low_dose_tamoxifen_repressed", "high_dose_tamoxifen_repressed",
               "TamR_up", "TamR_down")
FEATURES <- c("ESR1_percentile", GENE_SETS)

# Recurrence type of each patient (3 local, 5 distant metastases)
RECURRENCE_TYPE <- c(PI_R_048 = "Metastasis", PI_R_090 = "Metastasis", PI_R_104 = "Metastasis",
                     PI_R_105 = "Metastasis", PI_R_111 = "Metastasis",
                     PI_R_106 = "Local recurrence", PI_R_113 = "Local recurrence",
                     PI_R_117 = "Local recurrence")

# ---- Helpers ---------------------------------------------------------------------
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# Paired test on one feature: rows = patient x timepoint, columns = patient, timepoint, value
paired_test <- function(df) {
  w <- df %>%
    select(patient, timepoint, value) %>%
    pivot_wider(names_from = timepoint, values_from = value) %>%
    drop_na()
  tibble(
    n_pairs      = nrow(w),
    median_delta = median(w$Recurrence - w$Primary),
    mean_delta   = mean(w$Recurrence - w$Primary),
    n_up         = sum(w$Recurrence > w$Primary),
    p_wilcoxon   = suppressWarnings(wilcox.test(w$Recurrence, w$Primary, paired = TRUE, exact = TRUE)$p.value),
    p_ttest      = t.test(w$Recurrence, w$Primary, paired = TRUE)$p.value
  )
}

# ---- 1. Read the scores from script 2 --------------------------------------------
if (!file.exists(SCORES_FILE)) stop("Run script 2 first; not found: ", SCORES_FILE)
scores <- read.csv(SCORES_FILE, stringsAsFactors = FALSE)
stopifnot(all(FEATURES %in% colnames(scores)), all(scores$patient %in% names(RECURRENCE_TYPE)))
scores$timepoint       <- factor(scores$timepoint, levels = c("Primary", "Recurrence"))
scores$recurrence_type <- RECURRENCE_TYPE[scores$patient]

long <- scores %>%
  select(sample, patient, timepoint, recurrence_type, all_of(FEATURES)) %>%
  pivot_longer(all_of(FEATURES), names_to = "feature", values_to = "value")

# ---- 2. Paired tests (exploratory; no feature is selected from this table) --------
stats <- long %>%
  group_by(feature) %>% group_modify(~ paired_test(.x)) %>% ungroup() %>%
  mutate(FDR_wilcoxon_BH_12 = p.adjust(p_wilcoxon, "BH"),    # over the 12 features shown
         feature = factor(feature, levels = FEATURES)) %>%
  arrange(feature)
write.csv(stats, file.path(OUT_DIR, "geneset_panel_exploratory_paired_stats.csv"), row.names = FALSE)
print(as.data.frame(stats))

# ---- 3. Panel figure (same style as the ESR1 activity figure) --------------------
panel_labels <- stats %>%
  mutate(panel = sprintf("%s\nWilcoxon p = %.3f | %d of %d up", feature, p_wilcoxon, n_up, n_pairs)) %>%
  select(feature, panel)
long_plot <- long %>%
  left_join(panel_labels, by = "feature") %>%
  mutate(panel = factor(panel, levels = panel_labels$panel))

p_panel <- ggplot(long_plot, aes(timepoint, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = patient), colour = "grey60") +
  geom_point(aes(colour = patient, shape = recurrence_type), size = 2.6) +
  facet_wrap(~panel, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "Score (ESR1 percentile or GSVA)", colour = "Patient", shape = "Recurrence type") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))
save_fig(p_panel, "geneset_panel_exploratory", w = 16, h = 11)

# ---- 3b. Single panel: ESR1 expression (same style as the decoupleR ESR1 activity figure) ----
esr1_row <- stats %>% filter(feature == "ESR1_percentile")
p_esr1_expr <- long %>%
  filter(feature == "ESR1_percentile") %>%
  ggplot(aes(timepoint, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = patient), colour = "grey60") +
  geom_point(aes(colour = patient), size = 3) +
  labs(x = NULL, y = "ESR1 expression (percentile among the 16 samples)", colour = "Patient",
       title = sprintf("ESR1 expression, Wilcoxon p = %.3f", esr1_row$p_wilcoxon)) +
  theme_bw(base_size = 12)
save_fig(p_esr1_expr, "ESR1_expression_paired", w = 5, h = 5)

# ---- 4. Aids for choosing sets, independent of the recurrence result -------------
# 4a. Correlation of the features across the 16 samples (descriptive: each patient
#     contributes two samples, so the p-values are not used)
wide <- scores %>% select(all_of(FEATURES))
cor_mat <- cor(wide, method = "spearman")
write.csv(as.data.frame(cor_mat) %>% rownames_to_column("feature"),
          file.path(OUT_DIR, "geneset_panel_feature_correlation.csv"), row.names = FALSE)

cor_long <- as.data.frame(cor_mat) %>%
  rownames_to_column("feature_a") %>%
  pivot_longer(-feature_a, names_to = "feature_b", values_to = "rho") %>%
  mutate(feature_a = factor(feature_a, levels = FEATURES),
         feature_b = factor(feature_b, levels = rev(FEATURES)))
p_cor <- ggplot(cor_long, aes(feature_a, feature_b, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_cor, "geneset_panel_feature_correlation", w = 9, h = 7.5)

# 4b. Gene overlap between the sets (genes as listed in the workbook)
gs_long <- read_excel(GENESET_FILE, sheet = "Genes (long)")
colnames(gs_long) <- c("set", "category", "gene")
gene_sets <- lapply(split(gs_long$gene, gs_long$set), function(g) unique(na.omit(g)))
stopifnot(all(GENE_SETS %in% names(gene_sets)))
overlap <- sapply(GENE_SETS, function(a) sapply(GENE_SETS, function(b)
  length(intersect(gene_sets[[a]], gene_sets[[b]]))))
write.csv(as.data.frame(overlap) %>% rownames_to_column("gene_set"),
          file.path(OUT_DIR, "geneset_panel_gene_overlap.csv"), row.names = FALSE)

# 4c. One summary table: biology-based selection criteria next to the paired result
set_summary <- tibble(
  feature = FEATURES,
  n_genes = c(NA, lengths(gene_sets[GENE_SETS])),
  rho_with_ESR1_percentile = cor_mat[FEATURES, "ESR1_percentile"],
  max_overlap_with_er_upreg_or_er_pos_corr =
    c(NA, sapply(GENE_SETS, function(s) max(overlap[s, c("er_upreg", "er_pos_corr")])))
) %>%
  left_join(stats %>% mutate(feature = as.character(feature)) %>%
              select(feature, median_delta, n_up, p_wilcoxon), by = "feature")
write.csv(set_summary, file.path(OUT_DIR, "geneset_panel_selection_table.csv"), row.names = FALSE)
print(as.data.frame(set_summary))

sessionInfo()