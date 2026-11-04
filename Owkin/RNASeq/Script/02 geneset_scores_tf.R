# ==============================================================================
# Script 2: lab ER gene sets, composite scores and TF screen
#           8 Primary vs 8 Recurrence (paired)
#   A. GSVA scores of the lab gene sets (ER_scoring_gene_sets.xlsx)
#   B. Composite scores: ER activity index (4 components), Resistance index (lab definition),
#      ER TF-module score (ESR1 / FOXA1 / GATA3 activity)
#   C. TF screen: all CollecTRI TFs, paired limma (patient as blocking factor)
#   D. TFs whose activity tracks the ER activity index / Resistance index
# Output folder: .../Okwin/RNASeq/ER_signaling
# Install once if needed: BiocManager::install(c("GSVA", "limma")); install.packages("readxl")
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)          # also loads limma
  library(decoupleR)
  library(GSVA)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readxl)
  # tidyverse packages are loaded last so their functions are not masked
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})
select <- dplyr::select   # AnnotationDbi masks dplyr::select
rename <- dplyr::rename   # S4Vectors masks dplyr::rename

SEED <- 42
set.seed(SEED)

# ---- 0. Settings (CHANGE THESE) --------------------------------------------------
BASE         <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/RNASeq"
COUNTS_FILE  <- file.path(BASE, "counts/matrix/gene_counts.tsv")
META_FILE    <- file.path(BASE, "sample_metadata.csv")              # block, patient, sample_type
OUT_DIR      <- file.path(BASE, "ER_signaling")   # all results of script 2 and later scripts
GENESET_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/ER_scoring_gene_sets.xlsx"
NET_CANDIDATES <- c("/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Decoupler/collectri_human.tsv",
                    file.path(BASE, "collectri_human.tsv"))         # converted CollecTRI (source, target, mor)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ER_TFS <- c("ESR1", "FOXA1", "GATA3")   # TFs combined into the ER TF-module score

# Pre-specified ER-axis TF panel (ER cofactors, luminal TFs, endocrine-resistance TFs); edit freely
ER_PANEL <- c("ESR1", "FOXA1", "GATA3", "TFAP2C", "GRHL2", "XBP1", "MYB", "SPDEF", "PBX1",
              "NR2F2", "FOXM1", "E2F1", "MYC", "JUN", "FOS", "RELA", "NFKB1", "STAT3", "SOX9", "HIF1A")

# ---- Helpers ---------------------------------------------------------------------
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
zs <- function(x) as.numeric(scale(x))

# Exact paired Wilcoxon + paired t-test on (Recurrence - Primary); df: patient, timepoint, value
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

# If y_labels is given (named vector: feature -> label), each panel gets its own Y-axis label
# (the facet strip is moved to the left and replaces the shared Y title)
# If top_labels is given (named vector: feature -> label), it replaces the strip text above each panel
paired_plot <- function(df, ylab, ncol = NULL, y_labels = NULL, top_labels = NULL) {
  p <- ggplot(df, aes(timepoint, value, group = patient, colour = patient)) +
    geom_line(alpha = 0.8) + geom_point(size = 2) +
    labs(x = NULL, y = ylab, colour = "Patient") +
    theme_bw(base_size = 11)
  if (is.null(y_labels)) {
    if (is.null(top_labels)) p + facet_wrap(~feature, scales = "free_y", ncol = ncol)
    else p + facet_wrap(~feature, scales = "free_y", ncol = ncol, labeller = as_labeller(top_labels))
  } else {
    p + facet_wrap(~feature, scales = "free_y", ncol = ncol, strip.position = "left",
                   labeller = as_labeller(y_labels)) +
      theme(strip.placement = "outside", strip.background = element_blank(),
            strip.text.y.left = element_text(angle = 90, size = 10), axis.title.y = element_blank())
  }
}

# GSVA works with both the old and the new GSVA interface
run_gsva <- function(mat, sets) {
  if (utils::packageVersion("GSVA") >= "1.50.0") {
    GSVA::gsva(GSVA::gsvaParam(mat, sets, kcdf = "Gaussian", minSize = 5), verbose = FALSE)
  } else {
    GSVA::gsva(mat, sets, method = "gsva", kcdf = "Gaussian", min.sz = 5, verbose = FALSE)
  }
}

# ---- 1. Metadata and expression (same preparation as script 1) -------------------
meta <- read.csv(META_FILE, stringsAsFactors = FALSE) %>%
  mutate(
    sample    = paste0(block, "_mRNA"),
    timepoint = factor(ifelse(sample_type == "Primary", "Primary", "Recurrence"),
                       levels = c("Primary", "Recurrence"))
  )
counts <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE)

meta <- meta %>%
  filter(sample %in% colnames(counts)) %>%
  group_by(patient) %>%
  filter(any(timepoint == "Primary"), any(timepoint == "Recurrence")) %>%
  ungroup() %>%
  mutate(patient = factor(patient)) %>%
  arrange(patient, timepoint)
stopifnot(all(table(meta$patient) == 2))
cat("Pairs used:", nlevels(meta$patient), "\n")

counts <- counts[, meta$sample]

# Ensembl ID -> symbol; keep the highest-expressed row per symbol
ens    <- sub("\\..*$", "", rownames(counts))
symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL",
                                column = "SYMBOL", multiVals = "first")
keep_g <- !is.na(symbol)
counts <- counts[keep_g, ]
symbol <- unname(symbol[keep_g])
ord    <- order(rowSums(counts), decreasing = TRUE)
counts <- counts[ord, ]
symbol <- symbol[ord]
dedup  <- !duplicated(symbol)
counts <- counts[dedup, ]
rownames(counts) <- symbol[dedup]

y <- DGEList(counts = counts, group = meta$timepoint)
y <- y[filterByExpr(y, group = meta$timepoint), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 2)
colnames(logcpm) <- meta$sample
stopifnot("ESR1" %in% rownames(logcpm))

saveRDS(logcpm, file.path(OUT_DIR, "logCPM_matrix.rds"))   # read by script 6 (exploratory programs)
keep   <- apply(logcpm, 1, sd) > 0
expr_z <- t(scale(t(logcpm[keep, ])))   # gene-wise z-score across the 16 samples

# ---- 2. Lab gene sets scored with GSVA --------------------------------------------
gs_long <- read_excel(GENESET_FILE, sheet = "Genes (long)")
colnames(gs_long) <- c("set", "category", "gene")        # Gene set, Category, Gene
gene_sets <- lapply(split(gs_long$gene, gs_long$set), function(g) unique(na.omit(g)))

coverage <- tibble(
  gene_set    = names(gene_sets),
  n_listed    = lengths(gene_sets),
  n_in_matrix = sapply(gene_sets, function(g) sum(g %in% rownames(logcpm)))
)
print(as.data.frame(coverage))
write.csv(coverage, file.path(OUT_DIR, "geneset_coverage.csv"), row.names = FALSE)

gs_scores <- run_gsva(logcpm, gene_sets)                 # gene sets x samples
scores <- as.data.frame(t(gs_scores)) %>% rownames_to_column("sample")
stopifnot(identical(scores$sample, meta$sample))

# ---- 3. TF activity for ALL CollecTRI TFs (ulm) -------------------------------------
net_path <- NET_CANDIDATES[file.exists(NET_CANDIDATES)][1]
if (is.na(net_path)) stop("CollecTRI network file not found; check NET_CANDIDATES")
net <- read_tsv(net_path, show_col_types = FALSE)
stopifnot(all(c("source", "target", "mor") %in% colnames(net)))

tf_long <- run_ulm(mat = expr_z, network = net, .source = source, .target = target,
                   .mor = mor, minsize = 5) %>%
  filter(statistic == "ulm")
tf_mat <- tf_long %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
tf_mat <- tf_mat[, meta$sample, drop = FALSE]
cat("TFs scored:", nrow(tf_mat), "\n")
write.csv(as_tibble(tf_mat, rownames = "TF"), file.path(OUT_DIR, "TF_activity_all_per_sample.csv"),
          row.names = FALSE)

# ---- 4. Composite scores -------------------------------------------------------------
# Lab definitions (workbook "About" sheet), z-scores computed across these 16 samples
esr1_pct <- rank(logcpm["ESR1", ]) / ncol(logcpm)        # ESR1 expression percentile
scores$ESR1_percentile   <- esr1_pct[scores$sample]
scores$Resistance_index  <- zs(scores$TamR_up) - zs(scores$TamR_down)
# ER activity index (4 components): the lab index (ESR1 percentile, er_upreg, er_pos_corr)
# plus er_neg_corr with its sign flipped, because er_neg_corr is the set most anti-correlated with
# ESR1 expression (rho = -0.73) and shares almost no genes with the other sets.
# This is an extension of the lab definition (decided after the exploratory panel); the lab
# 3-component version gave the same result (8 of 8 pairs down, exact Wilcoxon p = 0.0078).
scores$ER_activity_index <- (zs(scores$ESR1_percentile) + zs(scores$er_upreg) +
                             zs(scores$er_pos_corr) - zs(scores$er_neg_corr)) / 4

# ER TF-module: mean z-score of ESR1, FOXA1, GATA3 activity
er_tfs_present <- intersect(ER_TFS, rownames(tf_mat))
z_tf <- t(scale(t(tf_mat[er_tfs_present, , drop = FALSE])))
scores$ER_TF_module <- colMeans(z_tf)[scores$sample]

scores <- scores %>%
  left_join(meta %>% select(sample, block, patient, timepoint), by = "sample")
write.csv(scores, file.path(OUT_DIR, "geneset_composite_scores_per_sample.csv"), row.names = FALSE)

# ---- 5. Paired statistics (Recurrence - Primary) ---------------------------------------
composite_feats <- c("ER_activity_index", "Resistance_index", "ER_TF_module")
set_feats       <- names(gene_sets)
all_feats       <- c(composite_feats, set_feats)

long <- scores %>%
  select(patient, timepoint, all_of(all_feats)) %>%
  pivot_longer(all_of(all_feats), names_to = "feature", values_to = "value")

stats <- long %>%
  group_by(feature) %>% group_modify(~ paired_test(.x)) %>% ungroup() %>%
  mutate(family = case_when(feature == "ER_activity_index" ~ "primary",
                            TRUE                           ~ "secondary"))
# Primary hypothesis: ER_activity_index (no adjustment). BH within the secondary family only (13 tests).
sec <- stats$family == "secondary"
stats$FDR_wilcoxon_secondary <- NA_real_
stats$FDR_ttest_secondary    <- NA_real_
stats$FDR_wilcoxon_secondary[sec] <- p.adjust(stats$p_wilcoxon[sec], "BH")
stats$FDR_ttest_secondary[sec]    <- p.adjust(stats$p_ttest[sec], "BH")
stats <- stats %>% arrange(factor(family, levels = c("primary", "secondary")), p_wilcoxon)
write.csv(stats, file.path(OUT_DIR, "geneset_composite_paired_stats.csv"), row.names = FALSE)
print(as.data.frame(stats))

# Result (exact Wilcoxon p and number of pairs going up) written on each panel
result_text <- setNames(sprintf("Wilcoxon p = %.3f | %d of %d up", stats$p_wilcoxon, stats$n_up, stats$n_pairs),
                        stats$feature)
comp_base <- c(ER_activity_index = "ER activity index (z-score units)",
               Resistance_index  = "Resistance index (z-score units)",
               ER_TF_module      = "ER TF-module (z-score units)")
comp_labels <- setNames(paste0(comp_base, "\n", result_text[names(comp_base)]), names(comp_base))
# ER activity index: its own figure (opening slide), same style as the ESR1 activity figure
er_row <- stats %>% filter(feature == "ER_activity_index")
p_er <- long %>%
  filter(feature == "ER_activity_index") %>%
  ggplot(aes(timepoint, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = patient), colour = "grey60") +
  geom_point(aes(colour = patient), size = 3) +
  labs(x = NULL, y = "ER activity index (z-score units)", colour = "Patient",
       title = sprintf("ER activity index\nWilcoxon p = %.4f | %d of %d up", er_row$p_wilcoxon, er_row$n_up, er_row$n_pairs)) +
  theme_bw(base_size = 12)
save_fig(p_er, "ER_activity_index_paired", w = 5.5, h = 5)

# The other two composites together (own Y label and result on each panel)
other_feats <- c("Resistance_index", "ER_TF_module")
save_fig(paired_plot(long %>% filter(feature %in% other_feats) %>%
                       mutate(feature = factor(feature, levels = other_feats)),
                     "", ncol = 2, y_labels = comp_labels[other_feats]),
         "composite_scores_paired", w = 8, h = 4.4)
set_labels <- setNames(paste0(set_feats, "\n", result_text[set_feats]), set_feats)
save_fig(paired_plot(long %>% filter(feature %in% set_feats), "GSVA score", ncol = 4, top_labels = set_labels),
         "geneset_gsva_scores_paired", w = 13, h = 9)

# Change in ER index vs change in other scores (8 pairs, descriptive)
wide <- long %>%
  filter(feature %in% composite_feats) %>%
  pivot_wider(names_from = timepoint, values_from = value) %>%
  mutate(delta = Recurrence - Primary) %>%
  select(patient, feature, delta) %>%
  pivot_wider(names_from = feature, values_from = delta)
cmp <- c("Resistance_index", "ER_TF_module")
delta_cor <- tibble(
  delta_ER_index_vs = cmp,
  spearman_rho = unname(sapply(cmp, function(f) cor(wide$ER_activity_index, wide[[f]], method = "spearman"))),
  p = unname(sapply(cmp, function(f) suppressWarnings(
    cor.test(wide$ER_activity_index, wide[[f]], method = "spearman", exact = FALSE)$p.value)))
)
write.csv(delta_cor, file.path(OUT_DIR, "delta_ER_index_correlations.csv"), row.names = FALSE)
print(delta_cor)

# ---- 6. TF screen: paired limma on TF activity -------------------------------------------
design <- model.matrix(~ patient + timepoint, data = meta)   # patient = blocking factor
fit    <- eBayes(lmFit(tf_mat, design))
tf_res <- topTable(fit, coef = "timepointRecurrence", number = Inf, sort.by = "none") %>%
  rownames_to_column("TF")

delta_tf <- tf_mat[, meta$sample[meta$timepoint == "Recurrence"], drop = FALSE] -
            tf_mat[, meta$sample[meta$timepoint == "Primary"],    drop = FALSE]
n_targets <- net %>% filter(target %in% rownames(expr_z)) %>% count(source, name = "n_targets_used")

tf_res <- tf_res %>%
  mutate(n_up = rowSums(delta_tf > 0)[TF], n_pairs = ncol(delta_tf),
         median_delta = apply(delta_tf, 1, median)[TF],
         in_ER_panel = TF %in% ER_PANEL) %>%
  left_join(n_targets, by = c("TF" = "source")) %>%
  rename(mean_delta = logFC, FDR_all_TFs = adj.P.Val)
# Separate, small multiplicity family for the pre-specified panel
tf_res$FDR_ER_panel <- NA_real_
tf_res$FDR_ER_panel[tf_res$in_ER_panel] <- p.adjust(tf_res$P.Value[tf_res$in_ER_panel], "BH")
tf_res <- tf_res %>% arrange(P.Value)
write.csv(tf_res, file.path(OUT_DIR, "TF_screen_paired_limma.csv"), row.names = FALSE)
cat("Panel TFs not scored (too few targets):", paste(setdiff(ER_PANEL, tf_res$TF), collapse = ", "), "\n")
cat("TFs with FDR < 0.05 (all TFs):", sum(tf_res$FDR_all_TFs < 0.05), "| FDR < 0.10:", sum(tf_res$FDR_all_TFs < 0.10), "\n")
print(head(as.data.frame(tf_res), 15))

# Volcano
tf_plot <- tf_res %>%
  mutate(label = ifelse(in_ER_panel | rank(P.Value) <= 15, TF, ""))
p_volc <- ggplot(tf_plot, aes(mean_delta, -log10(P.Value))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(colour = in_ER_panel), alpha = 0.6) +
  geom_text(aes(label = label), size = 3, vjust = -0.6, check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), name = "ER-axis panel") +
  labs(x = "Mean change in TF activity (Recurrence - Primary)", y = "-log10 P (paired limma)") +
  theme_bw(base_size = 12)
save_fig(p_volc, "TF_screen_volcano", w = 8, h = 6)

# Top 15 up and top 15 down by moderated t
top_tf <- bind_rows(head(arrange(tf_res, -t), 15), head(arrange(tf_res, t), 15)) %>%
  mutate(direction = ifelse(mean_delta > 0, "Higher in recurrence", "Lower in recurrence"),
         TF = factor(TF, levels = TF[order(t)]))
p_top <- ggplot(top_tf, aes(t, TF, fill = direction)) +
  geom_col() +
  labs(x = "Moderated t (paired limma)", y = NULL, fill = NULL) +
  theme_bw(base_size = 11)
save_fig(p_top, "TF_screen_top30", w = 6, h = 8)

# Paired plot of the ER-axis panel
panel_present <- intersect(ER_PANEL, rownames(tf_mat))
panel_long <- as.data.frame(t(tf_mat[panel_present, , drop = FALSE])) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "feature", values_to = "value") %>%
  left_join(meta %>% select(sample, patient, timepoint), by = "sample") %>%
  mutate(feature = factor(feature, levels = panel_present))
tf_labels <- setNames(sprintf("%s\nlimma p = %.3f | %d of %d up", tf_res$TF, tf_res$P.Value, tf_res$n_up, tf_res$n_pairs),
                      tf_res$TF)
save_fig(paired_plot(panel_long, "TF activity (ulm)", ncol = 5, top_labels = tf_labels[panel_present]),
         "TF_ER_panel_paired", w = 15, h = 10)

# ---- 7. TFs whose activity tracks the lab ER index / Resistance index ------------------------
# Spearman across the 16 samples (samples of the same patient are not independent: exploratory)
stopifnot(identical(scores$sample, colnames(tf_mat)))
cor_vs <- function(yv) {
  t(apply(tf_mat, 1, function(x) {
    ct <- suppressWarnings(cor.test(x, yv, method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }))
}
r_er  <- cor_vs(scores$ER_activity_index)
r_res <- cor_vs(scores$Resistance_index)
tf_cor <- tibble(TF = rownames(tf_mat),
                 rho_ER_index = r_er[, "rho"], p_ER_index = r_er[, "p"],
                 rho_Resistance_index = r_res[, "rho"], p_Resistance_index = r_res[, "p"]) %>%
  mutate(FDR_ER_index = p.adjust(p_ER_index, "BH"),
         FDR_Resistance_index = p.adjust(p_Resistance_index, "BH")) %>%
  arrange(p_ER_index)
write.csv(tf_cor, file.path(OUT_DIR, "TF_correlation_with_ER_and_Resistance_index.csv"), row.names = FALSE)
cat("ESR1/FOXA1/GATA3 vs ER index:\n")
print(as.data.frame(tf_cor %>% filter(TF %in% ER_TFS)))

sessionInfo()