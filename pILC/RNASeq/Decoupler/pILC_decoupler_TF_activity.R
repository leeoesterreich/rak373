# =====================================================================
# pILC: CollecTRI network preparation + decoupleR TF activity by RNA cluster
# k = 3 clusters (50-run NMF, GenePattern numbering)
# =====================================================================
# Install once if needed:
# BiocManager::install(c("decoupleR", "ComplexHeatmap"))

suppressPackageStartupMessages({
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

SEED <- 42
set.seed(SEED)

# ---- 0. Settings ----------------------------------------------------
root    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
out_dir <- file.path(root, "Decoupler")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

vst_file <- file.path(
  root, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47",
  "htseq_VST_matrix_all_expression_filtered_genes.rds"
)
cluster_file <- file.path(
  root, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv"
)
net_gz   <- file.path(out_dir, "collectri-human.tsv.gz")  # downloaded file
net_file <- file.path(out_dir, "collectri_human.tsv")     # converted network

n_top   <- 10                          # top TFs highlighted per cluster
key_tfs <- c("AR", "ESR1", "FOXA1")    # always shown in the boxplot

# ---- 1. TF-target network (CollecTRI) -------------------------------
# Converted once from the downloaded table; later runs reuse the .tsv
if (!file.exists(net_file)) {
  if (!file.exists(net_gz)) stop("Downloaded network not found: ", net_gz)

  raw <- read_tsv(net_gz, show_col_types = FALSE,
                  col_types = cols(.default = "c"))
  stopifnot(all(c("source_genesymbol", "target_genesymbol",
                  "is_stimulation", "is_inhibition") %in% colnames(raw)))
  cat("CollecTRI rows read:", nrow(raw), "\n")

  is_true <- function(x) x %in% c("True", "TRUE", "true")

  net0 <- raw %>%
    transmute(
      source = source_genesymbol,
      target = target_genesymbol,
      stim   = is_true(is_stimulation),
      inhib  = is_true(is_inhibition)
    ) %>%
    mutate(mor = case_when(stim & !inhib ~  1,
                           inhib & !stim ~ -1,
                           TRUE          ~ NA_real_))

  # Report what is removed instead of dropping it silently
  cat("Unsigned/ambiguous (both or neither flag):", sum(is.na(net0$mor)), "\n")
  cat("Missing symbols:", sum(is.na(net0$source) | is.na(net0$target)), "\n")
  net0 <- net0 %>% filter(!is.na(mor), !is.na(source), !is.na(target))

  # Same TF-target pair with conflicting signs is ambiguous: remove it
  conflict <- net0 %>% group_by(source, target) %>%
    filter(n_distinct(mor) > 1) %>% ungroup() %>% distinct(source, target)
  cat("Pairs with conflicting signs removed:", nrow(conflict), "\n")

  net0 <- net0 %>%
    anti_join(conflict, by = c("source", "target")) %>%
    distinct(source, target, mor)
  write_tsv(net0, net_file)
  cat("Saved network:", net_file, "\n")
}

net <- read_tsv(net_file, show_col_types = FALSE)
stopifnot(all(c("source", "target", "mor") %in% colnames(net)))
cat("Network:", nrow(net), "interactions |",
    length(unique(net$source)), "TFs (incl. complexes) |",
    length(unique(net$target)), "targets\n")
print(table(net$mor))
cat("AR, ESR1, FOXA1 present as TFs:",
    paste(key_tfs %in% net$source, collapse = ", "), "\n")

# ---- 2. Load expression and cluster labels --------------------------
expr <- as.matrix(readRDS(vst_file))
stopifnot(is.numeric(expr), !is.null(rownames(expr)), !is.null(colnames(expr)))
stopifnot(!anyDuplicated(rownames(expr)))
cat("VST matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Example row names:", paste(head(rownames(expr), 5), collapse = ", "), "\n")

# Row names must be gene symbols to match the network
n_overlap <- length(intersect(net$target, rownames(expr)))
cat("Network targets found in expression matrix:", n_overlap, "\n")
if (n_overlap < 1000) {
  stop("Very low overlap: row names are probably not gene symbols. ",
       "Map them with gene_annotation_from_GTF.tsv before running.")
}

clu <- read_csv(cluster_file, show_col_types = FALSE)
stopifnot(all(c("sample", "cluster") %in% colnames(clu)))
cl_raw <- as.character(clu$cluster)
clu$cluster_lab <- ifelse(grepl("^[0-9]+$", cl_raw),
                          paste0("Cluster", cl_raw), cl_raw)

# Harmonize sample IDs (hyphen vs dot); stop if anything does not match
norm_id <- function(x) gsub("-", ".", x, fixed = TRUE)
colnames(expr) <- norm_id(colnames(expr))
clu$sample     <- norm_id(clu$sample)

missing_in_expr <- setdiff(clu$sample, colnames(expr))
if (length(missing_in_expr) > 0) {
  stop("Cluster samples not found in VST columns: ",
       paste(missing_in_expr, collapse = ", "))
}
expr <- expr[, clu$sample, drop = FALSE]
cl   <- factor(clu$cluster_lab, levels = sort(unique(clu$cluster_lab)))
names(cl) <- clu$sample
print(table(cl))

pal     <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
cl_cols <- setNames(pal[seq_along(levels(cl))], levels(cl))

# ---- 3. Z-score each gene across samples ----------------------------
# Removes baseline expression so scores reflect differences between samples
keep   <- apply(expr, 1, sd) > 0
expr_z <- t(scale(t(expr[keep, ])))

# ---- 4. TF activity per sample (univariate linear model) ------------
tf_long <- run_ulm(
  mat = expr_z, network = net,
  .source = source, .target = target, .mor = mor,
  minsize = 5
) %>% filter(statistic == "ulm")

tf_mat <- tf_long %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
tf_mat <- tf_mat[, names(cl), drop = FALSE]
cat("TFs scored:", nrow(tf_mat), "\n")

# ---- 5. Statistics across clusters ----------------------------------
# 5a. Global test across all clusters
kw_p  <- apply(tf_mat, 1, function(x) kruskal.test(x ~ cl)$p.value)
kw_df <- tibble(TF = rownames(tf_mat), KW_p = kw_p,
                KW_FDR = p.adjust(kw_p, "BH"))

# 5b. One-vs-rest for every cluster
ovr <- bind_rows(lapply(levels(cl), function(k) {
  in_k <- cl == k
  tibble(
    TF        = rownames(tf_mat),
    Cluster   = k,
    mean_in   = rowMeans(tf_mat[, in_k, drop = FALSE]),
    mean_rest = rowMeans(tf_mat[, !in_k, drop = FALSE]),
    p = apply(tf_mat, 1, function(x)
      suppressWarnings(wilcox.test(x[in_k], x[!in_k], exact = FALSE)$p.value))
  ) %>%
    mutate(diff = mean_in - mean_rest, FDR = p.adjust(p, "BH"))
})) %>% left_join(kw_df, by = "TF")

# 5c. Top TFs per cluster (higher activity in that cluster)
top <- ovr %>%
  filter(diff > 0) %>%
  group_by(Cluster) %>%
  arrange(p, desc(diff), .by_group = TRUE) %>%
  slice_head(n = n_top) %>%
  ungroup() %>%
  mutate(Significant = FDR < 0.05)

write_csv(as_tibble(tf_mat, rownames = "TF"),
          file.path(out_dir, "TF_activity_per_sample.csv"))
write_csv(ovr, file.path(out_dir, "TF_cluster_stats.csv"))
write_csv(top, file.path(out_dir, "TF_top_per_cluster.csv"))

# 5d. Lists of all significant TFs (not only the top 10 per cluster)
# - Kruskal-Wallis: difference somewhere among the clusters (BH over TFs)
# - One-vs-rest: each cluster against the others (BH over TFs), with direction
kw_sig <- kw_df %>% filter(KW_FDR < 0.05) %>% arrange(KW_FDR)
sig_ovr <- ovr %>%
  filter(FDR < 0.05) %>%
  mutate(direction = ifelse(diff > 0, "Higher in cluster", "Lower in cluster")) %>%
  arrange(Cluster, direction, desc(abs(diff)))
write_csv(kw_sig, file.path(out_dir, "TF_significant_Kruskal_Wallis.csv"))
write_csv(sig_ovr, file.path(out_dir, "TF_significant_one_vs_rest_all_clusters.csv"))
cat("TFs tested:", nrow(tf_mat), "| Kruskal-Wallis FDR < 0.05:", nrow(kw_sig), "\n")
print(count(sig_ovr, Cluster, direction))

# ---- 6. Plots (PNG + PDF from one helper) ---------------------------
save_both <- function(name, w, h, draw_fun) {
  pdf(file.path(out_dir, paste0(name, ".pdf")), width = w, height = h)
  draw_fun(); dev.off()
  png(file.path(out_dir, paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  draw_fun(); dev.off()
}

# 6a. Heatmap of the union of top TFs, samples grouped by cluster
ord    <- order(cl)
hm_tfs <- unique(top$TF)
hm_mat <- t(scale(t(tf_mat[hm_tfs, ord, drop = FALSE])))  # row z-score for display
ha <- HeatmapAnnotation(
  Cluster = cl[ord], col = list(Cluster = cl_cols),
  annotation_name_side = "left"
)
ht <- Heatmap(
  hm_mat, name = "Row z-score\nTF activity",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  top_annotation = ha, column_split = cl[ord],
  cluster_columns = FALSE, cluster_rows = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7), row_names_gp = gpar(fontsize = 8)
)
save_both("TF_activity_heatmap_top_TFs", 9, 7, function() draw(ht))

# 6b. Bar plot: top TFs per cluster (mean activity difference vs rest)
top_plot <- top %>%
  mutate(TF_f = factor(paste(Cluster, TF, sep = "|"),
                       levels = rev(paste(Cluster, TF, sep = "|"))))
p_bar <- ggplot(top_plot, aes(x = diff, y = TF_f, fill = Cluster)) +
  geom_col() +
  geom_text(aes(label = ifelse(Significant, "*", "")), hjust = -0.3) +
  facet_wrap(~Cluster, scales = "free_y", ncol = 3) +
  scale_y_discrete(labels = function(x) sub("^.*\\|", "", x)) +
  scale_fill_manual(values = cl_cols, guide = "none") +
  labs(x = "Mean TF activity difference (cluster vs rest)", y = NULL,
       caption = "* FDR < 0.05 (one-vs-rest Wilcoxon, BH)") +
  theme_bw(base_size = 11)
save_both("TF_activity_top_TFs_barplot", 11, 5, function() print(p_bar))

# 6c. Boxplot: key TFs plus the top TF of each cluster
box_tfs <- unique(c(key_tfs, top %>% group_by(Cluster) %>%
                      slice_head(n = 1) %>% pull(TF)))
box_tfs <- intersect(box_tfs, rownames(tf_mat))

# Kruskal-Wallis FDR (across clusters, BH over all scored TFs) shown in each panel header
fdr_txt <- function(x) ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))
lab_df <- kw_df %>%
  filter(TF %in% box_tfs) %>%
  mutate(panel = paste0(TF, "\nKruskal-Wallis FDR = ", fdr_txt(KW_FDR)))
panel_levels <- lab_df$panel[match(box_tfs, lab_df$TF)]

box_df <- as_tibble(tf_mat[box_tfs, , drop = FALSE], rownames = "TF") %>%
  pivot_longer(-TF, names_to = "sample", values_to = "score") %>%
  mutate(Cluster = cl[sample]) %>%
  left_join(select(lab_df, TF, panel), by = "TF") %>%
  mutate(panel = factor(panel, levels = panel_levels))
p_box <- ggplot(box_df, aes(x = Cluster, y = score, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2) +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = cl_cols, guide = "none") +
  labs(y = "TF activity (ulm score)", x = NULL,
       caption = "FDR: Kruskal-Wallis across clusters, BH-adjusted over all scored TFs") +
  theme_bw(base_size = 11)
save_both("TF_activity_key_TFs_boxplot", 10, 7, function() print(p_box))

cat("Done. Results in:", out_dir, "\n")