# =====================================================================
# PDO: decoupleR TF activity (CollecTRI, ulm), shown against the pILC tumors
# 4 PDOs, no RNA clusters. Tumor clusters = k3 (50-run NMF, GenePattern numbering)
#
# Run first: pILC_decoupler_TF_activity.R (creates the CollecTRI network file)
#            PDO_counts_to_VST_GENCODEv47.R   (creates the PDO VST matrix)
#
# Design
#  * Each gene is z-scored with the mean and SD of the 31 TUMORS (not across the 4 PDOs,
#    which would be unstable). PDO values are therefore "SDs from the tumor average".
#  * Tumor TF scores are recomputed here on the same gene set, so PDOs and tumors are
#    directly comparable (they can differ slightly from the main tumor run).
#  * Statistics: only PDOs versus all tumors (Wilcoxon, BH). EXPLORATORY: the two groups
#    also differ in culture conditions, VST fit and absence of stromal/immune cells.
# =====================================================================

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
pdo_dir <- file.path(root, "PDO")
out_dir <- file.path(pdo_dir, "Decoupler")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

net_file     <- file.path(root, "Decoupler", "collectri_human.tsv")   # made by the tumor script
vst_t_file   <- file.path(root, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47",
                          "htseq_VST_matrix_all_expression_filtered_genes.rds")
vst_p_file   <- file.path(pdo_dir, "VST_data_GENCODEv47",
                          "htseq_VST_matrix_all_expression_filtered_genes.rds")
cluster_file <- file.path(root, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
for (f in c(net_file, vst_t_file, vst_p_file, cluster_file))
  if (!file.exists(f)) stop("Input not found: ", f)

n_top   <- 15                          # TFs shown per direction (higher / lower in PDOs)
key_tfs <- c("AR", "ESR1", "FOXA1")    # always shown in the heatmap

# ---- 1. Network (same file as the tumor analysis) -------------------
net <- read_tsv(net_file, show_col_types = FALSE)
stopifnot(all(c("source", "target", "mor") %in% colnames(net)))
cat("Network:", nrow(net), "interactions |",
    length(unique(net$source)), "TFs (incl. complexes) |",
    length(unique(net$target)), "targets\n")

# ---- 2. Expression and tumor cluster labels -------------------------
expr_t <- as.matrix(readRDS(vst_t_file))
expr_p <- as.matrix(readRDS(vst_p_file))
stopifnot(is.numeric(expr_t), is.numeric(expr_p),
          !anyDuplicated(rownames(expr_t)), !anyDuplicated(rownames(expr_p)))
if (ncol(expr_p) != 4) stop("Expected 4 PDO samples, found ", ncol(expr_p))
cat("Tumor VST:", nrow(expr_t), "genes x", ncol(expr_t), "samples\n")
cat("PDO VST  :", nrow(expr_p), "genes x", ncol(expr_p), "samples (", paste(colnames(expr_p), collapse = ", "), ")\n")

n_overlap <- length(intersect(net$target, rownames(expr_p)))
cat("Network targets found in the PDO matrix:", n_overlap, "\n")
if (n_overlap < 1000) stop("Very low overlap: row names are probably not gene symbols.")

# Tumor samples: harmonize IDs (hyphen vs dot) exactly as in the tumor script
norm_id <- function(x) gsub("-", ".", x, fixed = TRUE)
clu <- read_csv(cluster_file, show_col_types = FALSE)
stopifnot(all(c("sample", "cluster") %in% colnames(clu)))
cl_raw <- as.character(clu$cluster)
clu$cluster_lab <- ifelse(grepl("^[0-9]+$", cl_raw), paste0("Cluster", cl_raw), cl_raw)
colnames(expr_t) <- norm_id(colnames(expr_t))
clu$sample       <- norm_id(clu$sample)
missing_in_expr <- setdiff(clu$sample, colnames(expr_t))
if (length(missing_in_expr) > 0)
  stop("Cluster samples not found in tumor VST columns: ", paste(missing_in_expr, collapse = ", "))
expr_t <- expr_t[, clu$sample, drop = FALSE]
cl <- factor(clu$cluster_lab, levels = sort(unique(clu$cluster_lab)))
names(cl) <- clu$sample
print(table(cl))

pal     <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
cl_cols <- setNames(pal[seq_along(levels(cl))], levels(cl))
pdo_cols <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[seq_len(ncol(expr_p))], colnames(expr_p))

# ---- 3. z-score every gene on the TUMOR scale -----------------------
common <- intersect(rownames(expr_t), rownames(expr_p))
cat("Genes shared by tumor and PDO matrices:", length(common), "\n")
et  <- expr_t[common, , drop = FALSE]
ep  <- expr_p[common, , drop = FALSE]
mu  <- rowMeans(et)
sdv <- apply(et, 1, sd)
ok  <- is.finite(sdv) & sdv > 0
Z_t <- (et[ok, , drop = FALSE] - mu[ok]) / sdv[ok]
Z_p <- (ep[ok, , drop = FALSE] - mu[ok]) / sdv[ok]
cat("Genes used for scoring:", sum(ok), "\n")

# Sanity check: PDO z-values should be centred near 0 across genes if the two VST fits are on a similar scale
cat("Median PDO z-score per sample (expect close to 0):\n")
print(round(apply(Z_p, 2, median), 3))
cat("Share of PDO z-values with |z| > 5:\n")
print(round(apply(Z_p, 2, function(x) mean(abs(x) > 5)), 4))

# ---- 4. TF activity (univariate linear model) -----------------------
run_tf <- function(Z) {
  run_ulm(mat = Z, network = net, .source = source, .target = target, .mor = mor, minsize = 5) %>%
    filter(statistic == "ulm")
}
to_mat <- function(res, cols) {
  m <- res %>%
    select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()
  m[, cols, drop = FALSE]
}
res_t <- run_tf(Z_t)
res_p <- run_tf(Z_p)
tf_t  <- to_mat(res_t, names(cl))
tf_p  <- to_mat(res_p, colnames(expr_p))
tfs   <- intersect(rownames(tf_t), rownames(tf_p))
tf_t  <- tf_t[tfs, , drop = FALSE]
tf_p  <- tf_p[tfs, , drop = FALSE]
cat("TFs scored:", length(tfs), "\n")

# Per-sample ulm p-values for the PDOs (BH within each PDO)
p_tab <- res_p %>%
  select(TF = source, Sample = condition, score, p_value) %>%
  group_by(Sample) %>%
  mutate(FDR_within_sample = p.adjust(p_value, "BH")) %>%
  ungroup()
write_csv(p_tab, file.path(out_dir, "TF_activity_PDO_ulm_pvalues.csv"))

write_csv(as_tibble(tf_p, rownames = "TF"), file.path(out_dir, "TF_activity_per_sample_PDO.csv"))
write_csv(as_tibble(tf_t, rownames = "TF"), file.path(out_dir, "TF_activity_per_sample_tumors_common_genes.csv"))

# ---- 5. PDOs versus tumors ------------------------------------------
tum_mean <- rowMeans(tf_t)
tum_sd   <- apply(tf_t, 1, sd)
cl_means <- sapply(levels(cl), function(k) rowMeans(tf_t[, cl == k, drop = FALSE]))
colnames(cl_means) <- paste0("tumor_mean_", levels(cl))

cmp <- bind_cols(
  tibble(TF = tfs),
  as_tibble(tf_p),
  tibble(PDO_mean = rowMeans(tf_p), tumor_mean_all = tum_mean, tumor_sd_all = tum_sd),
  as_tibble(cl_means)
) %>%
  mutate(PDO_mean_z_vs_tumors   = (PDO_mean - tumor_mean_all) / tumor_sd_all,
         PDO_percentile_in_tumors = sapply(seq_along(tfs), function(i) 100 * mean(tf_t[i, ] <= PDO_mean[i])))

# Exploratory test: 4 PDOs vs 31 tumors (Wilcoxon rank-sum, BH over TFs)
p_w <- sapply(seq_along(tfs), function(i)
  suppressWarnings(wilcox.test(tf_p[i, ], tf_t[i, ], exact = FALSE)$p.value))
stats <- tibble(TF = tfs, mean_PDO = rowMeans(tf_p), mean_tumor = tum_mean,
                diff = rowMeans(tf_p) - tum_mean, p = p_w) %>%
  mutate(FDR = p.adjust(p, "BH")) %>%
  arrange(p)
cmp <- left_join(cmp, select(stats, TF, p_PDO_vs_tumors = p, FDR_PDO_vs_tumors = FDR), by = "TF") %>%
  arrange(desc(abs(PDO_mean_z_vs_tumors)))

write_csv(cmp,   file.path(out_dir, "TF_PDO_vs_tumors_comparison.csv"))
write_csv(stats, file.path(out_dir, "TF_PDO_vs_tumors_stats.csv"))      # input of the volcano script
cat("TFs with FDR < 0.05 (PDO vs tumors, exploratory):", sum(stats$FDR < 0.05), "\n")

# ---- 6. Plots (PNG + PDF from one helper) ---------------------------
save_both <- function(name, w, h, draw_fun) {
  pdf(file.path(out_dir, paste0(name, ".pdf")), width = w, height = h); draw_fun(); invisible(dev.off())
  png(file.path(out_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300); draw_fun(); invisible(dev.off())
}

# 6a. Heatmap: tumors by cluster, then PDOs; values = TF activity as SDs from the tumor mean
up   <- cmp %>% arrange(desc(PDO_mean_z_vs_tumors)) %>% slice_head(n = n_top) %>% pull(TF)
down <- cmp %>% arrange(PDO_mean_z_vs_tumors)       %>% slice_head(n = n_top) %>% pull(TF)
hm_tfs <- unique(c(up, down, intersect(key_tfs, tfs)))

ord <- order(cl)
hm_raw <- cbind(tf_t[hm_tfs, ord, drop = FALSE], tf_p[hm_tfs, , drop = FALSE])
hm_mat <- (hm_raw - tum_mean[hm_tfs]) / tum_sd[hm_tfs]
grp    <- factor(c(as.character(cl[ord]), rep("PDO", ncol(tf_p))), levels = c(levels(cl), "PDO"))

ha <- HeatmapAnnotation(
  Group = grp, col = list(Group = c(cl_cols, PDO = "#333333")),
  annotation_name_side = "left"
)
ht <- Heatmap(
  hm_mat, name = "TF activity\n(SD from tumor mean)",
  col = colorRamp2(c(-3, 0, 3), c("navy", "white", "firebrick")),
  top_annotation = ha, column_split = grp,
  cluster_columns = FALSE, cluster_rows = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 7), row_names_gp = gpar(fontsize = 8)
)
save_both("TF_activity_heatmap_PDO_vs_tumors", 11, 8, function() draw(ht))

# 6b. Bar plot: TFs most different in the PDOs (SDs from the tumor mean)
bar_df <- cmp %>%
  filter(TF %in% c(up, down)) %>%
  mutate(Direction = ifelse(PDO_mean_z_vs_tumors > 0, "Higher in PDOs", "Lower in PDOs"),
         TF = factor(TF, levels = TF[order(PDO_mean_z_vs_tumors)]))
p_bar <- ggplot(bar_df, aes(x = PDO_mean_z_vs_tumors, y = TF, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = c("Higher in PDOs" = "#b2182b", "Lower in PDOs" = "#2166ac"), name = NULL) +
  labs(x = "Mean PDO TF activity (SD from the tumor mean)", y = NULL,
       caption = paste0("Top ", n_top, " TFs per direction. Descriptive; PDOs also differ from tumors in culture conditions and cell composition.")) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
save_both("TF_activity_PDO_top_TFs_barplot", 7, 7, function() print(p_bar))

cat("Done. Results in:", out_dir, "\n")
