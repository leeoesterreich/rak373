###############################################################################
# PDO RNA-seq | Projection of the 4 PDOs onto the pILC k = 3 NMF (best run, matrix W)
#
# Four views:
#  A. NMF projection: with W fixed, solve X_sample ~ W %*% h (h >= 0) by NNLS -> loading of
#     each PDO on each metagene (H). Tumors are projected the same way for validation.
#  B. Top-weighted genes: top N genes per metagene -> signature score per sample, heatmap, boxplots.
#  C. PCA on the union of the top-weighted genes: axes from the 31 tumors only, PDOs projected.
#  D. Correlation of each PDO with the tumor cluster centroids (on the top-weighted genes).
#
# Validation: tumors are projected the same way (leave-one-out nearest centroid) and the script
# reports how many return to their own cluster. If this is low, the expression input differs from
# what the NMF was run on, and the PDO placement should NOT be interpreted.
# Metagene numbers are NOT assumed to equal cluster numbers (clusters use GenePattern labels);
# the metagene -> cluster map is learned from the tumors and printed.
#
# n = 4 PDOs: descriptive only, no p-values. PDOs lack stromal/immune cells, so metagenes driven
# by those cell types will look low in the PDOs for compositional reasons.
# Run first: PDO_counts_to_VST_GENCODEv47.R     Output names carry no dates.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})
if (!requireNamespace("nnls", quietly = TRUE)) stop("Install the nnls package: install.packages('nnls')")

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
ROOT    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR <- file.path(ROOT, "PDO")

W_FILE       <- file.path(ROOT, "NMF/R_NMFConsensus_50runs/NMF_best_run_k3_W.csv")
CNT_T_FILE   <- file.path(ROOT, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv")
CNT_P_FILE   <- file.path(PDO_DIR, "VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv")
CLUSTER_FILE <- file.path(ROOT, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
for (f in c(W_FILE, CNT_T_FILE, CNT_P_FILE, CLUSTER_FILE))
  if (!file.exists(f)) stop("Missing input: ", f)

OUT <- file.path(PDO_DIR, "NMF_projection")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TOP_N         <- 50     # heatmap, signature scores and centroid correlation: top-weighted genes per metagene (= per cluster)
PCA_N         <- 2000   # PCA only: the top-weighted genes overall (ranked by the highest weight over the 3 metagenes)
ASSIGN_UNIQUE <- TRUE   # TRUE: a gene can only be a top gene of the metagene where its weight is highest
MIN_COUNT     <- 5      # low-expression filter of the VST scripts: count >= 5 ...
MIN_FRAC      <- 0.10   # ... in at least 10% of samples
PCS           <- c(1, 2)  # principal components to plot
Z_CLIP        <- 3      # heatmap colour range (SDs from the tumor mean)

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
group_levels   <- c(cluster_levels, "PDO")
group_colors   <- c(cluster_colors, PDO = "#333333")

# ------------------------------- 2. HELPERS ----------------------------------
save_plot <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
}
save_pheatmap <- function(draw_fun, stem, w, h) {
  png(paste0(stem, ".png"), width = w * 300, height = h * 300, res = 300); draw_fun(); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = w, height = h);                        draw_fun(); invisible(dev.off())
}
norm_id <- function(x) {                        # TP09.9_S34 / TP09-9 -> TP09-09 (tumors only)
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}
read_counts <- function(f) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("Gene" %in% colnames(d))
  m <- as.matrix(d[, setdiff(colnames(d), "Gene"), drop = FALSE])
  mode(m) <- "numeric"; rownames(m) <- make.unique(as.character(d$Gene))
  if (anyNA(m)) stop("Missing values in ", f)
  round(m)
}

# ------------------------------- 3. JOINT VST (tumors + PDOs) ----------------
# Same construction as the ssGSEA script: one VST on the raw counts of both sets (blind, intercept only).
cnt_t <- read_counts(CNT_T_FILE); colnames(cnt_t) <- norm_id(colnames(cnt_t))
cnt_p <- read_counts(CNT_P_FILE)
if (ncol(cnt_t) != 31) stop("Expected 31 tumors, found ", ncol(cnt_t))
if (ncol(cnt_p) != 4)  stop("Expected 4 PDOs, found ", ncol(cnt_p))
if (any(colnames(cnt_p) %in% colnames(cnt_t))) stop("A PDO name equals a tumor name.")
tum_keys <- colnames(cnt_t); pdo_keys <- colnames(cnt_p)

common_cnt <- intersect(rownames(cnt_t), rownames(cnt_p))
cnt <- cbind(cnt_t[common_cnt, , drop = FALSE], cnt_p[common_cnt, , drop = FALSE])
keep <- rowSums(cnt >= MIN_COUNT) >= ceiling(ncol(cnt) * MIN_FRAC)
message("Genes kept after the low-expression filter: ", sum(keep), " of ", nrow(cnt))
storage.mode(cnt) <- "integer"
dds <- DESeqDataSetFromMatrix(countData = cnt[keep, ],
                              colData = data.frame(sample = colnames(cnt), row.names = colnames(cnt)),
                              design = ~ 1)
vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) { message("vst() failed (", conditionMessage(e), "); using varianceStabilizingTransformation()")
                                      varianceStabilizingTransformation(dds, blind = TRUE) })
expr <- assay(vsd)

sym <- sub("\\|.*$", "", rownames(expr))                      # "SYMBOL|ENSG" -> "SYMBOL"
if (mean(grepl("^ENSG", sym)) > 0.5) stop("Row names look like Ensembl IDs; the count files should carry gene symbols.")
okr <- !is.na(sym) & sym != ""
expr <- expr[okr, ]; sym <- sym[okr]
o <- order(rowMeans(expr), decreasing = TRUE)
expr <- expr[o, ]; sym <- sym[o]
expr <- expr[!duplicated(sym), ]; rownames(expr) <- sym[!duplicated(sym)]
message(sprintf("Expression input: %d genes x %d samples", nrow(expr), ncol(expr)))

# ------------------------------- 4. Cluster labels of the tumors -------------
clus <- read.csv(CLUSTER_FILE); clus$Key <- norm_id(clus$sample)
cl <- factor(paste("Cluster", clus$cluster[match(tum_keys, clus$Key)]), levels = cluster_levels)
if (anyNA(cl)) stop("No cluster for: ", paste(tum_keys[is.na(cl)], collapse = ", "))
print(table(cl))
n_cl <- as.integer(table(cl))

# ------------------------------- 5. Read W -----------------------------------
W_raw <- fread(W_FILE, data.table = FALSE, check.names = FALSE)
message("W file: ", nrow(W_raw), " rows x ", ncol(W_raw), " columns. Columns: ", paste(colnames(W_raw), collapse = ", "))
is_num   <- vapply(W_raw, is.numeric, logical(1))
gene_col <- which(!is_num)[1]
if (is.na(gene_col)) stop("No text column with gene names found in W; first rows:\n", paste(capture.output(print(head(W_raw, 3))), collapse = "\n"))
mg_cols <- which(is_num)
if (length(mg_cols) != 3) stop("Expected 3 numeric metagene columns (k = 3), found ", length(mg_cols), ": ", paste(colnames(W_raw)[mg_cols], collapse = ", "))

w_genes <- sub("\\|.*$", "", as.character(W_raw[[gene_col]]))
if (mean(grepl("^ENSG", w_genes)) > 0.5) stop("Genes in W are Ensembl IDs; they must be mapped to gene symbols first.")
W <- as.matrix(W_raw[, mg_cols, drop = FALSE])
mg_names <- paste0("MG", 1:3)
mg_orig  <- colnames(W)
colnames(W) <- mg_names; rownames(W) <- w_genes
if (any(W < 0, na.rm = TRUE)) warning("W has negative values; NMF weights should be non-negative.")
if (anyNA(W)) stop("W has missing values.")
if (anyDuplicated(rownames(W))) {                           # duplicated gene names: keep the row with the largest total weight
  W <- W[order(-rowSums(W)), , drop = FALSE]
  W <- W[!duplicated(rownames(W)), , drop = FALSE]
}
message("Metagene columns in W: ", paste(mg_orig, collapse = ", "), " -> ", paste(mg_names, collapse = ", "))

common <- intersect(rownames(W), rownames(expr))
cov_w  <- length(common) / nrow(W)
message(sprintf("Genes of W found in the expression input: %d of %d (%.1f%%)", length(common), nrow(W), 100 * cov_w))
if (cov_w < 0.5) stop("Fewer than 50% of the W genes are in the expression input; check the gene naming.")
if (cov_w < 0.8) warning("Less than 80% of the W genes are in the expression input; the projection is less reliable.")

gsd    <- apply(expr[common, tum_keys, drop = FALSE], 1, sd)  # genes without tumor variance cannot be z-scored
common_var <- common[gsd > 0]
W_c <- W[common, , drop = FALSE]
X   <- expr[common, , drop = FALSE]

# ------------------------------- 6. A. NMF projection (NNLS) -----------------
H <- apply(X, 2, function(x) nnls::nnls(W_c, x)$x)           # 3 metagenes x samples
rownames(H) <- mg_names; colnames(H) <- colnames(X)
H_prop <- sweep(H, 2, colSums(H), "/")

# VST values are not centred, so every metagene carries a large common baseline and the raw
# proportions are nearly equal across samples. Loadings are therefore compared on the tumor scale.
tH_mean <- rowMeans(H[, tum_keys, drop = FALSE])
tH_sd   <- apply(H[, tum_keys, drop = FALSE], 1, sd)
Hz      <- (H - tH_mean) / tH_sd                              # SDs from the tumor mean

cat("\nNNLS loadings H (raw): range", signif(range(H), 3), "| SD across tumors per metagene:", signif(tH_sd, 3), "\n")
if (any(!is.finite(Hz))) stop("Standardised loadings contain NA/Inf: a metagene has zero variance across tumors or all-zero loadings. ",
                              "Check W and the expression input (SD per metagene printed above).")

# Metagene -> cluster map, learned from the tumors
map_mat <- sapply(cluster_levels, function(k) rowMeans(Hz[, tum_keys[cl == k], drop = FALSE]))
rownames(map_mat) <- mg_names
map_cl  <- cluster_levels[max.col(map_mat, ties.method = "first")]
one_to_one <- !anyDuplicated(map_cl)
mg_label <- if (one_to_one) paste0(mg_names, " (", map_cl, ")") else mg_names
names(mg_label) <- mg_names
if (!one_to_one) warning("Two metagenes map to the same cluster; labels are left as MG1-MG3. Check map_mat.")
cat("\nMean tumor loading (SDs) per cluster; rows = metagenes\n"); print(round(map_mat, 2))
cat("Metagene -> cluster:", paste(mg_names, map_cl, sep = " -> ", collapse = "; "), "\n")
write.csv(data.frame(Metagene = mg_names, W_column = mg_orig, Mapped_cluster = map_cl, map_mat, check.names = FALSE),
          file.path(OUT, "NMF_metagene_to_cluster_map.csv"), row.names = FALSE)

# Nearest centroid in the loading space; tumors are classified leave-one-out
cent_full <- map_mat                                          # metagenes x clusters
loo_call <- vapply(seq_along(tum_keys), function(i) {
  cent <- cent_full
  k    <- as.character(cl[i])
  idx  <- which(cl == k & seq_along(cl) != i)
  cent[, k] <- rowMeans(Hz[, tum_keys[idx], drop = FALSE])
  d <- sqrt(colSums((Hz[, tum_keys[i]] - cent)^2))
  names(which.min(d))
}, character(1))
acc  <- mean(loo_call == as.character(cl))
conf <- table(Given = cl, Projected = factor(loo_call, levels = cluster_levels))
cat(sprintf("\nValidation: %d of %d tumors (%.0f%%) return to their own cluster after projection (leave-one-out).\n",
            sum(loo_call == as.character(cl)), length(cl), 100 * acc))
print(conf)
if (acc < 0.8) warning("Less than 80% of tumors return to their own cluster. The expression input probably differs ",
                       "from the NMF input (scale or gene set); do not interpret the PDO placement.")
write.csv(as.data.frame.matrix(conf), file.path(OUT, "NMF_projection_validation_tumors_confusion.csv"))

# PDO assignment
d_pdo   <- sapply(pdo_keys, function(s) sqrt(colSums((Hz[, s] - cent_full)^2)))   # clusters x PDOs
nearest <- cluster_levels[apply(d_pdo, 2, which.min)]
margin  <- apply(d_pdo, 2, function(v) sort(v)[2] - sort(v)[1])

tr <- function(m, prefix, cols) { x <- t(m); colnames(x) <- paste0(prefix, cols); x }   # samples x features, prefixed names
pdo_tab <- data.frame(Sample = pdo_keys,
                      tr(H[, pdo_keys, drop = FALSE],      "H_",             mg_names),
                      tr(H_prop[, pdo_keys, drop = FALSE], "H_share_",       mg_names),
                      tr(Hz[, pdo_keys, drop = FALSE],     "H_z_vs_tumors_", mg_names),
                      tr(d_pdo,                            "Distance_to_",   gsub(" ", "", cluster_levels)),
                      Nearest_cluster = nearest, Distance_margin = margin,
                      row.names = NULL, check.names = FALSE)
write.csv(pdo_tab, file.path(OUT, "NMF_projection_PDO.csv"), row.names = FALSE)

tum_tab <- data.frame(Sample = tum_keys, Given_cluster = as.character(cl), LOO_projected_cluster = loo_call,
                      tr(Hz[, tum_keys, drop = FALSE], "H_z_vs_tumors_", mg_names),
                      row.names = NULL, check.names = FALSE)
write.csv(tum_tab, file.path(OUT, "NMF_projection_tumors.csv"), row.names = FALSE)
cat("\nPDO projection:\n"); print(pdo_tab[, c("Sample", paste0("H_z_vs_tumors_", mg_names), "Nearest_cluster", "Distance_margin")], digits = 3)

# Figure A: loadings (SDs from the tumor mean), tumors by cluster then PDOs
col_ord <- c(tum_keys[order(cl)], pdo_keys)
grp     <- factor(c(as.character(cl[order(cl)]), rep("PDO", length(pdo_keys))), levels = group_levels)
ann_col <- data.frame(Group = grp, row.names = col_ord)
Hz_plot <- Hz[, col_ord, drop = FALSE]; rownames(Hz_plot) <- mg_label[rownames(Hz_plot)]
Hz_plot[Hz_plot >  Z_CLIP] <-  Z_CLIP
Hz_plot[Hz_plot < -Z_CLIP] <- -Z_CLIP
draw_A <- function() {
  pheatmap(Hz_plot, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           breaks = seq(-Z_CLIP, Z_CLIP, length.out = 102),
           annotation_col = ann_col, annotation_colors = list(Group = group_colors),
           cluster_rows = FALSE, cluster_cols = FALSE, gaps_col = cumsum(n_cl),
           fontsize_row = 10, fontsize_col = 7, border_color = NA,
           main = "NMF metagene loadings (NNLS projection, SDs from the tumor mean)")
}
save_pheatmap(draw_A, file.path(OUT, "NMF_projection_loadings_heatmap"), 11, 5)

# ------------------------------- 7. B. Top-weighted genes --------------------
Wv <- W[common_var, , drop = FALSE]
best <- max.col(Wv, ties.method = "first")
top_list <- lapply(1:3, function(k) {
  g <- if (ASSIGN_UNIQUE) rownames(Wv)[best == k] else rownames(Wv)
  head(g[order(-Wv[g, k])], TOP_N)
})
names(top_list) <- mg_names

# Gene set for the PCA: the PCA_N genes with the highest weight in W overall.
# Each gene is ranked by its largest weight over the 3 metagenes; there is no split per metagene.
w_max     <- apply(Wv, 1, max)
pca_genes <- names(sort(w_max, decreasing = TRUE))[seq_len(min(PCA_N, length(w_max)))]
pca_src   <- factor(mg_names[best[match(pca_genes, rownames(Wv))]], levels = mg_names)   # metagene where the gene is highest
message("PCA gene set: ", length(pca_genes), " genes; from ", paste(mg_names, as.integer(table(pca_src)), sep = " = ", collapse = ", "))
fwrite(data.table(Rank = seq_along(pca_genes), Gene = pca_genes, Max_weight = w_max[pca_genes],
                  Metagene = as.character(pca_src), Mapped_cluster = map_cl[as.integer(pca_src)]),
       file.path(OUT, "NMF_PCA_weighted_genes.csv"))
top_tab <- rbindlist(lapply(1:3, function(k) data.table(
  Metagene = mg_names[k], Mapped_cluster = map_cl[k], Rank = seq_along(top_list[[k]]), Gene = top_list[[k]],
  Weight = Wv[top_list[[k]], k], Share_of_gene_weight = Wv[top_list[[k]], k] / rowSums(Wv[top_list[[k]], , drop = FALSE]))))
fwrite(top_tab, file.path(OUT, "NMF_top_weighted_genes.csv"))
message("Top genes per metagene: ", paste(lengths(top_list), collapse = ", "))

union_g <- unique(unlist(top_list))
t_mean  <- rowMeans(expr[union_g, tum_keys, drop = FALSE])
t_sd    <- apply(expr[union_g, tum_keys, drop = FALSE], 1, sd)
gz      <- (expr[union_g, , drop = FALSE] - t_mean) / t_sd     # gene z-scores on the tumor scale

scores <- sapply(1:3, function(k) colMeans(gz[top_list[[k]], , drop = FALSE]))
colnames(scores) <- mg_names
fwrite(data.table(Sample = rownames(scores), Group = ifelse(rownames(scores) %in% pdo_keys, "PDO",
                  as.character(cl)[match(rownames(scores), tum_keys)]), scores),
       file.path(OUT, "NMF_top_gene_signature_scores.csv"))

# Figure B1: heatmap of the top genes
mg_order  <- if (one_to_one) order(match(map_cl, cluster_levels)) else 1:3   # row blocks in the same order as the column clusters
tl_ord    <- top_list[mg_order]
genes_ord <- unname(unlist(tl_ord))
mat_h <- gz[genes_ord, col_ord, drop = FALSE]
rownames(mat_h) <- make.unique(genes_ord)
mat_h[mat_h >  Z_CLIP] <-  Z_CLIP
mat_h[mat_h < -Z_CLIP] <- -Z_CLIP
ann_row <- data.frame(Metagene = factor(rep(mg_label[mg_order], lengths(tl_ord)), levels = unname(mg_label[mg_order])),
                      row.names = rownames(mat_h))
mg_cols_col <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A"), unname(mg_label))
draw_B <- function() {
  pheatmap(mat_h, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           breaks = seq(-Z_CLIP, Z_CLIP, length.out = 102),
           annotation_col = ann_col, annotation_row = ann_row,
           annotation_colors = list(Group = group_colors, Metagene = mg_cols_col),
           cluster_rows = FALSE, cluster_cols = FALSE,
           gaps_col = cumsum(n_cl), gaps_row = head(cumsum(lengths(tl_ord)), -1),
           fontsize_row = 4.5, fontsize_col = 7, border_color = NA,
           main = paste0("Top ", TOP_N, " weighted genes per metagene (SDs from the tumor mean)"))
}
save_pheatmap(draw_B, file.path(OUT, "NMF_top_genes_heatmap"), 11, max(10, 0.075 * length(genes_ord) + 1.5))   # height grows with the number of genes

# Figure B2: signature scores, tumors by cluster with PDOs overlaid
sc_long <- rbindlist(lapply(1:3, function(k) data.table(
  Sample = rownames(scores), Score = scores[, k], Metagene = mg_label[k],
  Group = ifelse(rownames(scores) %in% pdo_keys, "PDO", as.character(cl)[match(rownames(scores), tum_keys)]))))
sc_long[, Metagene := factor(Metagene, levels = unname(mg_label))]
tum_sc <- sc_long[Group != "PDO"]; pdo_sc <- sc_long[Group == "PDO"]
pdo_cols <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[seq_along(pdo_keys)], pdo_keys)
p_sc <- ggplot() +
  geom_boxplot(data = tum_sc, aes(Group, Score, fill = Group), outlier.shape = NA, alpha = 0.5, colour = "black") +
  geom_point(data = tum_sc, aes(Group, Score), colour = "grey30", size = 1.2,
             position = position_jitter(width = 0.15, height = 0, seed = SEED)) +
  geom_point(data = pdo_sc, aes(Group, Score, colour = Sample), size = 3,
             position = position_jitter(width = 0.1, height = 0, seed = SEED)) +
  facet_wrap(~ Metagene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_colour_manual(values = pdo_cols, name = "PDO") +
  scale_x_discrete(limits = group_levels) +
  labs(x = NULL, y = paste0("Mean z-score of the top ", TOP_N, " genes (tumor scale)"),
       title = "NMF top-weighted gene signatures: tumors by cluster and PDOs",
       caption = "Boxes: tumors by RNA cluster; grey points: tumors; coloured points: PDOs. Descriptive (n = 4).") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        plot.caption = element_text(size = 9, hjust = 0))
save_plot(p_sc, file.path(OUT, "NMF_top_gene_signature_scores_boxplot"), 11, 4.5)

# ------------------------------- 8. C. PCA on the top-weighted genes ---------
# Axes are built from the tumors only; the PDOs are projected with the tumor centring/scaling.
pca  <- prcomp(t(expr[pca_genes, tum_keys, drop = FALSE]), center = TRUE, scale. = TRUE)
Xp   <- scale(t(expr[pca_genes, pdo_keys, drop = FALSE]), center = pca$center, scale = pca$scale)
proj <- Xp %*% pca$rotation
ve   <- pca$sdev^2 / sum(pca$sdev^2)
npc  <- min(5, ncol(pca$x))
pc_all <- rbind(pca$x[, 1:npc, drop = FALSE], proj[, 1:npc, drop = FALSE])
pc_df  <- data.frame(Sample = rownames(pc_all), pc_all, check.names = FALSE,
                     Group = c(as.character(cl), rep("PDO", length(pdo_keys))), row.names = NULL)
write.csv(pc_df, file.path(OUT, "NMF_PCA_weighted_genes_coordinates.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(ve)), Variance_explained = ve),
          file.path(OUT, "NMF_PCA_weighted_genes_variance.csv"), row.names = FALSE)

px <- paste0("PC", PCS[1]); py <- paste0("PC", PCS[2])
p_pca <- ggplot() +
  geom_point(data = pc_df[pc_df$Group != "PDO", ], aes(.data[[px]], .data[[py]], colour = Group), size = 2.8, alpha = 0.9) +
  geom_point(data = pc_df[pc_df$Group == "PDO", ], aes(.data[[px]], .data[[py]]),
             shape = 24, size = 3, fill = alpha("black", 0.6), colour = "white", stroke = 0.6) +   # smaller, semi-transparent: tumors underneath stay visible
  scale_colour_manual(values = cluster_colors, name = "Tumor cluster") +
  labs(x = sprintf("PC%d (%.1f%% of tumor variance)", PCS[1], 100 * ve[PCS[1]]),
       y = sprintf("PC%d (%.1f%% of tumor variance)", PCS[2], 100 * ve[PCS[2]]),
       title = paste0("PCA on the top ", length(pca_genes), " NMF-weighted genes"),
       caption = "Axes from the 31 tumors only; PDOs (semi-transparent black triangles) are projected. Descriptive (n = 4).") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.caption = element_text(size = 9, hjust = 0))
p_pca <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_pca + ggrepel::geom_text_repel(data = pc_df[pc_df$Group == "PDO", ], aes(.data[[px]], .data[[py]], label = Sample),
                                   size = 3.4, fontface = "bold", point.padding = 0.5, box.padding = 0.6,
                                   min.segment.length = 0.3, seed = SEED)   # labels kept off the markers
} else {
  p_pca + geom_text(data = pc_df[pc_df$Group == "PDO", ], aes(.data[[px]], .data[[py]], label = Sample),
                    size = 3.4, fontface = "bold", vjust = -1.4)
}
save_plot(p_pca, file.path(OUT, "NMF_PCA_weighted_genes_tumors_with_PDO_projected"), 7.5, 6)

# ------------------------------- 9. D. Correlation with cluster centroids ----
cen <- sapply(cluster_levels, function(k) rowMeans(gz[, tum_keys[cl == k], drop = FALSE]))
cor_pdo <- t(sapply(pdo_keys, function(s) sapply(cluster_levels, function(k) cor(gz[, s], cen[, k], method = "spearman"))))
cor_loo <- t(sapply(seq_along(tum_keys), function(i) {
  sapply(cluster_levels, function(k) {
    idx <- which(cl == k & seq_along(cl) != i)
    cor(gz[, tum_keys[i]], rowMeans(gz[, tum_keys[idx], drop = FALSE]), method = "spearman")
  })
}))
cor_call <- cluster_levels[max.col(cor_loo, ties.method = "first")]
acc_cor  <- mean(cor_call == as.character(cl))
cat(sprintf("\nCentroid correlation (top genes), tumors leave-one-out: %.0f%% return to their own cluster.\n", 100 * acc_cor))
cor_tab <- data.frame(Sample = pdo_keys, cor_pdo, Highest_correlation = cluster_levels[max.col(cor_pdo, ties.method = "first")],
                      row.names = NULL, check.names = FALSE)
colnames(cor_tab)[2:4] <- paste0("Spearman_to_", gsub(" ", "", cluster_levels), "_centroid")
write.csv(cor_tab, file.path(OUT, "NMF_top_genes_centroid_correlation_PDO.csv"), row.names = FALSE)
cat("\nPDO correlation to tumor cluster centroids:\n"); print(cor_tab, digits = 3)

message("\nDone. Outputs in: ", OUT)