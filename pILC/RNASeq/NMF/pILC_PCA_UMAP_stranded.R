# pILC RNA-seq dimensionality plots (stranded / corrected counts)
# Outputs: PCA (cluster), PCA (cluster + grade), UMAP (cluster) as PNG + PDF,
#          plus tables with PCA scores, variance explained and UMAP coordinates.

suppressPackageStartupMessages({
  library(ggplot2)
  library(uwot)
})

SEED <- 42
set.seed(SEED)

# ---------------------------------------------------------------- 1. Paths
root    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"
vst_f   <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_VST_matrix_all_expression_filtered_genes.rds")
clus_f  <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")  # 50-run clusters, GenePattern numbering (columns: sample, cluster)
if (!file.exists(clus_f)) stop("Cluster file not found: ", clus_f)
clin_f  <- file.path(root, "Clinical_correlation/pILC_ClinicoPathological_details.csv")
out_dir <- file.path(root, "RNASeq/Dimensionality_plot/Stranded")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- 2. Settings
N_TOP_GENES <- 2000  # most variable genes used for PCA/UMAP (500 = DESeq2 plotPCA default)
N_NEIGHBORS <- 10    # UMAP neighbours (small cohort, n ~ 31)
clus_col    <- "cluster"   # columns in Cluster_membership_k3.csv: sample, cluster
clus_id_col <- "sample"
grade_col   <- "Grade/Differentiation-Desc"   # column in the clinical file
clin_id_col <- "Sample"
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")  # Dark2
grade_shapes <- c("II" = 16, "III" = 17, "Unknown" = 4)

# Harmonize sample IDs across assays (e.g. TP09.9_S34 / TP09-9 / TP09-09T -> TP09-09)
norm_id <- function(x) {
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("\\.", "-", x)
  x <- sub("^(TP[0-9]+)-([0-9]+).*$", "\\1-\\2", x)
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}

# ---------------------------------------------------------------- 3. Load data
vst <- as.matrix(readRDS(vst_f))   # expected: genes x samples
if (nrow(vst) < ncol(vst)) vst <- t(vst)   # safety: make sure genes are in rows
if (is.null(colnames(vst))) stop("VST matrix has no sample names; take them from VST_sample_names.csv")
clus <- read.csv(clus_f, check.names = FALSE)
clin <- read.csv(clin_f, check.names = FALSE)

# Show the real column names so a wrong guess is easy to fix
message("Cluster file columns: ", paste(colnames(clus), collapse = ", "))
message("Clinical file columns: ", paste(colnames(clin), collapse = ", "))
if (!clus_id_col %in% colnames(clus)) {
  message("'", clus_id_col, "' not found in cluster file; using first column '", colnames(clus)[1], "'")
  clus_id_col <- colnames(clus)[1]
}
if (!clus_col %in% colnames(clus)) stop("Set clus_col to one of: ", paste(colnames(clus), collapse = ", "))
if (!clin_id_col %in% colnames(clin)) stop("Set clin_id_col to one of: ", paste(colnames(clin), collapse = ", "))
if (!grade_col %in% colnames(clin)) stop("Set grade_col to one of: ", paste(colnames(clin), collapse = ", "))

meta <- data.frame(RNA_ID = colnames(vst), Key = norm_id(colnames(vst)))
clus$Key <- norm_id(clus[[clus_id_col]])
clin$Key <- norm_id(clin[[clin_id_col]])
meta$Cluster <- factor(clus[[clus_col]][match(meta$Key, clus$Key)], levels = names(cluster_cols))
g_raw <- as.character(clin[[grade_col]][match(meta$Key, clin$Key)])
message("Raw grade values found:")
print(table(g_raw, useNA = "ifany"))
# File values are "Grade II", "Grade III" or "Cell type not determined" (-> Unknown)
g_raw <- trimws(g_raw)
g <- rep("Unknown", length(g_raw))
g[g_raw %in% "Grade II"]  <- "II"
g[g_raw %in% "Grade III"] <- "III"
meta$Grade <- factor(g, levels = names(grade_shapes))

# Safety check: every sample must have a cluster
if (anyNA(meta$Cluster)) {
  stop("Samples without cluster: ", paste(meta$RNA_ID[is.na(meta$Cluster)], collapse = ", "))
}

# ---------------------------------------------------------------- 4. Gene selection
vars <- apply(vst, 1, var)
top  <- names(sort(vars, decreasing = TRUE))[seq_len(min(N_TOP_GENES, length(vars)))]
mat  <- t(vst[top, ])   # samples x genes

# ---------------------------------------------------------------- 5. PCA
pca <- prcomp(mat, center = TRUE, scale. = FALSE)  # VST is already variance-stabilized
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pc_df <- cbind(meta, PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# ---------------------------------------------------------------- 6. UMAP
um <- umap(mat, n_neighbors = N_NEIGHBORS, min_dist = 0.3, metric = "euclidean",
           n_components = 2, seed = SEED)  # seed -> reproducible layout
um_df <- cbind(meta, UMAP1 = um[, 1], UMAP2 = um[, 2])

# ---------------------------------------------------------------- 7. Plot theme
base_theme <- theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20),
        axis.line = element_line(linewidth = 0.8),
        legend.title = element_text(size = 16))

p_pca <- ggplot(pc_df, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 4, alpha = 0.95) +
  scale_color_manual(values = cluster_cols) +
  labs(title = "PCA of pILC Samples (vst_transformed)") + base_theme

p_pca_grade <- ggplot(pc_df, aes(PC1, PC2, color = Cluster, shape = Grade)) +
  geom_point(size = 4, alpha = 0.95, stroke = 1) +
  scale_color_manual(values = cluster_cols) +
  scale_shape_manual(values = grade_shapes) +
  labs(title = "PCA of pILC Samples by Cluster and Grade") + base_theme

p_umap <- ggplot(um_df, aes(UMAP1, UMAP2, color = Cluster)) +
  geom_point(size = 4, alpha = 0.95) +
  scale_color_manual(values = cluster_cols) +
  labs(title = "UMAP of pILC Samples (vst_transformed)") + base_theme

# ---------------------------------------------------------------- 8. Save (PNG + PDF, no dates in names)
save_both <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
}
save_both(p_pca,       "PCA_pILC_vst_Cluster")
save_both(p_pca_grade, "PCA_pILC_vst_Cluster_Grade", w = 8.5)
save_both(p_umap,      "UMAP_pILC_vst_Cluster")

# ---------------------------------------------------------------- 9. Tables paired with the figures
write.csv(cbind(pc_df, PC3 = pca$x[, 3], PC4 = pca$x[, 4]),
          file.path(out_dir, "PCA_scores_pILC_vst.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(pve)), Variance_explained_pct = pve),
          file.path(out_dir, "PCA_variance_explained_pILC_vst.csv"), row.names = FALSE)
write.csv(um_df, file.path(out_dir, "UMAP_coordinates_pILC_vst.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 10. NMF consensus: which cluster 3 samples lean toward cluster 2?
# Consensus matrix = fraction of NMF runs in which two samples fall in the same cluster (1 = always together)
# OPTIONAL extra analysis (not needed for PCA_pILC_vst_Cluster / UMAP_pILC_vst_Cluster). It reads the 20-RUN consensus
# matrix, so it only matches the 20-run labels; keep FALSE when clusters come from the 50-run labels above.
RUN_BOUNDARY_CHECK <- FALSE
if (RUN_BOUNDARY_CHECK) {
suppressPackageStartupMessages(library(ggrepel))
cons_f <- file.path(root, "RNASeq/NMF/R_NMFConsensus_20runs/Consensus_matrix_k3.csv")
cons <- as.matrix(read.csv(cons_f, row.names = 1, check.names = FALSE))
rownames(cons) <- norm_id(rownames(cons)); colnames(cons) <- norm_id(colnames(cons))
stopifnot(all(meta$Key %in% rownames(cons)))
cons <- cons[meta$Key, meta$Key]
diag(cons) <- NA   # ignore the sample's consensus with itself

# Mean consensus of every sample with the members of each cluster (itself excluded)
mean_cons <- sapply(levels(meta$Cluster), function(k)
  rowMeans(cons[, meta$Cluster == k, drop = FALSE], na.rm = TRUE))
colnames(mean_cons) <- paste0("mean_consensus_cluster", levels(meta$Cluster))

idx_own <- cbind(seq_len(nrow(mean_cons)), as.integer(as.character(meta$Cluster)))
own     <- mean_cons[idx_own]
m       <- mean_cons
m[idx_own] <- -Inf                                   # hide own cluster
best_other     <- max.col(m, ties.method = "first")  # column index = cluster number
best_other_val <- m[cbind(seq_len(nrow(m)), best_other)]

cons_tab <- data.frame(
  Sample = meta$RNA_ID, Key = meta$Key, Cluster = meta$Cluster, Grade = meta$Grade,
  PC1 = pc_df$PC1, PC2 = pc_df$PC2, mean_cons,
  own_cluster_consensus     = own,
  strongest_other_cluster   = best_other,
  strongest_other_consensus = best_other_val,
  other_minus_own           = best_other_val - own,
  # Ambiguous: another cluster's consensus is within 0.1 of (or above) the sample's own cluster
  Ambiguous = (best_other_val - own) > -0.1
)
cons_tab <- cons_tab[order(cons_tab$Cluster, -cons_tab$other_minus_own), ]
write.csv(cons_tab, file.path(out_dir, "Sample_NMF_consensus_by_cluster_pILC.csv"), row.names = FALSE)

c3_to_c2 <- subset(cons_tab, Cluster == "3" & Ambiguous & strongest_other_cluster == 2)
print(c3_to_c2[, c("Key", "own_cluster_consensus", "mean_consensus_cluster2", "other_minus_own")])
write.csv(c3_to_c2, file.path(out_dir, "Cluster3_samples_near_cluster2_pILC.csv"), row.names = FALSE)

# PCA with the ambiguous cluster 3 samples labelled
pc_df$Label <- ifelse(pc_df$RNA_ID %in% c3_to_c2$Sample, pc_df$Key, NA)
p_lab <- ggplot(pc_df, aes(PC1, PC2, color = Cluster, shape = Grade)) +
  geom_point(size = 4, alpha = 0.95, stroke = 1) +
  geom_text_repel(aes(label = Label), color = "black", size = 3.5, na.rm = TRUE,
                  max.overlaps = Inf, seed = SEED) +
  scale_color_manual(values = cluster_cols) +
  scale_shape_manual(values = grade_shapes) +
  labs(title = "PCA: cluster 3 samples leaning toward cluster 2 (NMF consensus)") + base_theme
save_both(p_lab, "PCA_pILC_vst_Cluster3_near_Cluster2_labeled", w = 9.5)
}   # end of RUN_BOUNDARY_CHECK

message("Done. Variance explained: PC1 = ", round(pve[1], 1), "%, PC2 = ", round(pve[2], 1), "%")
