# ============================================================
# Top cancer-epithelial genes from the PBCA InstaPrism reference + heatmap in the pILC bulk RNA-seq
#
# READS saved files only (reference .Rds, TPM, cluster labels, InstaPrism fractions, ESTIMATE scores).
# It does NOT rerun InstaPrism and does NOT overwrite any existing file; new files go to Downstream_stats/Cancer_epithelial_genes/.
#
# Steps
#   1. Build one "cancer epithelial" profile = mean of the cancer-cell profiles of the reference (phi.cs)
#   2. Rank genes by SPECIFICITY: log2(cancer profile / highest value in any non-epithelial cell type)
#   3. Keep genes that are well expressed in the cancer profile and in the bulk TPM; take the top N
#   4. Heatmap of those genes (log2 TPM+1, z-score per gene), samples grouped by GP-labelled cluster
#   5. Sanity check: signature score vs InstaPrism cancer fraction and ESTIMATE tumour purity
# ============================================================

suppressPackageStartupMessages({
  library(InstaPrism)   # only so readRDS() can rebuild the S4 reference object
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)   # heatmap with thin annotation strips (same style as the Hallmark ssGSEA heatmap)
  library(circlize)         # colorRamp2()
  library(grid)
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths and settings ----------------
RNA_DIR   <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
INSTA_DIR <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/InstaPrism_PBCA_reference/Stranded")
OUT_DIR   <- file.path(INSTA_DIR, "Downstream_stats", "Cancer_epithelial_genes")   # own subfolder, not mixed with the other outputs
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PBCA_REF_FILE <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/PBCA_InstaPrism_reference.Rds")
TPM_FILE      <- file.path(RNA_DIR, "4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt")
CLUSTER_FILE  <- file.path(RNA_DIR, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")   # GP numbering
FRAC_FILE     <- file.path(INSTA_DIR, "pILC_InstaPrism_PBCA_fraction_with_summary.csv")
ESTIMATE_FILE <- file.path(RNA_DIR, "immunedeconv/Stranded/TPM/ESTIMATE_scores.csv")

N_TOP                <- 100    # genes shown in the heatmap (50 / 100 / 150 are sensible; labels get smaller above 60)
ROW_ORDER            <- "cluster_mean"   # "cluster_mean" = genes grouped by the cluster where their mean z-score is highest
                                         #                  (Cluster 1 group, then 2, then 3; strongest gene first inside each group)
                                         # "cluster"      = genes ordered by hierarchical clustering of their expression across samples
                                         # "rank"         = genes ordered by specificity rank (rank 1 at the top)
SHOW_RANK_IN_LABEL   <- TRUE        # TRUE = row labels like "KRT19 (#3)", where #3 is the specificity rank
MIN_TPM_MEAN         <- 5      # gene must have mean TPM >= this across the pILC samples
MIN_CANCER_QUANTILE  <- 0.75   # gene must be in the top 25% of genes by value in the cancer profile
PSEUDO               <- 1e-6   # pseudocount added to reference values before the log2 ratio
EXCLUDE_EPITHELIAL_FROM_OTHERS <- TRUE   # TRUE = compare cancer vs NON-epithelial cells (gives epithelial-lineage genes)
                                         # FALSE = also compare against normal/myoepithelial cells (cancer-specific genes)

CLUSTER_LEVELS <- c("1", "2", "3")
cluster_cols   <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")   # same Dark2 colours as the other figures

# The 17 named (non-cancer) cell types of the PBCA reference; every other column is a cancer-cell profile
NAMED_TYPES <- c("CD8+ T Cells", "CD4+ T Cells", "Regulatory T Cells", "Fibroblasts", "Epithelial Cells",
                 "Endothelial Cells", "Perivascular-like (PVL) Cells", "Macrophages", "NK Cells",
                 "Dendritic Cells", "MDSCs", "Myoepithelial Cells", "B Cells", "Plasma Cells",
                 "Monocytes", "Mast Cells", "Neutrophils")
EPITHELIAL_TYPES <- c("Epithelial Cells", "Myoepithelial Cells")

# ---------------- 2. Helpers ----------------
need <- function(path, what) if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)

# Same sample-ID harmonisation as the other scripts: TP09-9 / TP09-09 / TP09.9_S34 / TP09-09T -> TP09_9
make_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- sub("([0-9])T$", "\\1", x)
  gsub("_0+([0-9])", "_\\1", x)
}

fmt_p <- function(p) ifelse(is.na(p), "n/a",
                            ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3))))


for (f in c(PBCA_REF_FILE, TPM_FILE, CLUSTER_FILE)) need(f, "Input file")

# ---------------- 3. Cancer-epithelial profile from the reference ----------------
ref <- readRDS(PBCA_REF_FILE)
phi <- as.matrix(ref@phi.cs)                       # genes x profiles (cell types + cancer-cell profiles)
message("Reference: ", nrow(phi), " genes x ", ncol(phi), " profiles")
if (mean(grepl("^ENSG", rownames(phi))) > 0.5)
  stop("Reference row names look like Ensembl IDs; gene symbols are needed to match the TPM matrix.")
message("Column sums of the reference (first 5): ", paste(round(colSums(phi)[1:5], 3), collapse = ", "),
        "  (values are each gene's share of the profile)")

missing_named <- setdiff(NAMED_TYPES, colnames(phi))
if (length(missing_named) > 0) stop("Expected cell types not in the reference: ", paste(missing_named, collapse = ", "))

cancer_cols <- setdiff(colnames(phi), NAMED_TYPES)
other_cols  <- if (EXCLUDE_EPITHELIAL_FROM_OTHERS) setdiff(NAMED_TYPES, EPITHELIAL_TYPES) else NAMED_TYPES
message("Cancer-cell profiles used: ", length(cancer_cols), " (expected about 61, including 'Other')")
message("Non-cancer cell types used as comparison: ", length(other_cols))
if (length(cancer_cols) < 3) stop("Too few cancer-cell profiles found; check the column names printed above.")

cancer_mean <- rowMeans(phi[, cancer_cols, drop = FALSE])
other_max   <- apply(phi[, other_cols, drop = FALSE], 1, max)
other_which <- other_cols[apply(phi[, other_cols, drop = FALSE], 1, which.max)]   # cell type with the highest value

rank_tab <- data.frame(
  gene                 = rownames(phi),
  cancer_mean_value    = cancer_mean,
  highest_other_value  = other_max,
  highest_other_type   = other_which,
  log2_specificity     = log2((cancer_mean + PSEUDO) / (other_max + PSEUDO)),
  stringsAsFactors = FALSE
)
cancer_cutoff <- quantile(rank_tab$cancer_mean_value, MIN_CANCER_QUANTILE)

# ---------------- 4. Bulk TPM, clusters, candidate genes ----------------
tpm_dt <- fread(TPM_FILE, data.table = FALSE)
tpm <- as.matrix(tpm_dt[, -1]); rownames(tpm) <- tpm_dt[[1]]; storage.mode(tpm) <- "numeric"
colnames(tpm) <- make_key(colnames(tpm))
if (anyDuplicated(colnames(tpm))) stop("Duplicated sample keys in the TPM matrix.")

cl_raw <- fread(CLUSTER_FILE, data.table = FALSE)
message("Cluster file columns: ", paste(colnames(cl_raw), collapse = " | "))
s_hit <- grep("^sample$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
s_col <- if (length(s_hit) > 0) s_hit[1] else colnames(cl_raw)[1]
c_hit <- grep("^cluster$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(cl_raw)[ncol(cl_raw)]
cl <- data.frame(key = make_key(cl_raw[[s_col]]),
                 num = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cl_raw[[c_col]])))))
if (anyNA(cl$num) || !all(cl$num %in% 1:3)) stop("Cluster values must be 1, 2 or 3 in ", CLUSTER_FILE)
if (anyDuplicated(cl$key)) stop("Duplicated sample keys in the cluster file.")
cl$Cluster <- factor(as.character(cl$num), levels = CLUSTER_LEVELS)

samples <- intersect(colnames(tpm), cl$key)
message("Samples matched between TPM and clusters: ", length(samples), " of ", ncol(tpm))
if (length(samples) < 25) stop("Too few matched samples. Check the sample IDs.")
tpm <- tpm[, samples, drop = FALSE]
cl  <- cl[match(samples, cl$key), ]
message("Samples per cluster:"); print(table(cl$Cluster))

shared <- intersect(rank_tab$gene, rownames(tpm))
message("Genes shared by the reference and the TPM matrix: ", length(shared))
if (length(shared) < 5000) stop("Very few shared genes; gene symbols probably do not match.")

rank_tab$mean_TPM_pILC <- rowMeans(tpm)[match(rank_tab$gene, rownames(tpm))]
rank_tab$passes_filters <- !is.na(rank_tab$mean_TPM_pILC) &
  rank_tab$mean_TPM_pILC >= MIN_TPM_MEAN & rank_tab$cancer_mean_value >= cancer_cutoff
cand <- rank_tab[rank_tab$passes_filters, ]
cand <- cand[order(-cand$log2_specificity), ]
cand$rank <- seq_len(nrow(cand))
message("Genes passing the filters: ", nrow(cand))
if (nrow(cand) < N_TOP) stop("Fewer genes pass the filters than N_TOP; relax MIN_TPM_MEAN or MIN_CANCER_QUANTILE.")

write.csv(cand, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_gene_ranking.csv"), row.names = FALSE)
top_genes <- cand$gene[seq_len(N_TOP)]
message("Top genes: ", paste(top_genes, collapse = ", "))

# ---------------- 5. Heatmap data ----------------
expr <- log2(tpm[top_genes, , drop = FALSE] + 1)
Z <- t(scale(t(expr)))                              # z-score per gene across the samples
Z[is.na(Z)] <- 0                                    # a gene with zero variance would give NaN
score <- colMeans(Z)                                # signature score = mean z-score of the top genes

# In which cluster is each top gene highest? (mean z-score per cluster; the gene goes to the cluster with the largest mean)
zc <- sapply(CLUSTER_LEVELS, function(k) rowMeans(Z[, cl$Cluster == k, drop = FALSE]))
colnames(zc) <- paste0("mean_z_cluster_", CLUSTER_LEVELS)
highest <- CLUSTER_LEVELS[max.col(zc, ties.method = "first")]
names(highest) <- rownames(zc)
cand$cluster_highest <- highest[cand$gene]                       # NA for genes outside the top N
cand[paste0("mean_z_cluster_", CLUSTER_LEVELS)] <- zc[match(cand$gene, rownames(zc)), , drop = FALSE]
by_cl <- data.frame(Cluster = CLUSTER_LEVELS,
                    n_top_genes = as.integer(table(factor(highest, levels = CLUSTER_LEVELS))))
by_cl$percent_of_top_genes <- round(100 * by_cl$n_top_genes / sum(by_cl$n_top_genes), 1)
write.csv(by_cl, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_genes_by_cluster_highest.csv"), row.names = FALSE)
message("Top ", N_TOP, " genes by the cluster where their mean z-score is highest:"); print(by_cl)

# InstaPrism cancer fraction and ESTIMATE tumour purity for the same samples (optional, for the checks)
cancer_frac <- rep(NA_real_, length(samples)); purity <- rep(NA_real_, length(samples))
if (file.exists(FRAC_FILE)) {
  fr <- fread(FRAC_FILE, data.table = FALSE)
  id_col <- if ("Sample_ID" %in% colnames(fr)) "Sample_ID" else colnames(fr)[1]
  if ("Cancer Epithelial Cells" %in% colnames(fr))
    cancer_frac <- fr[["Cancer Epithelial Cells"]][match(samples, make_key(fr[[id_col]]))]
} else message("InstaPrism fraction file not found; cancer-fraction check skipped.")
if (file.exists(ESTIMATE_FILE)) {
  es <- fread(ESTIMATE_FILE, data.table = FALSE)
  if (all(c("Sample", "TumorPurity") %in% colnames(es)))
    purity <- es$TumorPurity[match(samples, make_key(es$Sample))]
} else message("ESTIMATE file not found; purity check skipped.")

ann <- data.frame(Sample = samples, Cluster = cl$Cluster, Signature_score = score,
                  InstaPrism_cancer_fraction = cancer_frac, ESTIMATE_TumorPurity = purity)

# Sample order: by cluster, then by signature score (high to left); gene order: hierarchical clustering
ann <- ann[order(ann$Cluster, -ann$Signature_score), ]
ann$Sample <- factor(ann$Sample, levels = ann$Sample)
# Row order: by cluster of highest mean z-score, by hierarchical clustering (1 - Pearson correlation, average linkage), or by rank
gene_order <- if (ROW_ORDER == "rank") {
  top_genes
} else if (ROW_ORDER == "cluster_mean") {
  rownames(zc)[order(match(highest, CLUSTER_LEVELS), -apply(zc, 1, max))]   # group by cluster of highest mean z, strongest first
} else {
  hc_rows <- hclust(as.dist(1 - cor(t(Z))), method = "average")      # kept, so the figure uses exactly this tree
  rownames(Z)[hc_rows$order]
}
# Save the final row position of each heatmap gene next to its specificity rank
cand$heatmap_row <- match(cand$gene, gene_order)          # 1 = top row; NA = gene not in the heatmap
write.csv(cand, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_gene_ranking.csv"), row.names = FALSE)

# ---------------- 6. Checks: signature score vs cluster, cancer fraction, purity ----------------
kw_p <- kruskal.test(Signature_score ~ Cluster, data = ann)$p.value
sp <- function(y) {
  ok <- is.finite(ann$Signature_score) & is.finite(y)
  if (sum(ok) < 8) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(ann$Signature_score[ok], y[ok], method = "spearman", exact = FALSE))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
r1 <- sp(ann$InstaPrism_cancer_fraction); r2 <- sp(ann$ESTIMATE_TumorPurity)
summary_tab <- data.frame(
  test = c("Kruskal-Wallis: signature score across clusters",
           "Spearman: signature score vs InstaPrism cancer epithelial fraction",
           "Spearman: signature score vs ESTIMATE tumour purity"),
  statistic = c(NA, r1["rho"], r2["rho"]),
  p_value   = c(kw_p, r1["p"], r2["p"]),
  n         = c(nrow(ann), r1["n"], r2["n"]),
  row.names = NULL)
write.csv(summary_tab, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_signature_summary.csv"), row.names = FALSE)
write.csv(ann, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_signature_scores.csv"), row.names = FALSE)
print(summary_tab)

# ---------------- 6b. Does a longer gene list change the separation between clusters? ----------------
# The gene lists are always the top-n of the SAME reference ranking (chosen without looking at your clusters).
N_GRID <- c(10, 25, 50, 100, 200, 300)
N_GRID <- N_GRID[N_GRID <= nrow(cand)]
sens <- do.call(rbind, lapply(N_GRID, function(n) {
  g  <- cand$gene[seq_len(n)]
  zz <- t(scale(t(log2(tpm[g, , drop = FALSE] + 1)))); zz[is.na(zz)] <- 0
  d  <- data.frame(sc = colMeans(zz), Cluster = cl$Cluster, cf = cancer_frac)
  ok <- is.finite(d$cf)
  data.frame(n_genes = n,
             KW_p_across_clusters   = kruskal.test(sc ~ Cluster, data = d)$p.value,
             rho_vs_cancer_fraction = if (sum(ok) >= 8) suppressWarnings(cor(d$sc[ok], d$cf[ok], method = "spearman")) else NA_real_)
}))
write.csv(sens, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_signature_by_gene_number.csv"), row.names = FALSE)
message("Signature score vs number of genes:"); print(sens)

# ---------------- 7. Figure (ComplexHeatmap) ----------------
row_font <- if (N_TOP <= 60) 9 else if (N_TOP <= 120) 6 else 4.5
n_by_cl  <- table(ann$Cluster)
col_titles <- paste0("Cluster ", CLUSTER_LEVELS, " (n = ", n_by_cl[CLUSTER_LEVELS], ")")

# Matrix in the final order: rows = genes (gene_order), columns = samples (cluster, then signature score)
mat <- if (ROW_ORDER == "cluster") Z[, as.character(ann$Sample), drop = FALSE] else Z[gene_order, as.character(ann$Sample), drop = FALSE]
row_labs <- if (SHOW_RANK_IN_LABEL) paste0(rownames(mat), " (#", cand$rank[match(rownames(mat), cand$gene)], ")") else rownames(mat)
rownames(mat) <- row_labs

# Thin annotation strips above the columns: RNA cluster and InstaPrism cancer fraction
top_ann <- HeatmapAnnotation(
  Cluster = ann$Cluster,
  `Cancer fraction` = ann$InstaPrism_cancer_fraction,
  col = list(Cluster = cluster_cols,
             `Cancer fraction` = colorRamp2(c(0, 0.5, 1), c("white", "#74C476", "#00441B"))),
  na_col = "grey85",
  simple_anno_size = unit(4, "mm"),
  gap = unit(1, "mm"),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontface = "bold", fontsize = 11),
  annotation_legend_param = list(Cluster = list(title = "Cluster", at = CLUSTER_LEVELS, labels = paste("Cluster", CLUSTER_LEVELS)),
                                 `Cancer fraction` = list(title = "InstaPrism\ncancer fraction")))

# Row groups: by the cluster of highest mean z-score (only for ROW_ORDER = "cluster_mean")
row_split <- NULL
if (ROW_ORDER == "cluster_mean") {
  grp_lab <- paste0("Highest in\nCluster ", by_cl$Cluster, "\n", by_cl$n_top_genes, " genes (", by_cl$percent_of_top_genes, "%)")
  row_split <- factor(highest[gene_order], levels = by_cl$Cluster, labels = grp_lab)
}

# Row dendrogram only for ROW_ORDER = "cluster" (same tree that defined gene_order, so the CSV positions match the figure)
row_clust <- if (ROW_ORDER == "cluster") hc_rows else FALSE

ht <- Heatmap(
  mat,
  name = "z-score of\nlog2(TPM+1)",
  col = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
  cluster_rows = row_clust, row_dend_side = "left", row_dend_width = unit(15, "mm"),
  cluster_columns = FALSE,
  column_split = factor(as.character(ann$Cluster), levels = CLUSTER_LEVELS), column_title = col_titles,
  column_title_gp = gpar(fontface = "bold", fontsize = 12), column_gap = unit(2, "mm"),
  row_split = row_split, row_title = NULL, row_gap = unit(2, "mm"),     # gaps between the gene blocks, no block labels on the left
  row_names_side = "right", row_names_gp = gpar(fontsize = row_font),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_ann,
  border = FALSE, use_raster = FALSE)

row_txt <- switch(ROW_ORDER,
                  rank = "specificity rank",
                  cluster_mean = "cluster of highest mean z-score (strongest first within each group)",
                  "hierarchical clustering")
fig_title <- paste0("Top ", N_TOP, " cancer-epithelial genes from the PBCA reference in the pILC samples")
fig_cap <- paste0("Row order: ", row_txt, "; (#n) = specificity rank.",
                  if (ROW_ORDER == "cluster_mean") paste0(" Gene blocks, top to bottom: highest in Cluster ", paste0(by_cl$Cluster, " (", by_cl$n_top_genes, " genes)", collapse = ", "), ".") else "",
                  "\n",
                  "Genes ranked by log2(cancer-cell profile / highest non-", ifelse(EXCLUDE_EPITHELIAL_FROM_OTHERS, "epithelial", "cancer"),
                  " cell type); mean TPM >= ", MIN_TPM_MEAN, ". Samples ordered by signature score within cluster.\n",
                  "Signature score vs InstaPrism cancer fraction: Spearman rho = ", formatC(r1["rho"], format = "f", digits = 2),
                  ", p ", fmt_p(r1["p"]), "; across clusters: Kruskal-Wallis p ", fmt_p(kw_p), ".")

draw_fig <- function() {
  draw(ht, column_title = fig_title, column_title_gp = gpar(fontface = "bold", fontsize = 15),
       heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
       padding = unit(c(22, 4, 4, 4), "mm"))                       # extra space at the bottom for the caption
  grid.text(fig_cap, x = unit(4, "mm"), y = unit(10, "mm"), just = c("left", "center"), gp = gpar(fontsize = 9))
}
fig_name <- file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cancer_epithelial_top_genes_heatmap")
fig_w <- 13
fig_h <- if (N_TOP <= 60) 0.17 * N_TOP + 4.5 else 0.1 * N_TOP + 4.5
png(paste0(fig_name, ".png"), width = fig_w, height = fig_h, units = "in", res = 300); draw_fig(); invisible(dev.off())
pdf(paste0(fig_name, ".pdf"), width = fig_w, height = fig_h);                          draw_fig(); invisible(dev.off())

message("Done. Files written to: ", OUT_DIR)