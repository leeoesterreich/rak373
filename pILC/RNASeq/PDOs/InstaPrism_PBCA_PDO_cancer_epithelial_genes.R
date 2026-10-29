# ============================================================
# Top cancer-epithelial genes (PBCA reference) in the 4 PDOs, shown next to the pILC tumors
#
# READS saved files only (tumor gene ranking, TPM matrices, cluster labels, InstaPrism fractions).
# It does NOT rerun InstaPrism and does NOT overwrite any existing file.
#
# The gene list is the SAME one used for the tumors (ranking file written by the tumor script),
# so the PDO heatmap rows are identical to the tumor heatmap rows.
# Each gene is z-scored using the mean and SD of the 31 TUMORS and that scale is applied to the PDOs,
# so a PDO value reads as "SDs above/below the tumor average" (z-scoring across 4 PDOs alone would be unstable).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths and settings ----------------
RNA_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR <- file.path(RNA_DIR, "PDO")
INSTA_T <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/InstaPrism_PBCA_reference/Stranded")
INSTA_P <- file.path(PDO_DIR, "InstaPrism_PBCA_reference")
OUT_DIR <- file.path(INSTA_P, "Downstream", "Cancer_epithelial_genes")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Gene ranking written by the tumor script (contains gene, rank, cluster_highest, heatmap_row)
RANK_FILE <- file.path(INSTA_T, "Downstream_stats", "Cancer_epithelial_genes",
                       "pILC_InstaPrism_PBCA_cancer_epithelial_gene_ranking.csv")
TPM_T     <- file.path(RNA_DIR, "4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt")
TPM_P     <- file.path(PDO_DIR, "Max_TPM/tpm_matrix_PDO.txt")
CLUSTER_F <- file.path(RNA_DIR, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
FRAC_T    <- file.path(INSTA_T, "pILC_InstaPrism_PBCA_fraction_with_summary.csv")   # optional (annotation strip)
FRAC_P    <- file.path(INSTA_P, "PDO_InstaPrism_PBCA_fraction_with_summary.csv")    # optional (annotation strip)

N_TOP          <- 100    # used only if the ranking file has no heatmap_row column
CLUSTER_LEVELS <- c("1", "2", "3")
group_cols     <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "PDO" = "#333333")
Z_CLIP         <- 2.5

# ---------------- 2. Helpers ----------------
need <- function(path, what) if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)

make_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- sub("([0-9])T$", "\\1", x)
  gsub("_0+([0-9])", "_\\1", x)
}

read_tpm <- function(path) {
  d <- fread(path, data.table = FALSE, check.names = FALSE)
  m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; storage.mode(m) <- "numeric"
  colnames(m) <- make_key(colnames(m))
  if (anyDuplicated(colnames(m))) stop("Duplicated sample keys in ", path)
  m
}

fmt_p <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3))))

for (f in c(RANK_FILE, TPM_T, TPM_P, CLUSTER_F)) need(f, "Input file")

# ---------------- 3. Gene list from the tumor analysis ----------------
rank_tab <- fread(RANK_FILE, data.table = FALSE)
if (!all(c("gene", "rank") %in% colnames(rank_tab))) stop("Ranking file needs the columns 'gene' and 'rank'.")
if ("heatmap_row" %in% colnames(rank_tab) && any(!is.na(rank_tab$heatmap_row))) {
  top <- rank_tab[!is.na(rank_tab$heatmap_row), ]
  top <- top[order(top$heatmap_row), ]                      # same row order as the tumor heatmap
} else {
  top <- head(rank_tab[order(rank_tab$rank), ], N_TOP)
}
message("Genes from the tumor ranking: ", nrow(top))

# ---------------- 4. Expression and clusters ----------------
tpm_t <- read_tpm(TPM_T)
tpm_p <- read_tpm(TPM_P)
if (ncol(tpm_p) != 4) stop("Expected 4 PDO samples, found ", ncol(tpm_p))

cl_raw <- fread(CLUSTER_F, data.table = FALSE)
s_hit <- grep("^sample$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
s_col <- if (length(s_hit) > 0) s_hit[1] else colnames(cl_raw)[1]
c_hit <- grep("^cluster$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(cl_raw)[ncol(cl_raw)]
cl <- data.frame(key = make_key(cl_raw[[s_col]]),
                 num = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cl_raw[[c_col]])))))
if (anyNA(cl$num) || !all(cl$num %in% 1:3)) stop("Cluster values must be 1, 2 or 3 in ", CLUSTER_F)
if (anyDuplicated(cl$key)) stop("Duplicated sample keys in the cluster file.")

tum_keys <- intersect(colnames(tpm_t), cl$key)
message("Tumors matched between TPM and clusters: ", length(tum_keys), " of ", ncol(tpm_t))
if (length(tum_keys) < 25) stop("Too few matched tumors. Check the sample IDs.")
tpm_t <- tpm_t[, tum_keys, drop = FALSE]
tum_cl <- factor(as.character(cl$num[match(tum_keys, cl$key)]), levels = CLUSTER_LEVELS)

genes <- top$gene[top$gene %in% rownames(tpm_t) & top$gene %in% rownames(tpm_p)]
if (length(genes) < nrow(top))
  warning(nrow(top) - length(genes), " genes are missing in the tumor or PDO TPM matrix and are dropped: ",
          paste(setdiff(top$gene, genes), collapse = ", "))
top <- top[match(genes, top$gene), ]

# ---------------- 5. z-scores: PDOs on the tumor scale ----------------
e_t <- log2(tpm_t[genes, , drop = FALSE] + 1)
e_p <- log2(tpm_p[genes, , drop = FALSE] + 1)
mu  <- rowMeans(e_t)
sdv <- apply(e_t, 1, sd)
sdv[sdv == 0] <- NA                                           # a constant gene cannot be z-scored
Z_t <- (e_t - mu) / sdv; Z_t[!is.finite(Z_t)] <- 0
Z_p <- (e_p - mu) / sdv; Z_p[!is.finite(Z_p)] <- 0

score_t <- colMeans(Z_t)                                      # signature score = mean z of the genes
score_p <- colMeans(Z_p)

# Optional annotation: InstaPrism cancer epithelial fraction
cf_of <- function(path, keys) {
  if (!file.exists(path)) return(rep(NA_real_, length(keys)))
  fr <- fread(path, data.table = FALSE, check.names = FALSE)
  if (!"Cancer Epithelial Cells" %in% colnames(fr)) return(rep(NA_real_, length(keys)))
  fr[["Cancer Epithelial Cells"]][match(keys, make_key(fr$Sample_ID))]
}
cf_t <- cf_of(FRAC_T, colnames(Z_t))
cf_p <- cf_of(FRAC_P, colnames(Z_p))

# PDO signature score relative to the tumors (descriptive)
score_tab <- data.frame(
  Sample = colnames(Z_p), Signature_score = unname(score_p),
  InstaPrism_cancer_fraction = cf_p,
  tumor_median_score = median(score_t), tumor_min_score = min(score_t), tumor_max_score = max(score_t),
  PDO_percentile_in_tumors = sapply(score_p, function(s) 100 * mean(score_t <= s)),
  n_genes = length(genes), row.names = NULL)
write.csv(score_tab, file.path(OUT_DIR, "PDO_cancer_epithelial_signature_scores.csv"), row.names = FALSE)
print(score_tab)

# PDO z-score of every gene (for inspection)
z_p_df <- as.data.frame(Z_p, check.names = FALSE)
colnames(z_p_df) <- paste0("z_", colnames(Z_p))
gene_out <- data.frame(gene = genes, specificity_rank = top$rank,
                       cluster_highest = if ("cluster_highest" %in% colnames(top)) top$cluster_highest else NA,
                       tumor_mean_log2TPM = unname(mu), tumor_sd_log2TPM = unname(sdv),
                       z_p_df, check.names = FALSE, row.names = NULL)
write.csv(gene_out, file.path(OUT_DIR, "PDO_cancer_epithelial_genes_zscores.csv"), row.names = FALSE)

# ---------------- 6. Heatmap: tumors by cluster, then the PDOs ----------------
ord_t <- order(tum_cl, -score_t)                              # cluster, then signature score (high to left)
ord_p <- order(-score_p)
mat <- cbind(Z_t[, ord_t, drop = FALSE], Z_p[, ord_p, drop = FALSE])
col_group <- factor(c(as.character(tum_cl[ord_t]), rep("PDO", ncol(Z_p))), levels = c(CLUSTER_LEVELS, "PDO"))
cfrac     <- c(cf_t[ord_t], cf_p[ord_p])

n_by <- table(col_group)
col_titles <- c(paste0("Cluster ", CLUSTER_LEVELS, " (n = ", n_by[CLUSTER_LEVELS], ")"), paste0("PDO (n = ", n_by["PDO"], ")"))

rownames(mat) <- paste0(rownames(mat), " (#", top$rank, ")")

top_ann <- HeatmapAnnotation(
  Group = col_group,
  `Cancer fraction` = cfrac,
  col = list(Group = group_cols,
             `Cancer fraction` = colorRamp2(c(0, 0.5, 1), c("white", "#74C476", "#00441B"))),
  na_col = "grey85", simple_anno_size = unit(4, "mm"), gap = unit(1, "mm"),
  annotation_name_side = "right", annotation_name_gp = gpar(fontface = "bold", fontsize = 11),
  annotation_legend_param = list(Group = list(title = "Group", at = c(CLUSTER_LEVELS, "PDO"),
                                              labels = c(paste("Cluster", CLUSTER_LEVELS), "PDO")),
                                 `Cancer fraction` = list(title = "InstaPrism\ncancer fraction")))

# Row blocks by the cluster where each gene is highest in the tumors (same blocks as the tumor figure)
row_split <- NULL
if ("cluster_highest" %in% colnames(top) && !anyNA(top$cluster_highest)) {
  n_g <- table(factor(as.character(top$cluster_highest), levels = CLUSTER_LEVELS))
  row_split <- factor(as.character(top$cluster_highest), levels = CLUSTER_LEVELS,
                      labels = paste0("Highest in\nCluster ", CLUSTER_LEVELS, "\n", as.integer(n_g), " genes"))
}

N_G <- length(genes)
row_font <- if (N_G <= 60) 9 else if (N_G <= 120) 6 else 4.5

ht <- Heatmap(
  mat,
  name = "z-score vs tumors\n(log2(TPM+1))",
  col = colorRamp2(c(-Z_CLIP, 0, Z_CLIP), c("#2166AC", "white", "#B2182B")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_group, column_title = col_titles,
  column_title_gp = gpar(fontface = "bold", fontsize = 12), column_gap = unit(2, "mm"),
  row_split = row_split, row_title = NULL, row_gap = unit(2, "mm"),
  row_names_side = "right", row_names_gp = gpar(fontsize = row_font),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = top_ann, border = FALSE, use_raster = FALSE)

fig_title <- paste0("Top ", N_G, " cancer-epithelial genes (PBCA reference): pILC tumors and PDOs")
fig_cap <- paste0("Gene list, row order and (#n) = specificity rank are the same as in the tumor heatmap. ",
                  "z-score per gene uses the mean and SD of the ", ncol(Z_t), " tumors; PDOs are placed on that scale.\n",
                  "Columns: tumors by RNA cluster, then PDOs, each ordered by signature score (mean z of the genes). ",
                  "Descriptive only (n = ", ncol(Z_p), " PDOs).")

draw_fig <- function() {
  draw(ht, column_title = fig_title, column_title_gp = gpar(fontface = "bold", fontsize = 15),
       heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE,
       padding = unit(c(22, 4, 4, 4), "mm"))
  grid.text(fig_cap, x = unit(4, "mm"), y = unit(10, "mm"), just = c("left", "center"), gp = gpar(fontsize = 9))
}
fig_name <- file.path(OUT_DIR, "PDO_cancer_epithelial_top_genes_heatmap")
fig_w <- 14
fig_h <- if (N_G <= 60) 0.17 * N_G + 4.5 else 0.1 * N_G + 4.5
png(paste0(fig_name, ".png"), width = fig_w, height = fig_h, units = "in", res = 300); draw_fig(); invisible(dev.off())
pdf(paste0(fig_name, ".pdf"), width = fig_w, height = fig_h);                          draw_fig(); invisible(dev.off())

message("Done. Files written to: ", OUT_DIR)
