# pILC RNA-seq: sample-sample Pearson correlation heatmap (stranded VST), annotated by NMF cluster.
# Two versions: (1) VST values as they are (like the earlier figure), (2) after centering each gene
# across samples (shows sample differences more clearly). Each figure comes with its tables.

suppressPackageStartupMessages({
  library(pheatmap)        # install.packages("pheatmap")
  library(RColorBrewer)
  library(grid)
})
SEED <- 42
set.seed(SEED)

# ---------------------------------------------------------------- 1. Paths and settings
root    <- Sys.getenv("PILC_ROOT", "/ix1/alee/LO_LAB/Personal/Rahul/pILC")
vst_f   <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_VST_matrix_all_expression_filtered_genes.rds")
clus_f  <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
out_dir <- file.path(root, "RNASeq/Dimensionality_plot/Stranded")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
N_TOP_GENES <- NULL     # NULL = all expression-filtered genes; or e.g. 2000 for the most variable genes
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")   # Dark2, as in the other figures

norm_id <- function(x) {                       # TP09.9_S34 / TP09-9 -> TP09-09
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}

# ---------------------------------------------------------------- 2. Load VST and cluster labels
vst <- as.matrix(readRDS(vst_f))                 # expected genes x samples
if (nrow(vst) < ncol(vst)) vst <- t(vst)
if (is.null(colnames(vst))) stop("VST matrix has no sample names")
if (!is.null(N_TOP_GENES)) {
  top <- names(sort(apply(vst, 1, var), decreasing = TRUE))[seq_len(min(N_TOP_GENES, nrow(vst)))]
  vst <- vst[top, ]
}
clus <- read.csv(clus_f); clus$Key <- norm_id(clus$sample)
cl <- clus$cluster[match(norm_id(colnames(vst)), clus$Key)]
if (anyNA(cl)) stop("No cluster for: ", paste(colnames(vst)[is.na(cl)], collapse = ", "))
meta <- data.frame(Sample = colnames(vst), Cluster = factor(cl, levels = 1:3))
meta <- meta[order(meta$Cluster, meta$Sample), ]  # order: cluster, then sample name (as in the earlier figure)
vst  <- vst[, meta$Sample]
message("Using ", nrow(vst), " genes x ", ncol(vst), " samples")

# ---------------------------------------------------------------- 3. Correlation, heatmap and summary
run_version <- function(mat, tag, title, center_zero) {
  cmat <- cor(mat, method = "pearson")
  off  <- cmat[upper.tri(cmat)]
  write.csv(cmat, file.path(out_dir, paste0("Sample_sample_Pearson_", tag, ".csv")))

  # Summary per sample: mean correlation with the others, within and between clusters
  summ <- do.call(rbind, lapply(seq_len(ncol(cmat)), function(i) {
    others <- setdiff(seq_len(ncol(cmat)), i)
    same   <- others[meta$Cluster[others] == meta$Cluster[i]]
    between <- others[meta$Cluster[others] != meta$Cluster[i]]
    data.frame(Sample = meta$Sample[i], Cluster = meta$Cluster[i],
               mean_cor_all_others = mean(cmat[i, others]),
               mean_cor_same_cluster = mean(cmat[i, same]),
               mean_cor_other_clusters = mean(cmat[i, between]))
  }))
  summ$z_mean_cor <- as.numeric(scale(summ$mean_cor_all_others))
  summ$low_correlation_flag <- summ$z_mean_cor < -2          # unusually dissimilar to the rest
  write.csv(summ, file.path(out_dir, paste0("Sample_correlation_summary_", tag, ".csv")), row.names = FALSE)

  # Cluster-by-cluster mean correlation (off-diagonal pairs only)
  cc <- outer(levels(meta$Cluster), levels(meta$Cluster), Vectorize(function(a, b) {
    m <- cmat[meta$Cluster == a, meta$Cluster == b, drop = FALSE]
    if (a == b) mean(m[upper.tri(m)]) else mean(m)
  }))
  dimnames(cc) <- list(paste0("cluster", levels(meta$Cluster)), paste0("cluster", levels(meta$Cluster)))
  write.csv(cc, file.path(out_dir, paste0("Cluster_mean_correlation_", tag, ".csv")))

  # Heatmap (diagonal shown at the top of the colour scale)
  plot_mat <- cmat
  if (center_zero) {
    lim <- max(abs(off)); brks <- seq(-lim, lim, length.out = 101)
  } else {
    lim <- 1; brks <- seq(min(off), 1, length.out = 101)
  }
  diag(plot_mat) <- lim
  annot <- data.frame(Cluster = meta$Cluster, row.names = meta$Sample)
  ph <- pheatmap(plot_mat, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), breaks = brks,
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 annotation_col = annot, annotation_row = annot,
                 annotation_colors = list(Cluster = cluster_cols),
                 main = title, fontsize = 12, border_color = NA, silent = TRUE)
  png(file.path(out_dir, paste0("Sample_sample_Pearson_", tag, ".png")), width = 10, height = 9, units = "in", res = 300)
  grid.newpage(); grid.draw(ph$gtable); dev.off()
  pdf(file.path(out_dir, paste0("Sample_sample_Pearson_", tag, ".pdf")), width = 10, height = 9)
  grid.newpage(); grid.draw(ph$gtable); dev.off()
  print(summ[summ$low_correlation_flag, ])
  invisible(cmat)
}

run_version(vst, "VST", "Sample-Sample Pearson Correlation", center_zero = FALSE)
vst_c <- t(scale(t(vst), center = TRUE, scale = FALSE))          # center each gene across samples
vst_c <- vst_c[apply(vst_c, 1, var) > 0, ]
run_version(vst_c, "VST_gene_centered", "Sample-Sample Pearson Correlation (gene-centered VST)", center_zero = TRUE)

message("Done. Output in ", out_dir)
