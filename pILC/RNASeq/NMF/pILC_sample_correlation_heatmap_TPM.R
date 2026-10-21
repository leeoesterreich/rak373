# pILC RNA-seq: sample-sample Pearson correlation heatmaps from the TPM matrices, annotated by NMF cluster.
# Three versions (each figure comes with its tables, no dates in names):
#   TPM                      linear TPM  (Max_TPM/tpm_matrix_pILC.txt)
#   log2TPM                  log2 TPM    (Max_TPM/tpm_matrix_pILC_log2TPM.txt)
#   log2TPM_gene_centered    the log2 TPM matrix with every gene centered across samples (shows sample differences best)
# Genes: same genes as the VST matrix used in Sample_sample_Pearson_VST (so the comparison is like for like).

suppressPackageStartupMessages({
  library(pheatmap)        # install.packages("pheatmap")
  library(RColorBrewer)
  library(grid)
})
SEED <- 42
set.seed(SEED)

# ---------------------------------------------------------------- 1. Paths and settings
root     <- Sys.getenv("PILC_ROOT", "/ix1/alee/LO_LAB/Personal/Rahul/pILC")
tpm_f    <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt")        # linear TPM
ltpm_f   <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC_log2TPM.txt") # log2 TPM
vst_f    <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_VST_matrix_all_expression_filtered_genes.rds")
clus_f   <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")  # 50-run clusters, GenePattern numbering (columns: sample, cluster)
out_dir  <- file.path(root, "RNASeq/Dimensionality_plot/Stranded")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(tpm_f, ltpm_f, vst_f, clus_f)) if (!file.exists(f)) stop("File not found: ", f)
USE_VST_GENES <- TRUE    # TRUE = use only the genes present in the VST matrix (same gene set as the VST heatmap)
MIN_TPM <- 1; MIN_SAMPLES <- 5   # used only if the VST gene set cannot be matched: TPM >= MIN_TPM in >= MIN_SAMPLES samples
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")   # Dark2, as in the other figures

norm_id <- function(x) {                       # TP09.9_S34 / TP09-9 -> TP09-09
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}

# ---------------------------------------------------------------- 2. Load the matrices and cluster labels
# Reads a genes x samples table (first column = gene name; .csv = comma, anything else = tab separated)
read_expr <- function(f) {
  raw <- if (grepl("\\.csv$", f, ignore.case = TRUE)) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE) else
    read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  gene <- as.character(raw[[1]])
  num  <- vapply(raw[-1], is.numeric, logical(1))             # keep numeric columns only (drops e.g. a description column)
  m <- as.matrix(raw[-1][, num, drop = FALSE]); mode(m) <- "numeric"
  rownames(m) <- gene; colnames(m) <- norm_id(colnames(m))
  if (anyNA(m)) stop("Missing values in ", f)
  if (anyDuplicated(gene)) {                                   # repeated gene names: keep the row with the highest mean
    m <- m[order(rowMeans(m), decreasing = TRUE), ]
    message(basename(f), ": ", sum(duplicated(rownames(m))), " duplicated gene rows removed")
    m <- m[!duplicated(rownames(m)), ]
  }
  message(basename(f), ": ", nrow(m), " genes x ", ncol(m), " samples")
  m
}
tpm  <- read_expr(tpm_f)
ltpm <- read_expr(ltpm_f)
if (min(tpm) < 0) stop("Negative values in the linear TPM matrix: is this really TPM?")
if (any(!is.finite(ltpm))) stop("Infinite values in the log2 TPM matrix (log2 of 0 without a pseudo-count?)")
if (!setequal(colnames(tpm), colnames(ltpm))) stop("The two matrices do not have the same samples")
genes_both <- intersect(rownames(tpm), rownames(ltpm))
tpm <- tpm[genes_both, ]; ltpm <- ltpm[genes_both, colnames(tpm)]

# Check how the log2 file was made (informational: which pseudo-count?)
chk <- max(abs(ltpm - log2(tpm + 1)))
message("Largest difference between the log2 TPM file and log2(TPM + 1) from the linear file: ", signif(chk, 3),
        if (chk < 1e-4) "  -> log2TPM file = log2(TPM + 1)" else "  -> log2TPM file was made differently (other pseudo-count?)")

clus <- read.csv(clus_f); clus$Key <- norm_id(clus$sample)
cl <- clus$cluster[match(colnames(tpm), clus$Key)]
if (anyNA(cl)) stop("No cluster for: ", paste(colnames(tpm)[is.na(cl)], collapse = ", "))
meta <- data.frame(Sample = colnames(tpm), Cluster = factor(cl, levels = 1:3))
meta <- meta[order(meta$Cluster, meta$Sample), ]              # order: cluster, then sample name (as in the VST figure)
tpm  <- tpm[, meta$Sample]; ltpm <- ltpm[, meta$Sample]
if (ncol(tpm) != 31) warning("Expected 31 samples, found ", ncol(tpm))

# Gene set: same genes as the VST heatmap, if the gene names can be matched
if (USE_VST_GENES) {
  v <- as.matrix(readRDS(vst_f)); if (nrow(v) < ncol(v)) v <- t(v)
  common <- intersect(rownames(tpm), rownames(v))
  message("Genes shared with the VST matrix: ", length(common), " of ", nrow(v), " VST genes")
  if (length(common) >= 1000) {
    keep <- common
  } else {
    warning("Fewer than 1000 shared genes (different gene IDs?). Using the TPM filter instead: TPM >= ", MIN_TPM,
            " in >= ", MIN_SAMPLES, " samples.")
    keep <- rownames(tpm)[rowSums(tpm >= MIN_TPM) >= MIN_SAMPLES]
  }
} else {
  keep <- rownames(tpm)[rowSums(tpm >= MIN_TPM) >= MIN_SAMPLES]
}
tpm <- tpm[keep, ]; ltpm <- ltpm[keep, ]
message("Using ", nrow(tpm), " genes x ", ncol(tpm), " samples")

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
  message(tag, ": correlation range between different samples = ", round(min(off), 3), " to ", round(max(off), 3))
  print(summ[summ$low_correlation_flag, ])
  invisible(cmat)
}

run_version(tpm, "TPM", "Sample-Sample Pearson Correlation (TPM)", center_zero = FALSE)
run_version(ltpm, "log2TPM", "Sample-Sample Pearson Correlation (log2 TPM)", center_zero = FALSE)
ltpm_c <- t(scale(t(ltpm), center = TRUE, scale = FALSE))        # center each gene across samples
ltpm_c <- ltpm_c[apply(ltpm_c, 1, var) > 0, ]
run_version(ltpm_c, "log2TPM_gene_centered", "Sample-Sample Pearson Correlation (gene-centered log2 TPM)", center_zero = TRUE)

message("Done. Output in ", out_dir)
