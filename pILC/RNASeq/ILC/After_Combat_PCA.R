# ============================================================
# PCA of ILC RNA-seq before and after ComBat batch correction
# Batch = source (Bionano vs Chelsea)
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)        # ComBat
  library(ggplot2)
  library(ggrepel)
  library(patchwork)  # side-by-side plot (install.packages("patchwork") if missing)
})
set.seed(42)

# ---------- 0. Settings ----------
top_n     <- 2000   # number of most variable genes used for PCA
mean_only <- TRUE   # TRUE = adjust mean only (safer with 3 samples in one batch)

# ---------- 1. Paths and sample sheet ----------
base_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_ILC"
out_dir  <- file.path(base_dir, "Stranded_No")

bionano_ids <- c("TP10-1655", "TP12-813", "TP14-79")
chelsea_ids <- c("TP17_M378", "TP17_M383", "TP17_M384", "TP17_M510",
                 "TP17_M559", "TP17_M605", "TP17_M622", "TP17_M802",
                 "TP17_M882", "TP17_M899", "TP17_M901",
                 "TP18_M95",  "TP18_M121", "TP18_M202", "TP18_M251",
                 "TP18_M278", "TP18_M302", "TP18_M329", "TP18_M355",
                 "TP18_M372", "TP18_M403")

meta <- rbind(
  data.frame(sample_id = bionano_ids, source = "Bionano",
             file = file.path(out_dir, paste0(bionano_ids, ".txt"))),
  data.frame(sample_id = chelsea_ids, source = "Chelsea",
             file = file.path(base_dir, paste0(chelsea_ids, ".txt")))
)
rownames(meta) <- meta$sample_id

# ---------- 2. Read count files ----------
# ASSUMPTION: gene ID = first column, count = last column
read_counts <- function(f) {
  x <- read.delim(f, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[!grepl("^__", x[[1]]), ]
  x <- x[!x[[1]] %in% c("Geneid", "gene_id", "Gene"), ]
  v <- suppressWarnings(as.numeric(x[[ncol(x)]]))
  names(v) <- x[[1]]
  v[!is.na(v)]
}

count_list <- lapply(meta$file, read_counts)
names(count_list) <- meta$sample_id
common_genes <- Reduce(intersect, lapply(count_list, names))
counts <- round(sapply(count_list, function(v) v[common_genes]))
cat("Genes:", nrow(counts), " Samples:", ncol(counts), "\n")

# ---------- 3. Filter and VST ----------
counts <- counts[rowSums(counts >= 10) >= 3, ]
cat("Genes after filtering:", nrow(counts), "\n")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta[colnames(counts), ],
                              design    = ~ 1)
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else
  varianceStabilizingTransformation(dds, blind = TRUE)
mat <- assay(vsd)

# ---------- 4. ComBat batch correction ----------
batch <- meta[colnames(mat), "source"]

# ComBat fails if a gene has zero variance inside any batch -> remove those genes
zero_var <- apply(mat, 1, function(g) any(tapply(g, batch, var) == 0 |
                                          is.na(tapply(g, batch, var))))
mat <- mat[!zero_var, ]
cat("Genes removed (zero variance in a batch):", sum(zero_var), "\n")

mat_cb <- ComBat(dat = mat, batch = batch, mod = NULL,
                 par.prior = TRUE, mean.only = mean_only)

# ---------- 5. PCA helper (same settings before and after) ----------
run_pca <- function(m, label) {
  top_n_use <- min(top_n, nrow(m))
  m_top <- m[head(order(apply(m, 1, var), decreasing = TRUE), top_n_use), ]
  pca   <- prcomp(t(m_top), center = TRUE, scale. = FALSE)
  ve    <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  df <- data.frame(sample_id = rownames(pca$x),
                   PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
                   source = meta[rownames(pca$x), "source"],
                   stage = label)
  # R2 of PC1 explained by batch (how much PC1 is still a batch effect)
  r2 <- summary(lm(pca$x[, 1] ~ df$source))$r.squared
  list(df = df, ve = ve, n = top_n_use, r2 = r2)
}

res_before <- run_pca(mat,    "Before ComBat")
res_after  <- run_pca(mat_cb, "After ComBat")
cat(sprintf("PC1 variance explained by batch: before = %.1f%%, after = %.1f%%\n",
            100 * res_before$r2, 100 * res_after$r2))

# ---------- 6. Plot helper ----------
plot_pca <- function(res, title) {
  ggplot(res$df, aes(PC1, PC2, color = source)) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text_repel(aes(label = sample_id), size = 3, show.legend = FALSE,
                    max.overlaps = Inf) +
    scale_color_manual(values = c(Bionano = "#D55E00", Chelsea = "#0072B2")) +
    labs(x = sprintf("PC1 (%.1f%%)", res$ve[1]),
         y = sprintf("PC2 (%.1f%%)", res$ve[2]),
         title = title, color = "Source") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

p_before <- plot_pca(res_before, sprintf("Before ComBat (top %d genes)", res_before$n))
p_after  <- plot_pca(res_after,  sprintf("After ComBat (top %d genes)",  res_after$n))
p_both   <- p_before + p_after + plot_layout(guides = "collect")

# ---------- 7. Save (gene number in name) ----------
tag <- paste0("top", res_after$n)

ggsave(file.path(out_dir, paste0("PCA_ILC_after_ComBat_", tag, ".pdf")),
       p_after, width = 8, height = 6)
ggsave(file.path(out_dir, paste0("PCA_ILC_after_ComBat_", tag, ".png")),
       p_after, width = 8, height = 6, dpi = 300)

ggsave(file.path(out_dir, paste0("PCA_ILC_before_vs_after_ComBat_", tag, ".pdf")),
       p_both, width = 14, height = 6)
ggsave(file.path(out_dir, paste0("PCA_ILC_before_vs_after_ComBat_", tag, ".png")),
       p_both, width = 14, height = 6, dpi = 300)

write.csv(rbind(res_before$df, res_after$df),
          file.path(out_dir, paste0("PCA_ILC_coordinates_ComBat_", tag, ".csv")),
          row.names = FALSE)
write.csv(data.frame(stage = c("Before ComBat", "After ComBat"),
                     PC1_pct = c(res_before$ve[1], res_after$ve[1]),
                     PC2_pct = c(res_before$ve[2], res_after$ve[2]),
                     PC1_batch_R2 = c(res_before$r2, res_after$r2)),
          file.path(out_dir, paste0("PCA_ILC_ComBat_summary_", tag, ".csv")),
          row.names = FALSE)

# Corrected matrix, in case you need it for other analyses
saveRDS(mat_cb, file.path(out_dir, "VST_ComBat_corrected_matrix.rds"))