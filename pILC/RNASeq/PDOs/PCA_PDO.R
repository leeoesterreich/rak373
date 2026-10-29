###############################################################################
# PCA of VST-normalised RNA-seq (4 PDO samples)
# Purpose: quick overview of how the 4 PDO samples relate to each other.
# Note: with n = 4 samples a PCA has only 3 non-trivial PCs, and no batch
#       variable is defined for the PDOs, so no batch statistics are computed.
#
# Input : RNASeq/PDO/VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv
#         (column "Gene" + one column per sample; RAW counts)
# Output: RNASeq/PDO/Dimensionality_plot/
#           PCA_PDO.png / .pdf
#           PCA_PDO_coordinates.csv   (PC1-PC3 per sample)
#           PCA_PDO_variance.csv      (variance explained per PC)
###############################################################################
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})
SEED <- 42; set.seed(SEED)

# ------------------------------ 1. Settings ----------------------------------
PDO_ROOT   <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/PDO"
COUNT_FILE <- file.path(PDO_ROOT, "VST_data_GENCODEv47", "htseq_raw_counts_with_gene_names.csv")
OUT_DIR    <- file.path(PDO_ROOT, "Dimensionality_plot")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PCA_NTOP    <- 2000   # PCA uses the PCA_NTOP most variable genes
MIN_COUNT   <- 10     # a gene must have >= MIN_COUNT counts ...
MIN_SAMPLES <- 2      # ... in at least MIN_SAMPLES samples (5 is not possible with 4 samples)

# ------------------------------ 2. Read raw counts ---------------------------
if (!file.exists(COUNT_FILE)) stop("Count file not found: ", COUNT_FILE)
raw <- read.csv(COUNT_FILE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("Gene" %in% colnames(raw))
genes <- make.unique(as.character(raw$Gene))
cnt <- as.matrix(raw[, setdiff(colnames(raw), "Gene"), drop = FALSE])
mode(cnt) <- "numeric"
rownames(cnt) <- genes
if (anyNA(cnt)) stop("Missing values found in the count matrix")
cnt <- round(cnt)
message("Count matrix: ", nrow(cnt), " genes x ", ncol(cnt), " samples")
if (ncol(cnt) != 4) stop("Expected 4 samples, found ", ncol(cnt))

# ------------------------------ 3. VST (blind) -------------------------------
keep <- rowSums(cnt >= MIN_COUNT) >= MIN_SAMPLES
message("Genes kept after low-count filter: ", sum(keep), " of ", nrow(cnt))
dds <- DESeqDataSetFromMatrix(countData = cnt[keep, ],
                              colData = data.frame(sample = colnames(cnt), row.names = colnames(cnt)),
                              design = ~ 1)
vsd <- assay(vst(dds, blind = TRUE))

# ------------------------------ 4. PCA on the most variable genes ------------
rv  <- apply(vsd, 1, var)
top <- names(sort(rv, decreasing = TRUE))[seq_len(min(PCA_NTOP, length(rv)))]
pca <- prcomp(t(vsd[top, ]), center = TRUE, scale. = FALSE)
pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
n_pc <- min(3, ncol(pca$x))                       # 4 samples -> at most 3 informative PCs

scores <- as.data.frame(pca$x[, seq_len(n_pc), drop = FALSE])
scores$sample <- rownames(scores)
write.csv(scores[, c("sample", paste0("PC", seq_len(n_pc)))],
          file.path(OUT_DIR, "PCA_PDO_coordinates.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(pve)), Variance_explained_pct = round(pve, 2)),
          file.path(OUT_DIR, "PCA_PDO_variance.csv"), row.names = FALSE)
print(round(pve, 2))

# ------------------------------ 5. Plot --------------------------------------
pdo_cols <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
p <- ggplot(scores, aes(PC1, PC2, color = sample)) +
  geom_point(size = 4.5, alpha = 0.9) +
  scale_color_manual(values = setNames(pdo_cols[seq_len(nrow(scores))], scores$sample)) +
  labs(title = sprintf("PCA of VST data, PDOs (top %d variable genes)", length(top)),
       x = sprintf("PC1 (%.1f%% variance)", pve[1]), y = sprintf("PC2 (%.1f%% variance)", pve[2])) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 16),
        axis.title = element_text(face = "bold"), legend.title = element_blank())
ggsave(file.path(OUT_DIR, "PCA_PDO.png"), p, width = 7, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, "PCA_PDO.pdf"), p, width = 7, height = 6)
message("Done. Output in ", OUT_DIR)
