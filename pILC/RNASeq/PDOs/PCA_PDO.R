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

# ------------------------------ 4b. Genes driving PC1 and PC2 -----------------
suppressPackageStartupMessages(library(pheatmap))
N_LOAD <- 15                                    # genes shown per direction per PC

load_df <- data.frame(Gene = rownames(pca$rotation),
                      PC1 = pca$rotation[, "PC1"], PC2 = pca$rotation[, "PC2"])
load_df$PC1_contribution_pct <- 100 * load_df$PC1^2    # rotation columns have unit length
load_df$PC2_contribution_pct <- 100 * load_df$PC2^2
load_df <- load_df[order(-abs(load_df$PC1)), ]
write.csv(load_df, file.path(OUT_DIR, "PCA_PDO_gene_loadings.csv"), row.names = FALSE)

# top genes per PC: N_LOAD positive and N_LOAD negative
pick <- function(pc) {
  v <- sort(setNames(load_df[[pc]], load_df$Gene))
  g <- c(head(v, N_LOAD), tail(v, N_LOAD))
  data.frame(PC = pc, Gene = names(g), Loading = as.numeric(g),
             Direction = ifelse(g > 0, "Higher in PC+ samples", "Higher in PC- samples"))
}
top_load <- rbind(pick("PC1"), pick("PC2"))
write.csv(top_load, file.path(OUT_DIR, "PCA_PDO_top_loading_genes.csv"), row.names = FALSE)

# bar plot of top loadings
top_load$Gene_ord <- paste(top_load$PC, top_load$Gene, sep = "_")
top_load$Gene_ord <- factor(top_load$Gene_ord, levels = top_load$Gene_ord[order(top_load$Loading)])
p_ld <- ggplot(top_load, aes(Loading, Gene_ord, fill = Direction)) +
  geom_col() +
  facet_wrap(~ PC, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^PC[12]_", "", x)) +
  scale_fill_manual(values = c("Higher in PC+ samples" = "#B2182B", "Higher in PC- samples" = "#2166AC")) +
  labs(title = "Genes with the largest PCA loadings (PDOs)", x = "Loading", y = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", plot.title.position = "plot", legend.title = element_blank())
ggsave(file.path(OUT_DIR, "PCA_PDO_top_loading_genes_barplot.png"), p_ld, width = 10, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "PCA_PDO_top_loading_genes_barplot.pdf"), p_ld, width = 10, height = 7)

# heatmap of the same genes (VST, centred per gene across the 4 PDOs)
hm_genes <- unique(top_load$Gene)
hm <- t(scale(t(vsd[hm_genes, , drop = FALSE]), center = TRUE, scale = FALSE))
png(file.path(OUT_DIR, "PCA_PDO_top_loading_genes_heatmap.png"), width = 5,
    height = 0.17 * length(hm_genes) + 2, units = "in", res = 300)
pheatmap(hm, cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA,
         main = "Top PC1/PC2 loading genes (VST, centred)")
dev.off()

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