# ============================================================
# PCA of ILC RNA-seq: Bionano (3 samples) vs Chelsea (21 samples)
# Top N most variable genes (set top_n below)
# ============================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})
set.seed(42)

# ---------- 0. Settings ----------
top_n <- 2000   # number of most variable genes used for PCA

# ---------- 1. Paths and sample sheet ----------
base_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_ILC"
out_dir  <- file.path(base_dir, "Stranded_NO")

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
  x <- x[!grepl("^__", x[[1]]), ]                         # drop HTSeq summary lines
  x <- x[!x[[1]] %in% c("Geneid", "gene_id", "Gene"), ]   # drop header row if present
  v <- suppressWarnings(as.numeric(x[[ncol(x)]]))
  names(v) <- x[[1]]
  v[!is.na(v)]
}

count_list <- lapply(meta$file, read_counts)
names(count_list) <- meta$sample_id

common_genes <- Reduce(intersect, lapply(count_list, names))
counts <- round(sapply(count_list, function(v) v[common_genes]))
cat("Genes:", nrow(counts), " Samples:", ncol(counts), "\n")
print(colSums(counts))   # QC: total counts per sample

# ---------- 3. Filter low-expressed genes ----------
counts <- counts[rowSums(counts >= 10) >= 3, ]
cat("Genes after filtering:", nrow(counts), "\n")

# ---------- 4. VST normalization (unsupervised) ----------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta[colnames(counts), ],
                              design    = ~ 1)
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else
  varianceStabilizingTransformation(dds, blind = TRUE)
mat <- assay(vsd)

# ---------- 5. Top N most variable genes ----------
top_n_use <- min(top_n, nrow(mat))   # safety: cannot use more genes than available
top_idx   <- head(order(apply(mat, 1, var), decreasing = TRUE), top_n_use)
mat_top   <- mat[top_idx, ]
cat("Genes used for PCA:", nrow(mat_top), "\n")

# ---------- 6. PCA ----------
pca     <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  source = meta[rownames(pca$x), "source"]
)

# ---------- 7. Plot ----------
p <- ggplot(pca_df, aes(PC1, PC2, color = source)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = sample_id), size = 3, show.legend = FALSE,
                  max.overlaps = Inf) +
  scale_color_manual(values = c(Bionano = "#D55E00", Chelsea = "#0072B2")) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       title = sprintf("PCA of ILC RNA-seq (top %d variable genes, VST)", top_n_use),
       color = "Source") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# ---------- 8. Save figure and tables (gene number in name, no dates) ----------
tag <- paste0("top", top_n_use)
ggsave(file.path(out_dir, paste0("PCA_ILC_Bionano_vs_Chelsea_", tag, ".pdf")),
       p, width = 8, height = 6)
ggsave(file.path(out_dir, paste0("PCA_ILC_Bionano_vs_Chelsea_", tag, ".png")),
       p, width = 8, height = 6, dpi = 300)

write.csv(pca_df,
          file.path(out_dir, paste0("PCA_ILC_coordinates_", tag, ".csv")),
          row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(var_exp)),
                     variance_explained_pct = var_exp),
          file.path(out_dir, paste0("PCA_ILC_variance_explained_", tag, ".csv")),
          row.names = FALSE)