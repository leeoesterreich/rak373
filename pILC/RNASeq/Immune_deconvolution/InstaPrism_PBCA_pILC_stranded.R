# ============================================================
# InstaPrism deconvolution of pILC bulk RNA-seq using the PBCA reference
# Input : TPM matrix from the corrected (htseq-count -s reverse) counts
#         (made by TPM_pILC_stranded.R; run that script first)
# Output: fractions, summary table, stacked barplot, heatmap, total immune plot
# ============================================================

suppressPackageStartupMessages({
  library(InstaPrism)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(patchwork)   # stack the barplot and the cluster strip
})

# ---------------- 1. Paths ----------------
BASE_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Immune_deconvolution_BC_Atlas"

# New output folder so the earlier (old-count) results are not overwritten
OUT_DIR <- file.path(BASE_DIR, "InstaPrism_PBCA_reference", "Stranded")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# TPM from the corrected stranded counts
BULK_TPM_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt"
if (!file.exists(BULK_TPM_FILE))
  stop("TPM file not found. Run TPM_pILC_stranded.R first: ", BULK_TPM_FILE)

# RNA cluster membership (R consensus NMF, 50 runs, k = 3)
CLUSTER_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv"
if (!file.exists(CLUSTER_FILE)) stop("Cluster file not found: ", CLUSTER_FILE)
K <- 3
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "Unassigned" = "grey70")  # Dark2, as in the ESTIMATE plots

PBCA_REF_FILE <- file.path(BASE_DIR, "PBCA_InstaPrism_reference.Rds")
if (!file.exists(PBCA_REF_FILE)) stop("Reference not found: ", PBCA_REF_FILE)

# ---------------- 2. Load InstaPrism reference ----------------
refPhi_obj <- readRDS(PBCA_REF_FILE)

cat("Reference object class:\n")
print(class(refPhi_obj))

cat("\nReference matrix dimensions:\n")
print(dim(refPhi_obj@phi.cs))

cat("\nReference cell types:\n")
print(names(refPhi_obj@map))

cat("\nNumber of cell states per cell type:\n")
print(sapply(refPhi_obj@map, length))

# ---------------- 3. Read pILC TPM matrix ----------------
bulk_df <- fread(BULK_TPM_FILE, data.table = FALSE)

colnames(bulk_df)[1] <- "Gene"
bulk_df$Gene <- as.character(bulk_df$Gene)

# Collapse duplicated genes by mean TPM, if any
bulk_df <- bulk_df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

bulk_expr <- as.data.frame(bulk_df)
rownames(bulk_expr) <- bulk_expr$Gene
bulk_expr$Gene <- NULL

bulk_expr <- as.matrix(bulk_expr)
bulk_expr[is.na(bulk_expr)] <- 0
bulk_expr[bulk_expr < 0] <- 0

cat("\nBulk TPM dimensions:\n")
print(dim(bulk_expr))
print(bulk_expr[1:5, 1:min(5, ncol(bulk_expr))])

# ---------------- 4. Match genes between bulk and reference ----------------
ref_genes <- rownames(refPhi_obj@phi.cs)
bulk_genes <- rownames(bulk_expr)

common_genes <- intersect(ref_genes, bulk_genes)

cat("Reference genes:", length(ref_genes), "\n")
cat("Bulk genes:", length(bulk_genes), "\n")
cat("Common genes:", length(common_genes), "\n")

if (length(common_genes) < 10000) {
  warning("Low gene overlap. Please check gene symbols.")
}

bulk_expr_use <- bulk_expr[common_genes, , drop = FALSE]

refPhi_use <- refPhi_obj
phi_use <- refPhi_obj@phi.cs[common_genes, , drop = FALSE]

# Re-normalize reference columns after gene filtering
phi_colsum <- colSums(phi_use)
phi_use <- sweep(phi_use, 2, phi_colsum, FUN = "/")
phi_use[!is.finite(phi_use)] <- 1e-8

refPhi_use@phi.cs <- phi_use

cat("\nFinal bulk matrix:\n")
print(dim(bulk_expr_use))

cat("\nFinal reference phi matrix:\n")
print(dim(refPhi_use@phi.cs))

write.csv(
  data.frame(
    Metric = c("Reference genes", "Bulk genes", "Common genes", "Bulk samples",
               "Reference cell states", "Reference cell types"),
    N = c(
      length(ref_genes),
      length(bulk_genes),
      length(common_genes),
      ncol(bulk_expr_use),
      ncol(refPhi_use@phi.cs),
      length(refPhi_use@map)
    )
  ),
  file.path(OUT_DIR, "PBCA_InstaPrism_input_summary.csv"),
  row.names = FALSE
)

# ---------------- 5. Run InstaPrism ----------------
deconv_res <- InstaPrism(
  bulk_Expr = bulk_expr_use,
  refPhi_cs = refPhi_use
)

saveRDS(
  deconv_res,
  file.path(OUT_DIR, "pILC_InstaPrism_PBCA_deconv_result.rds")
)

# ---------------- 6. Extract cell-type fractions ----------------
estimated_frac <- t(deconv_res@Post.ini.ct@theta)

estimated_frac <- as.data.frame(estimated_frac)
estimated_frac$Sample_ID <- rownames(estimated_frac)
estimated_frac <- estimated_frac %>% relocate(Sample_ID)

write.csv(
  estimated_frac,
  file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cell_fractions.csv"),
  row.names = FALSE
)

print(estimated_frac[1:min(5, nrow(estimated_frac)), ])

# ---------------- 7. Immune / stromal / epithelial summary ----------------
immune_cols <- intersect(
  c("CD8+ T Cells", "CD4+ T Cells", "Regulatory T Cells", "NK Cells", "B Cells",
    "Plasma Cells", "Macrophages", "Monocytes", "Dendritic Cells", "MDSCs",
    "Mast Cells", "Neutrophils"),
  colnames(estimated_frac)
)

stromal_cols <- intersect(
  c("Fibroblasts", "Perivascular-like (PVL) Cells", "Endothelial Cells"),
  colnames(estimated_frac)
)

epithelial_cols <- intersect(
  c("Cancer Epithelial Cells", "Epithelial Cells", "Myoepithelial Cells"),
  colnames(estimated_frac)
)

summary_df <- estimated_frac %>%
  mutate(
    Immune_Fraction     = rowSums(across(all_of(immune_cols)), na.rm = TRUE),
    Stromal_Fraction    = rowSums(across(all_of(stromal_cols)), na.rm = TRUE),
    Epithelial_Fraction = rowSums(across(all_of(epithelial_cols)), na.rm = TRUE)
  )

write.csv(
  summary_df,
  file.path(OUT_DIR, "pILC_InstaPrism_PBCA_fraction_with_summary.csv"),
  row.names = FALSE
)

print(head(arrange(summary_df, desc(Immune_Fraction)), 10))

# ---------------- 8. Stacked barplot, samples grouped by RNA cluster ----------------
# Make sample names comparable between the TPM and the cluster file:
#   TP03-117 and TP03.117_S25 -> TP03_117
make_key <- function(x) {
  x <- sub("_S[0-9]+$", "", x)
  gsub("[^A-Za-z0-9]+", "_", x)
}

# Read the cluster membership CSV (one row per sample).
# Sample column = "Sample" if present, otherwise the first column.
# Cluster column = "Cluster" (any case) if present, otherwise the last column.
cl_raw <- fread(CLUSTER_FILE, data.table = FALSE)
message("Cluster file columns: ", paste(colnames(cl_raw), collapse = " | "))
s_col <- if ("Sample" %in% colnames(cl_raw)) "Sample" else colnames(cl_raw)[1]
c_hit <- grep("^cluster$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(cl_raw)[ncol(cl_raw)]
message("Using sample column '", s_col, "' and cluster column '", c_col, "'")

cluster_df <- data.frame(
  key     = make_key(as.character(cl_raw[[s_col]])),
  Cluster = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cl_raw[[c_col]]))))
)
if (anyNA(cluster_df$Cluster) || any(!cluster_df$Cluster %in% seq_len(K)))
  stop("Some cluster values are not in 1..", K, ". Please check ", CLUSTER_FILE)
if (anyDuplicated(cluster_df$key))
  stop("Duplicated sample keys in the cluster file after ID harmonization.")

sample_info <- data.frame(Sample_ID = estimated_frac$Sample_ID,
                          key       = make_key(estimated_frac$Sample_ID)) %>%
  left_join(cluster_df, by = "key") %>%
  mutate(Cluster = ifelse(is.na(Cluster), "Unassigned", as.character(Cluster)))

message("Samples per cluster:")
message("Samples per cluster:"); print(table(sample_info$Cluster))
if (any(sample_info$Cluster == "Unassigned"))
  warning("Samples with no cluster: ",
          paste(sample_info$Sample_ID[sample_info$Cluster == "Unassigned"], collapse = ", "))

# Order: cluster 1, 2, 3 (then Unassigned), samples sorted by name within each cluster
sample_info <- sample_info %>%
  mutate(Cluster = factor(Cluster, levels = c(as.character(seq_len(K)), "Unassigned"))) %>%
  arrange(Cluster, Sample_ID)
sample_info$Sample_ID <- factor(sample_info$Sample_ID, levels = sample_info$Sample_ID)
sample_info <- droplevels(sample_info)

write.csv(sample_info %>% dplyr::select(Sample_ID, Cluster),
          file.path(OUT_DIR, "pILC_InstaPrism_PBCA_sample_cluster_order.csv"),
          row.names = FALSE)

frac_long <- estimated_frac %>%
  pivot_longer(
    cols = -Sample_ID,
    names_to = "Cell_Type",
    values_to = "Fraction"
  ) %>%
  left_join(sample_info %>% dplyr::select(Sample_ID, Cluster), by = "Sample_ID") %>%
  mutate(Sample_ID = factor(Sample_ID, levels = levels(sample_info$Sample_ID)))

# Top panel: stacked fractions, one block per cluster
p_stack_main <- ggplot(frac_long, aes(x = Sample_ID, y = Fraction, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text   = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing.x = unit(0.4, "lines")
  ) +
  labs(
    title = "pILC deconvolution using InstaPrism PBCA reference",
    y = "Estimated cell fraction",
    fill = "Cell type"
  )

# Bottom strip: cluster annotation (colour + number), with the sample labels under it
p_stack_strip <- ggplot(sample_info, aes(x = Sample_ID, y = 1, fill = Cluster)) +
  geom_tile(height = 1, colour = "white") +
  geom_text(aes(label = ifelse(Cluster == "Unassigned", "NA", as.character(Cluster))),
            colour = "white", fontface = "bold", size = 3.5) +
  facet_grid(. ~ Cluster, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cluster_cols, name = "RNA cluster",
                    labels = function(x) ifelse(x == "Unassigned", x, paste("Cluster", x))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    strip.text   = element_blank(),
    strip.background = element_blank(),
    panel.grid   = element_blank(),
    panel.spacing.x = unit(0.4, "lines")
  ) +
  labs(x = "Sample")

# guides = "collect" stacks both legends in one column on the right, so they cannot overlap
p_stack <- p_stack_main / p_stack_strip +
  plot_layout(heights = c(20, 1), guides = "collect")

ggsave(file.path(OUT_DIR, "pILC_InstaPrism_PBCA_stacked_barplot.pdf"),
       plot = p_stack, width = 13, height = 7)
ggsave(file.path(OUT_DIR, "pILC_InstaPrism_PBCA_stacked_barplot.png"),
       plot = p_stack, width = 13, height = 7, dpi = 300)

# ---------------- 9. Heatmap ----------------
frac_mat <- estimated_frac
rownames(frac_mat) <- frac_mat$Sample_ID
frac_mat$Sample_ID <- NULL

frac_mat <- as.matrix(frac_mat)
mode(frac_mat) <- "numeric"
frac_mat[!is.finite(frac_mat)] <- 0

frac_mat_t <- t(frac_mat)

# Drop cell types with zero variance (cannot be row-scaled)
row_sd <- apply(frac_mat_t, 1, sd, na.rm = TRUE)
frac_mat_t_plot <- frac_mat_t[row_sd > 0, , drop = FALSE]

draw_heatmap <- function() {
  pheatmap(
    frac_mat_t_plot,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = "InstaPrism PBCA deconvolution fractions"
  )
}

pdf(file.path(OUT_DIR, "pILC_InstaPrism_PBCA_fraction_heatmap.pdf"), width = 8, height = 8)
draw_heatmap()
invisible(dev.off())

png(file.path(OUT_DIR, "pILC_InstaPrism_PBCA_fraction_heatmap.png"),
    width = 8, height = 8, units = "in", res = 300)
draw_heatmap()
invisible(dev.off())

# ---------------- 10. Total immune fraction plot, bars coloured by RNA cluster ----------------
immune_df <- summary_df %>%
  left_join(sample_info %>% mutate(Sample_ID = as.character(Sample_ID)) %>%
              dplyr::select(Sample_ID, Cluster),
            by = "Sample_ID") %>%
  mutate(Cluster = ifelse(is.na(Cluster), "Unassigned", as.character(Cluster)))

p_immune <- ggplot(
  immune_df,
  aes(x = reorder(Sample_ID, Immune_Fraction), y = Immune_Fraction, fill = Cluster)
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = cluster_cols, name = "RNA cluster",
                    labels = function(x) ifelse(x == "Unassigned", x, paste("Cluster", x))) +
  theme_bw() +
  labs(
    title = "Estimated total immune fraction from InstaPrism PBCA reference",
    x = "Sample",
    y = "Immune fraction"
  )

ggsave(file.path(OUT_DIR, "pILC_InstaPrism_PBCA_total_immune_fraction.pdf"),
       plot = p_immune, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "pILC_InstaPrism_PBCA_total_immune_fraction.png"),
       plot = p_immune, width = 10, height = 5, dpi = 300)

# ---------------- 11. Session information ----------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_sessionInfo.txt"))

message("Done. Outputs in: ", OUT_DIR)