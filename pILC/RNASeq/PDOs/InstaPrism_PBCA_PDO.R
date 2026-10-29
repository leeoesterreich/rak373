# ============================================================
# InstaPrism deconvolution of the 4 PDO bulk RNA-seq samples using the PBCA reference
# Input : PDO TPM matrix (made by TPM_PDO.R; run that script first)
# Output: fractions, summary table, stacked barplot, heatmap, total immune plot
# Note  : no RNA clusters for PDOs, so the cluster strip of the tumor script is removed.
# ============================================================

suppressPackageStartupMessages({
  library(InstaPrism)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths ----------------
RNA_DIR  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR  <- file.path(RNA_DIR, "PDO")
OUT_DIR  <- file.path(PDO_DIR, "InstaPrism_PBCA_reference")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BULK_TPM_FILE <- file.path(PDO_DIR, "Max_TPM", "tpm_matrix_PDO.txt")
if (!file.exists(BULK_TPM_FILE)) stop("TPM file not found. Run TPM_PDO.R first: ", BULK_TPM_FILE)

# Same PBCA reference as the tumor analysis
PBCA_REF_FILE <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/PBCA_InstaPrism_reference.Rds")
if (!file.exists(PBCA_REF_FILE)) stop("Reference not found: ", PBCA_REF_FILE)

# ---------------- 2. Load InstaPrism reference ----------------
refPhi_obj <- readRDS(PBCA_REF_FILE)

cat("Reference object class:\n");                 print(class(refPhi_obj))
cat("\nReference matrix dimensions:\n");          print(dim(refPhi_obj@phi.cs))
cat("\nReference cell types:\n");                 print(names(refPhi_obj@map))
cat("\nNumber of cell states per cell type:\n");  print(sapply(refPhi_obj@map, length))

# ---------------- 3. Read PDO TPM matrix ----------------
bulk_df <- fread(BULK_TPM_FILE, data.table = FALSE, check.names = FALSE)
colnames(bulk_df)[1] <- "Gene"
bulk_df$Gene <- as.character(bulk_df$Gene)
if (ncol(bulk_df) - 1 != 4) stop("Expected 4 PDO samples, found ", ncol(bulk_df) - 1)

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

cat("\nBulk TPM dimensions:\n"); print(dim(bulk_expr))
print(bulk_expr[1:5, , drop = FALSE])

# ---------------- 4. Match genes between bulk and reference ----------------
ref_genes  <- rownames(refPhi_obj@phi.cs)
bulk_genes <- rownames(bulk_expr)

common_genes <- intersect(ref_genes, bulk_genes)

cat("Reference genes:", length(ref_genes), "\n")
cat("Bulk genes:", length(bulk_genes), "\n")
cat("Common genes:", length(common_genes), "\n")

if (length(common_genes) < 10000) warning("Low gene overlap. Please check gene symbols.")

bulk_expr_use <- bulk_expr[common_genes, , drop = FALSE]

refPhi_use <- refPhi_obj
phi_use <- refPhi_obj@phi.cs[common_genes, , drop = FALSE]

# Re-normalize reference columns after gene filtering (same as the tumor script)
phi_colsum <- colSums(phi_use)
phi_use <- sweep(phi_use, 2, phi_colsum, FUN = "/")
phi_use[!is.finite(phi_use)] <- 1e-8

refPhi_use@phi.cs <- phi_use

cat("\nFinal bulk matrix:\n");          print(dim(bulk_expr_use))
cat("\nFinal reference phi matrix:\n"); print(dim(refPhi_use@phi.cs))

write.csv(
  data.frame(
    Metric = c("Reference genes", "Bulk genes", "Common genes", "Bulk samples",
               "Reference cell states", "Reference cell types"),
    N = c(length(ref_genes), length(bulk_genes), length(common_genes), ncol(bulk_expr_use),
          ncol(refPhi_use@phi.cs), length(refPhi_use@map))
  ),
  file.path(OUT_DIR, "PBCA_InstaPrism_input_summary_PDO.csv"),
  row.names = FALSE
)

# ---------------- 5. Run InstaPrism ----------------
deconv_res <- InstaPrism(
  bulk_Expr = bulk_expr_use,
  refPhi_cs = refPhi_use
)

saveRDS(deconv_res, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_deconv_result.rds"))

# ---------------- 6. Extract cell-type fractions ----------------
estimated_frac <- t(deconv_res@Post.ini.ct@theta)
estimated_frac <- as.data.frame(estimated_frac, check.names = FALSE)
estimated_frac$Sample_ID <- rownames(estimated_frac)
estimated_frac <- estimated_frac %>% relocate(Sample_ID)

write.csv(estimated_frac, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_cell_fractions.csv"), row.names = FALSE)
print(estimated_frac)

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

write.csv(summary_df, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_fraction_with_summary.csv"), row.names = FALSE)

# Quick look: organoids should be almost entirely cancer epithelial
print(summary_df[, c("Sample_ID", "Immune_Fraction", "Stromal_Fraction", "Epithelial_Fraction")])
if ("Cancer Epithelial Cells" %in% colnames(summary_df)) {
  low <- summary_df$Sample_ID[summary_df[["Cancer Epithelial Cells"]] < 0.8]
  if (length(low) > 0)
    message("Check: cancer epithelial fraction below 0.8 in ", paste(low, collapse = ", "),
            " (residual non-tumor cells, or model leakage).")
}

# ---------------- 8. Stacked barplot (one bar per PDO) ----------------
sample_levels <- sort(estimated_frac$Sample_ID)

frac_long <- estimated_frac %>%
  pivot_longer(cols = -Sample_ID, names_to = "Cell_Type", values_to = "Fraction") %>%
  mutate(Sample_ID = factor(Sample_ID, levels = sample_levels))

p_stack <- ggplot(frac_long, aes(x = Sample_ID, y = Fraction, fill = Cell_Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PDO deconvolution using InstaPrism PBCA reference",
       x = "Sample", y = "Estimated cell fraction", fill = "Cell type")

ggsave(file.path(OUT_DIR, "PDO_InstaPrism_PBCA_stacked_barplot.pdf"), p_stack, width = 8, height = 7)
ggsave(file.path(OUT_DIR, "PDO_InstaPrism_PBCA_stacked_barplot.png"), p_stack, width = 8, height = 7, dpi = 300)

# ---------------- 9. Heatmap of fractions ----------------
# Not row-scaled: with only 4 samples a per-row z-score would exaggerate tiny differences.
frac_mat <- estimated_frac
rownames(frac_mat) <- frac_mat$Sample_ID
frac_mat$Sample_ID <- NULL
frac_mat <- as.matrix(frac_mat)
mode(frac_mat) <- "numeric"
frac_mat[!is.finite(frac_mat)] <- 0

frac_mat_t <- t(frac_mat)
frac_mat_t_plot <- frac_mat_t[apply(frac_mat_t, 1, max) > 0, , drop = FALSE]   # drop all-zero cell types

draw_heatmap <- function() {
  pheatmap(
    frac_mat_t_plot,
    scale = "none",
    cluster_cols = FALSE,
    cluster_rows = nrow(frac_mat_t_plot) > 1,
    display_numbers = TRUE, number_format = "%.3f", fontsize_number = 7,
    main = "InstaPrism PBCA deconvolution fractions (PDOs)"
  )
}

pdf(file.path(OUT_DIR, "PDO_InstaPrism_PBCA_fraction_heatmap.pdf"), width = 7, height = 8)
draw_heatmap()
invisible(dev.off())

png(file.path(OUT_DIR, "PDO_InstaPrism_PBCA_fraction_heatmap.png"),
    width = 7, height = 8, units = "in", res = 300)
draw_heatmap()
invisible(dev.off())

# ---------------- 10. Total immune fraction plot ----------------
p_immune <- ggplot(summary_df, aes(x = reorder(Sample_ID, Immune_Fraction), y = Immune_Fraction)) +
  geom_col(fill = "grey40", width = 0.6) +
  coord_flip() +
  theme_bw() +
  labs(title = "Estimated total immune fraction (InstaPrism PBCA) in PDOs",
       x = "Sample", y = "Immune fraction")

ggsave(file.path(OUT_DIR, "PDO_InstaPrism_PBCA_total_immune_fraction.pdf"), p_immune, width = 8, height = 4)
ggsave(file.path(OUT_DIR, "PDO_InstaPrism_PBCA_total_immune_fraction.png"), p_immune, width = 8, height = 4, dpi = 300)

# ---------------- 11. Session information ----------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_sessionInfo_PDO.txt"))

message("Done. Outputs in: ", OUT_DIR)
