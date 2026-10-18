#!/usr/bin/env Rscript
# =============================================================================
# Waterfall plots (one separate figure per gene): CDH1, ERBB2, ESR1
# Bars = log2TPM per sample, ranked high -> low, filled with a colour gradient
# Input : log2TPM matrix (genes x samples)
# Output: per gene PDF + PNG, and one table with all plotted values
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)

# ---- 1. Settings ------------------------------------------------------------
tpm_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC_log2TPM.txt"
out_dir  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Specific_gene_exp/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genes <- c("CDH1", "ERBB2", "ESR1")

# Optional footnote under each plot, e.g. "21 ILC samples are from Chelsea dataset"
caption_text <- NULL

# ---- 2. Read the matrix (genes in rows, samples in columns) -----------------
mat <- fread(tpm_file, data.table = FALSE)
rownames(mat) <- make.unique(as.character(mat[[1]]))  # first column = gene symbol
mat <- mat[, -1, drop = FALSE]

cat("Matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
if (ncol(mat) != 31) {
  warning("Expected 31 samples but found ", ncol(mat), ". Please check the matrix.")
}

missing_genes <- setdiff(genes, rownames(mat))
if (length(missing_genes) > 0) {
  stop("Gene(s) not found in matrix: ", paste(missing_genes, collapse = ", "))
}

# ---- 3. One waterfall plot per gene -----------------------------------------
all_values <- list()

for (g in genes) {

  x <- unlist(mat[g, ])
  d <- data.table(sample = names(x), log2TPM = as.numeric(x))

  if (anyNA(d$log2TPM)) {
    warning(g, ": ", sum(is.na(d$log2TPM)),
            " sample(s) have missing values (shown without a bar).")
  }

  # Rank high -> low; missing values go to the end. All samples stay on the axis.
  d <- d[order(-log2TPM)]
  d[, rank := seq_len(.N)]
  d[, sample := factor(sample, levels = sample)]
  all_values[[g]] <- cbind(gene = g, d)

  p <- ggplot(d, aes(x = sample, y = log2TPM, fill = log2TPM)) +
    geom_col(width = 0.8) +
    scale_fill_gradient2(low = "blue", mid = "grey85", high = "#E8836F",
                         midpoint = median(d$log2TPM, na.rm = TRUE),
                         limits = c(0, max(d$log2TPM, na.rm = TRUE)),
                         name = "log2TPM") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste(g, "Expression Waterfall Plot"),
         x = "Sample ID",
         y = paste0(g, " Expression (log2TPM)"),
         caption = caption_text) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

  ggsave(file.path(out_dir, paste0("Waterfall_", g, ".pdf")), p, width = 10, height = 6)
  ggsave(file.path(out_dir, paste0("Waterfall_", g, ".png")), p, width = 10, height = 6, dpi = 300)
}

# ---- 4. Save the values used in the figures ---------------------------------
fwrite(rbindlist(all_values),
       file.path(out_dir, "Waterfall_CDH1_ERBB2_ESR1_values.csv"))

cat("Done. Outputs written to:", out_dir, "\n")
