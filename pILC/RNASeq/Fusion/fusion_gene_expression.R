#!/usr/bin/env Rscript
# =============================================================================
# Expression of recurrent fusion genes: carrier samples vs rest (pILC, log2TPM)
# Primary analysis : robust z-score + percentile of each carrier within "rest"
# Secondary        : exploratory Wilcoxon (genes with >= 2 carrier samples)
#                    and a pooled test of carrier z-scores against 0
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)  # reproducible jitter / any random step

# ---- 1. Paths ---------------------------------------------------------------
tpm_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC_log2TPM.txt"
out_dir  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/Fusion_Gene_Exp/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Carrier table (from your slide; unique samples per gene) ------------
carriers <- list(
  SLC25A13 = c("TP14-490"),
  ANO1     = c("TP13-1150", "TP13-820"),
  FCHSD2   = c("TP13-820",  "TP13-1150"),
  UVRAG    = c("TP13-820",  "TP13-1150"),
  RAB6A    = c("TP13-1150", "TP05-1063"),
  NEDD9    = c("TP05-668"),
  MYEOV    = c("TP04-104",  "TP13-820"),
  CARMIL1  = c("TP13-820",  "TP05-668"),
  MFSD11   = c("TP11-1516"),
  PITPNC1  = c("TP11-1516")
)

# ---- 3. Read the matrix (genes in rows, samples in columns) -----------------
mat <- fread(tpm_file, data.table = FALSE)
rownames(mat) <- make.unique(as.character(mat[[1]]))  # first column = gene symbol
mat <- mat[, -1, drop = FALSE]

# Make sample names comparable ("TP13_820" / "TP13.820" -> "TP13-820")
colnames(mat) <- gsub("[._]", "-", colnames(mat))

# Check that every carrier sample exists in the matrix
all_carriers <- unique(unlist(carriers))
missing <- setdiff(all_carriers, colnames(mat))
if (length(missing) > 0) {
  stop("Carrier samples not found in matrix columns: ", paste(missing, collapse = ", "),
       "\nExample column names: ", paste(head(colnames(mat), 5), collapse = ", "))
}

# The file is already log2(TPM [+1]); do NOT log-transform again.
cat("Matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
cat("Value range:", signif(range(as.matrix(mat), na.rm = TRUE), 3), "\n")

# ---- 4. Per-gene carrier vs rest statistics ---------------------------------
robust_z <- function(x, ref) {
  # Median/MAD of the reference are resistant to outliers in the rest group
  (x - median(ref)) / (1.4826 * mad(ref, constant = 1))
}

res_list <- list()
plot_df  <- list()

for (g in names(carriers)) {
  if (!g %in% rownames(mat)) { warning(g, " not in matrix; skipped"); next }

  expr    <- unlist(mat[g, ])
  is_carr <- names(expr) %in% carriers[[g]]
  rest    <- expr[!is_carr]
  carr    <- expr[is_carr]

  res_list[[g]] <- data.table(
    gene          = g,
    sample        = names(carr),
    log2TPM       = as.numeric(carr),
    rest_median   = median(rest),
    rest_n        = length(rest),
    robust_z      = robust_z(as.numeric(carr), rest),
    percentile    = vapply(carr, function(v) mean(rest <= v), numeric(1)),  # rank in rest
    n_carriers    = length(carr)
  )

  plot_df[[g]] <- data.table(
    gene = g, sample = names(expr), log2TPM = as.numeric(expr),
    group = ifelse(is_carr, "Fusion carrier", "Rest")
  )
}

res  <- rbindlist(res_list)
pdat <- rbindlist(plot_df)

# ---- 5. Exploratory tests ---------------------------------------------------
# (a) Per gene Wilcoxon, only when >= 2 carriers (very low power: exploratory!)
wil <- rbindlist(lapply(names(carriers), function(g) {
  d <- pdat[gene == g]
  nc <- sum(d$group == "Fusion carrier")
  if (nc < 2) return(data.table(gene = g, n_carriers = nc, wilcox_p = NA_real_))
  data.table(gene = g, n_carriers = nc,
             wilcox_p = wilcox.test(log2TPM ~ group, data = d, exact = FALSE)$p.value)
}))
wil[, wilcox_FDR := p.adjust(wilcox_p, method = "BH")]  # one family: the 10 genes

# (b) Pooled: are carrier z-scores shifted from 0 across all gene-sample pairs?
pooled_p <- wilcox.test(res$robust_z, mu = 0, exact = FALSE)$p.value
cat("Pooled one-sample Wilcoxon on carrier robust z (exploratory): p =",
    signif(pooled_p, 3), "\n")

# ---- 6. Save the statistics table -------------------------------------------
stats_out <- merge(res, wil[, .(gene, wilcox_p, wilcox_FDR)], by = "gene")
fwrite(stats_out, file.path(out_dir, "Fusion_gene_expression_stats.csv"))

# ---- 7. Figure: rest as boxplot, carriers highlighted -----------------------
pdat[, gene := factor(gene, levels = names(carriers))]

p <- ggplot(pdat[group == "Rest"], aes(x = gene, y = log2TPM)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey40") +
  geom_jitter(width = 0.15, size = 0.8, colour = "grey60") +
  geom_point(data = pdat[group == "Fusion carrier"],
             aes(colour = sample), size = 3) +
  labs(x = NULL, y = "log2 TPM", colour = "Carrier sample",
       title = "Expression of recurrent fusion genes",
       subtitle = "Grey = rest of cohort; coloured points = fusion-carrier samples") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "Fusion_gene_expression_boxplot.pdf"), p, width = 9, height = 5)
ggsave(file.path(out_dir, "Fusion_gene_expression_boxplot.png"), p, width = 9, height = 5, dpi = 300)

# ---- 8. Figure: carrier robust z-scores -------------------------------------
p2 <- ggplot(res, aes(x = gene, y = robust_z, colour = sample)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  labs(x = NULL, y = "Robust z-score vs rest", colour = "Carrier sample") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "Fusion_gene_carrier_zscore.pdf"), p2, width = 8, height = 4.5)
ggsave(file.path(out_dir, "Fusion_gene_carrier_zscore.png"), p2, width = 8, height = 4.5, dpi = 300)

cat("Done. Outputs written to:", out_dir, "\n")
