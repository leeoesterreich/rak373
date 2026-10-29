# ============================================================
# ESTIMATE (immunedeconv) on the 4 PDO RNA-seq samples
# Inputs : PDO TPM matrix (primary) and PDO VST matrix (sensitivity check)
#          Optional: tumor ESTIMATE scores (for the PDO vs tumor figure)
# Outputs: scores CSV, PDO bar figure, PDO-vs-tumor figure (PNG + PDF) per input
# Note   : n = 4 and no RNA clusters, so NO group statistics are computed.
#          Run TPM_PDO.R and PDO_counts_to_VST_GENCODEv47.R first.
# ============================================================

suppressPackageStartupMessages({
  library(immunedeconv)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths and settings ----------------
RNA_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR <- file.path(RNA_DIR, "PDO")

tpm_file <- file.path(PDO_DIR, "Max_TPM", "tpm_matrix_PDO.txt")
vst_file <- file.path(PDO_DIR, "VST_data_GENCODEv47", "htseq_VST_gene_names_all_expression_filtered_genes.csv")
for (f in c(tpm_file, vst_file)) if (!file.exists(f)) stop("Input not found: ", f)

# Tumor ESTIMATE scores (columns Sample, Cluster, StromalScore, ...), used only for the comparison figure
tumor_est_file <- file.path(RNA_DIR, "immunedeconv/Stranded/TPM/ESTIMATE_scores.csv")

out_root <- file.path(PDO_DIR, "immunedeconv")

score_names  <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")     # Dark2, as in the tumor figures
pdo_cols     <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")            # same PDO colours as the PCA

# ---------------- 2. Helper functions ----------------

# Read a genes x samples table; first column = gene symbol. Sample names are kept as in the file (BO-046 ...).
read_expr <- function(path) {
  dt <- fread(path, data.table = FALSE, check.names = FALSE)
  m  <- as.matrix(dt[, -1])
  rownames(m) <- dt[[1]]
  storage.mode(m) <- "numeric"

  # ESTIMATE needs HGNC symbols, not Ensembl IDs
  if (mean(grepl("^ENSG", rownames(m))) > 0.5)
    stop("Row names look like Ensembl IDs. Convert to gene symbols first: ", path)

  # Collapse duplicated symbols by maximum expression
  if (anyDuplicated(rownames(m))) {
    message("  Collapsing duplicated gene symbols (max): ", path)
    idx <- split(seq_len(nrow(m)), rownames(m))
    m   <- do.call(rbind, lapply(idx, function(i) apply(m[i, , drop = FALSE], 2, max)))
  }
  m
}

# Run ESTIMATE directly (same call as the tumor analysis); one row per sample
run_estimate <- function(expr) {
  res <- immunedeconv::deconvolute_estimate(expr)   # scores x samples
  res <- as.data.frame(res)
  message("  Score labels returned: ", paste(rownames(res), collapse = " | "))

  # Standardize labels (e.g. "Stromal score" -> "StromalScore")
  std <- c(stromalscore = "StromalScore", immunescore = "ImmuneScore",
           estimatescore = "ESTIMATEScore", tumorpurity = "TumorPurity")
  key <- tolower(gsub("[^A-Za-z]", "", rownames(res)))
  rownames(res) <- ifelse(key %in% names(std), std[key], rownames(res))

  missing <- setdiff(score_names, rownames(res))
  if (length(missing) > 0)
    stop("Scores not found: ", paste(missing, collapse = ", "), ". Labels returned are printed above.")

  res[score_names, , drop = FALSE] %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("Sample")
}

# Tumor scores for the comparison figure (NULL if the file is missing)
read_tumor_scores <- function() {
  if (!file.exists(tumor_est_file)) { message("Tumor ESTIMATE file not found; comparison figure skipped."); return(NULL) }
  t <- fread(tumor_est_file, data.table = FALSE)
  if (!all(c("Sample", "Cluster", score_names) %in% colnames(t))) {
    message("Tumor ESTIMATE file lacks the expected columns; comparison figure skipped."); return(NULL)
  }
  t$Cluster <- factor(as.character(t$Cluster), levels = c("1", "2", "3"))
  t
}

# Full analysis for one expression matrix
analyze <- function(expr, outdir, tumor) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("  Samples: ", paste(colnames(expr), collapse = ", "))
  if (ncol(expr) != 4) stop("Expected 4 PDO samples, found ", ncol(expr))

  scores <- run_estimate(expr)
  write.csv(scores, file.path(outdir, "ESTIMATE_scores_PDO.csv"), row.names = FALSE)

  col_map <- setNames(pdo_cols[seq_len(nrow(scores))], scores$Sample)
  long_p  <- scores %>%
    pivot_longer(all_of(score_names), names_to = "Score", values_to = "Value") %>%
    mutate(Score = factor(Score, levels = score_names))

  # (a) PDO-only figure: one bar per PDO, one panel per score
  p_bar <- ggplot(long_p, aes(x = Sample, y = Value, fill = Sample)) +
    geom_col(width = 0.7, colour = "grey20") +
    facet_wrap(~ Score, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = col_map) +
    labs(title = "ESTIMATE scores of the PDOs", x = NULL, y = "Score") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 18))
  ggsave(file.path(outdir, "ESTIMATE_scores_PDO.png"), p_bar, width = 14, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "ESTIMATE_scores_PDO.pdf"), p_bar, width = 14, height = 5)

  # (b) PDO vs tumors: tumor boxplots by RNA cluster, PDOs as labelled points (descriptive only)
  if (!is.null(tumor)) {
    lv <- c("1", "2", "3", "PDO")
    n_t <- table(tumor$Cluster)
    x_labs <- c(setNames(paste0("Cluster ", names(n_t), "\n(n = ", as.integer(n_t), ")"), names(n_t)),
                PDO = paste0("PDO\n(n = ", nrow(scores), ")"))
    long_t <- tumor %>%
      dplyr::select(Sample, Cluster, all_of(score_names)) %>%
      pivot_longer(all_of(score_names), names_to = "Score", values_to = "Value") %>%
      mutate(Score = factor(Score, levels = score_names), Group = as.character(Cluster))
    long_pp <- long_p %>% mutate(Group = "PDO")

    p_cmp <- ggplot() +
      geom_boxplot(data = long_t, aes(x = Group, y = Value, fill = Group),
                   width = 0.6, outlier.shape = NA, colour = "grey20") +
      geom_point(data = long_t, aes(x = Group, y = Value),
                 position = position_jitter(width = 0.12, height = 0, seed = SEED),
                 size = 1.4, alpha = 0.6, colour = "grey30") +
      geom_point(data = long_pp, aes(x = Group, y = Value, colour = Sample),
                 position = position_jitter(width = 0.1, height = 0, seed = SEED), size = 3.5) +
      facet_wrap(~ Score, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = cluster_cols, guide = "none") +
      scale_colour_manual(values = col_map, name = "PDO") +
      scale_x_discrete(limits = lv, labels = x_labs) +
      labs(title = "ESTIMATE scores: PDOs versus the 31 tumors",
           x = NULL, y = "Score",
           caption = "Boxes: tumors by RNA cluster (median and interquartile range). Coloured points: individual PDOs. Descriptive only (n = 4 PDOs).") +
      theme_classic(base_size = 13) +
      theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 18),
            plot.caption = element_text(size = 10, hjust = 0))
    ggsave(file.path(outdir, "ESTIMATE_scores_PDO_vs_tumors.png"), p_cmp, width = 18, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(outdir, "ESTIMATE_scores_PDO_vs_tumors.pdf"), p_cmp, width = 18, height = 6)
  }
  scores
}

# ---------------- 3. Run ----------------
tumor <- read_tumor_scores()

message("TPM analysis")
scores_tpm <- analyze(read_expr(tpm_file), file.path(out_root, "TPM"), tumor)

message("VST analysis")
scores_vst <- analyze(read_expr(vst_file), file.path(out_root, "VST"), tumor)

# ---------------- 4. Sensitivity check: TPM vs VST ----------------
# Spearman rho across the 4 PDOs only describes agreement of the ordering (no p-value with n = 4).
shared <- inner_join(scores_tpm, scores_vst, by = "Sample", suffix = c("_TPM", "_VST"))
agreement <- data.frame(
  Score = score_names,
  n_samples = nrow(shared),
  Spearman_rho = sapply(score_names, function(s)
    cor(shared[[paste0(s, "_TPM")]], shared[[paste0(s, "_VST")]], method = "spearman"))
)
write.csv(agreement, file.path(out_root, "VST", "ESTIMATE_TPM_vs_VST_correlation_PDO.csv"), row.names = FALSE)
print(agreement)
message("Done. Outputs in: ", out_root)
