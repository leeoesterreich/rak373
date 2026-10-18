# ============================================================
# ESTIMATE (immunedeconv) on pILC RNA-seq, compared across NMF clusters
# Inputs : TPM matrix (primary) and VST matrix (sensitivity check)
# Outputs: scores CSV, stats CSV, 4-panel boxplot (PNG + PDF) per input
# Clusters: R consensus-NMF k=3 
# ============================================================

suppressPackageStartupMessages({
  library(immunedeconv)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)     # stat_pvalue_manual()
  library(rstatix)    # wilcox_test(), kruskal_test(), add_xy_position()
  library(patchwork)  # combine the 4 panels
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths and settings ----------------
base_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"

# New cluster membership (R consensus NMF, k = 3)
cluster_file <- file.path(base_dir, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")

# Corrected (htseq-count -s reverse) VST data, all expression-filtered genes (gene symbols)
# (not the top-5000 variable set: ESTIMATE needs its full signature gene list)
vst_file <- file.path(base_dir, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47",
                      "htseq_VST_gene_names_all_expression_filtered_genes.csv")
if (!file.exists(vst_file)) stop("VST file not found: ", vst_file)

# TPM matrix from the corrected stranded counts (made by TPM_pILC_stranded.R; run that script first)
tpm_file <- file.path(base_dir, "4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt")
if (!file.exists(tpm_file)) stop("TPM file not found. Run TPM_pILC_stranded.R first: ", tpm_file)

out_root <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/immunedeconv/Stranded"

K             <- 3
score_names   <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
cluster_cols  <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")  # Dark2, as in the example

# ---------------- 2. Helper functions ----------------

# Make sample names comparable across files:
#   TP03-117  and  TP03.117_S25  ->  TP03_117
#   TP09-9    and  TP09.9_S34    ->  TP09_9
make_key <- function(x) {
  x <- sub("_S[0-9]+$", "", x)          # remove the sequencing suffix (_S25, _S34, ...)
  gsub("[^A-Za-z0-9]+", "_", x)         # unify "-" and "." into "_"
}

# Read a genes x samples table; first column = gene symbol
read_expr <- function(path) {
  dt <- fread(path, data.table = FALSE)
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
  colnames(m) <- make_key(colnames(m))
  m
}

# Read the R consensus-NMF membership CSV (one row per sample).
# Sample column = "Sample" if present, otherwise the first column.
# Cluster column = "Cluster" (any case) if present, otherwise the last column.
read_clusters <- function(path, k) {
  g <- fread(path, data.table = FALSE)
  message("  Cluster file columns: ", paste(colnames(g), collapse = " | "))

  s_col <- if ("Sample" %in% colnames(g)) "Sample" else colnames(g)[1]
  c_hit <- grep("^cluster$", colnames(g), ignore.case = TRUE, value = TRUE)
  c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(g)[ncol(g)]
  message("  Using sample column '", s_col, "' and cluster column '", c_col, "'")

  # Keep only the digits so that "1" and "C1" both become 1
  cl_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(g[[c_col]]))))

  cl <- data.frame(Sample  = make_key(as.character(g[[s_col]])),
                   Cluster = factor(cl_num, levels = seq_len(k)))

  if (anyNA(cl$Cluster))
    stop("Some membership values are not in 1..", k, ". Please check ", path)
  if (anyDuplicated(cl$Sample))
    stop("Duplicated sample keys after ID harmonization.")
  cl
}

# Run ESTIMATE directly (same call as your earlier analysis); one row per sample
run_estimate <- function(expr) {
  res <- immunedeconv::deconvolute_estimate(expr)   # scores x samples
  res <- as.data.frame(res)
  message("  Score labels returned: ", paste(rownames(res), collapse = " | "))

  # Standardize labels (e.g. "Stromal score" -> "StromalScore") so downstream code is stable
  std <- c(stromalscore = "StromalScore", immunescore = "ImmuneScore",
           estimatescore = "ESTIMATEScore", tumorpurity = "TumorPurity")
  key <- tolower(gsub("[^A-Za-z]", "", rownames(res)))
  rownames(res) <- ifelse(key %in% names(std), std[key], rownames(res))

  missing <- setdiff(score_names, rownames(res))
  if (length(missing) > 0)
    stop("Scores not found: ", paste(missing, collapse = ", "),
         ". Labels returned are printed above.")

  res[score_names, , drop = FALSE] %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("Sample")
}

# Format a p-value for the figure: "= 0.012" or "< 0.001"
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3)))

# One boxplot panel + its statistics
#   - Kruskal-Wallis test across the clusters (panel subtitle)
#   - pairwise Wilcoxon rank-sum tests, BH-adjusted within this score (brackets)
plot_one_score <- function(long_df, score) {
  d <- long_df %>% filter(Score == score)

  # Pairwise Wilcoxon, BH-adjusted within this score; add bracket positions
  st <- d %>%
    wilcox_test(Value ~ Cluster, p.adjust.method = "BH") %>%
    add_xy_position(x = "Cluster")
  st$bracket_label <- paste0("adj. p ", fmt_p(st$p.adj))

  kw <- d %>% kruskal_test(Value ~ Cluster)

  # x-axis labels with the number of samples in each cluster
  n_per  <- table(d$Cluster)
  x_labs <- setNames(paste0("Cluster ", names(n_per), "\n(n = ", as.integer(n_per), ")"),
                     names(n_per))

  p <- ggplot(d, aes(x = Cluster, y = Value, fill = Cluster)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, colour = "grey20") +
    geom_point(position = position_jitter(width = 0.12, seed = SEED),
               size = 2, alpha = 0.7, colour = "black") +
    scale_fill_manual(values = cluster_cols) +
    scale_x_discrete(labels = x_labs) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    stat_pvalue_manual(st, label = "bracket_label", tip.length = 0.01,
                       size = 4.5, hide.ns = FALSE) +
    labs(title = score,
         subtitle = paste("Kruskal-Wallis p", fmt_p(kw$p)),
         x = NULL, y = "Score") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          plot.title    = element_text(face = "bold", size = 20),
          plot.subtitle = element_text(size = 14),
          axis.title.y  = element_text(face = "bold", size = 18))

  st_out <- st %>%
    transmute(Score = score, group1, group2, n1, n2,
              p_raw = p, p_adj_BH_within_score = p.adj,
              signif_label = p.adj.signif,
              kruskal_p = kw$p)
  list(plot = p, stats = st_out)
}

# Full analysis for one expression matrix
analyze <- function(expr, clusters, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  common <- intersect(colnames(expr), clusters$Sample)
  message("  Matched samples: ", length(common), " of ", ncol(expr),
          " (expression) / ", nrow(clusters), " (clusters)")
  if (length(common) < 6) stop("Too few matching samples. Check the sample IDs.")

  scores <- run_estimate(expr[, common]) %>%
    left_join(clusters, by = "Sample") %>%
    dplyr::select(Sample, Cluster, all_of(score_names))

  long_df <- scores %>%
    pivot_longer(all_of(score_names), names_to = "Score", values_to = "Value") %>%
    mutate(Score = factor(Score, levels = score_names))

  res   <- lapply(score_names, function(s) plot_one_score(long_df %>% mutate(Score = as.character(Score)), s))
  stats <- bind_rows(lapply(res, `[[`, "stats")) %>%
    mutate(FDR_BH_all_12_tests = p.adjust(p_raw, method = "BH"))

  # Save tables
  write.csv(scores, file.path(outdir, "ESTIMATE_scores.csv"), row.names = FALSE)
  write.csv(stats,  file.path(outdir, "ESTIMATE_scores_by_cluster_stats.csv"), row.names = FALSE)

  # Save the 4-panel figure (PNG + PDF)
  fig <- wrap_plots(lapply(res, `[[`, "plot"), nrow = 1) +
    plot_annotation(
      caption = paste0("Boxes show the median and interquartile range; points are individual samples.\n",
                       "Overall difference across clusters: Kruskal-Wallis test. ",
                       "Pairwise comparisons: Wilcoxon rank-sum test, Benjamini-Hochberg adjusted within each score."),
      theme = theme(plot.caption = element_text(size = 13, hjust = 0))
    )
  ggsave(file.path(outdir, "ESTIMATE_scores_by_cluster.png"), fig,
         width = 22, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "ESTIMATE_scores_by_cluster.pdf"), fig,
         width = 22, height = 7)

  scores
}

# ---------------- 3. Run ----------------
clusters <- read_clusters(cluster_file, K)
print(table(clusters$Cluster))   

message("TPM analysis")
scores_tpm <- analyze(read_expr(tpm_file), clusters, file.path(out_root, "TPM"))

message("VST analysis")
scores_vst <- analyze(read_expr(vst_file), clusters, file.path(out_root, "VST"))

# ---------------- 4. Sensitivity check: TPM vs VST ----------------
shared <- inner_join(scores_tpm, scores_vst, by = c("Sample", "Cluster"),
                     suffix = c("_TPM", "_VST"))

agreement <- data.frame(
  Score = score_names,
  Spearman_rho = sapply(score_names, function(s)
    cor(shared[[paste0(s, "_TPM")]], shared[[paste0(s, "_VST")]], method = "spearman"))
)
write.csv(agreement, file.path(out_root, "VST", "ESTIMATE_TPM_vs_VST_correlation.csv"),
          row.names = FALSE)
print(agreement)