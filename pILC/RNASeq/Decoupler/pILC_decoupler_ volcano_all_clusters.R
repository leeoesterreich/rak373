# =====================================================================
# pILC: TF activity volcano plots, every cluster vs the rest (3 panels)
# All TFs, colored only by direction and significance (no annotation)
# Needs TF_cluster_stats.csv in the decoupler folder
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

SEED <- 42
set.seed(SEED)

# ---- EDIT HERE ------------------------------------------------------
fdr_cut          <- 0.05   # FDR line and significance cutoff
n_label          <- 8      # TFs labeled per direction per panel (largest |difference|)
highlight_tfs    <- c("AR", "ESR1", "FOXA1")  # always labeled if present
label_complexes  <- FALSE  # FALSE: complexes (names with "_") are plotted but not labeled
# ---------------------------------------------------------------------

root <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
# All results are read from and saved to this one folder
out_dir <- file.path(root, "Decoupler")
if (!dir.exists(out_dir)) stop("Folder not found: ", out_dir)
cat("Using folder:", out_dir, "\n")

# ---- 1. Load TF statistics ------------------------------------------
stats <- read_csv(file.path(out_dir, "TF_cluster_stats.csv"),
                  show_col_types = FALSE)
stopifnot(all(c("TF", "Cluster", "diff", "p", "FDR") %in% colnames(stats)))

lv <- c("Not significant", "Lower in cluster", "Higher in cluster")
df <- stats %>%
  mutate(
    neg_log_fdr = -log10(FDR),
    status = case_when(FDR < fdr_cut & diff > 0 ~ "Higher in cluster",
                       FDR < fdr_cut & diff < 0 ~ "Lower in cluster",
                       TRUE ~ "Not significant"),
    status = factor(status, levels = lv)
  )

# ---- 2. Panel titles with counts ------------------------------------
counts <- df %>%
  group_by(Cluster) %>%
  summarise(n_TFs = n(),
            n_higher = sum(status == "Higher in cluster"),
            n_lower  = sum(status == "Lower in cluster"),
            .groups = "drop") %>%
  mutate(panel = paste0(sub("Cluster", "Cluster ", Cluster), " vs rest\n",
                        n_higher, " higher, ", n_lower, " lower (FDR < ",
                        fdr_cut, ")"))
print(counts)

df <- df %>%
  left_join(select(counts, Cluster, panel), by = "Cluster") %>%
  mutate(panel = factor(panel, levels = counts$panel[order(counts$Cluster)])) %>%
  arrange(status)   # grey drawn first, significant points on top

# ---- 3. Labels: largest effects per direction + chosen TFs ----------
lab_pool <- if (label_complexes) df else df %>% filter(!grepl("_", TF))
lab_sig <- lab_pool %>%
  filter(status != "Not significant") %>%
  group_by(Cluster, status) %>%
  arrange(desc(abs(diff)), .by_group = TRUE) %>%
  slice_head(n = n_label) %>%
  ungroup()
lab_hi <- lab_pool %>% filter(TF %in% highlight_tfs)
lab <- bind_rows(lab_sig, lab_hi) %>% distinct(Cluster, TF, .keep_all = TRUE)

# Matching tables for the figure
write_csv(counts %>% select(-panel),
          file.path(out_dir, "TF_activity_volcano_all_clusters_counts.csv"))
write_csv(lab %>% select(Cluster, TF, diff, p, FDR, status),
          file.path(out_dir, "TF_activity_volcano_all_clusters_labeled_TFs.csv"))

# All significant TFs behind the plot (every cluster, with direction)
sig_tbl <- df %>%
  filter(status != "Not significant") %>%
  transmute(Cluster, TF, direction = as.character(status), diff, p, FDR) %>%
  arrange(Cluster, direction, desc(abs(diff)))
write_csv(sig_tbl,
          file.path(out_dir, "TF_activity_volcano_all_clusters_significant_TFs.csv"))
cat("Significant TFs saved:", nrow(sig_tbl), "rows\n")

# ---- 4. Volcano plot, three panels ----------------------------------
cols <- c("Not significant" = "grey70",
          "Lower in cluster" = "#2166ac",
          "Higher in cluster" = "#b2182b")

p <- ggplot(df, aes(x = diff, y = neg_log_fdr, colour = status)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~panel, ncol = 3) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = "Mean TF activity difference (cluster vs rest)",
       y = "-log10 FDR (one-vs-rest Wilcoxon, BH)",
       caption = paste0("Dashed line: FDR = ", fdr_cut, ". Labels: up to ", n_label,
                        " TFs per direction with the largest |difference| among ",
                        "significant TFs, plus ", paste(highlight_tfs, collapse = ", "),
                        ". Same axes in all panels.")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  p + ggrepel::geom_text_repel(data = lab, aes(label = TF), size = 2.8,
                               max.overlaps = Inf, show.legend = FALSE,
                               colour = "grey20", segment.size = 0.2)
} else {
  p + geom_text(data = lab, aes(label = TF), size = 2.8, vjust = -0.6,
                check_overlap = TRUE, show.legend = FALSE, colour = "grey20")
}

ggsave(file.path(out_dir, "TF_activity_volcano_all_clusters.pdf"), p,
       width = 15, height = 6)
ggsave(file.path(out_dir, "TF_activity_volcano_all_clusters.png"), p,
       width = 15, height = 6, dpi = 300)

cat("Done. Results in:", out_dir, "\n")