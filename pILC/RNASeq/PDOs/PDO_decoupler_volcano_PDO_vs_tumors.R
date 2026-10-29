# =====================================================================
# PDO: TF activity volcano, PDOs versus all pILC tumors (one panel)
# All TFs, colored only by direction and significance
# EXPLORATORY: 4 PDOs vs 31 tumors; differences also reflect culture conditions and the
# absence of stromal/immune cells in organoids, not only tumor-intrinsic biology.
# Needs TF_PDO_vs_tumors_stats.csv from PDO_decoupler_TF_activity.R
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
n_label          <- 8      # TFs labeled per direction (largest |difference|)
highlight_tfs    <- c("AR", "ESR1", "FOXA1")  # always labeled if present
label_complexes  <- FALSE  # FALSE: complexes (names with "_") are plotted but not labeled
# ---------------------------------------------------------------------

root    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
out_dir <- file.path(root, "PDO", "Decoupler")
if (!dir.exists(out_dir)) stop("Folder not found: ", out_dir)
cat("Using folder:", out_dir, "\n")

# ---- 1. Load TF statistics ------------------------------------------
stats <- read_csv(file.path(out_dir, "TF_PDO_vs_tumors_stats.csv"), show_col_types = FALSE)
stopifnot(all(c("TF", "diff", "p", "FDR") %in% colnames(stats)))

lv <- c("Not significant", "Lower in PDOs", "Higher in PDOs")
df <- stats %>%
  mutate(
    neg_log_fdr = -log10(FDR),
    status = case_when(FDR < fdr_cut & diff > 0 ~ "Higher in PDOs",
                       FDR < fdr_cut & diff < 0 ~ "Lower in PDOs",
                       TRUE ~ "Not significant"),
    status = factor(status, levels = lv)
  ) %>%
  arrange(status)   # grey drawn first, significant points on top

counts <- tibble(n_TFs = nrow(df),
                 n_higher = sum(df$status == "Higher in PDOs"),
                 n_lower  = sum(df$status == "Lower in PDOs"))
print(counts)
panel_title <- paste0("PDOs (n = 4) vs tumors (n = 31)\n",
                      counts$n_higher, " higher, ", counts$n_lower, " lower (FDR < ", fdr_cut, ")")

# ---- 2. Labels: largest effects per direction + chosen TFs ----------
lab_pool <- if (label_complexes) df else df %>% filter(!grepl("_", TF))
lab_sig <- lab_pool %>%
  filter(status != "Not significant") %>%
  group_by(status) %>%
  arrange(desc(abs(diff)), .by_group = TRUE) %>%
  slice_head(n = n_label) %>%
  ungroup()
lab_hi <- lab_pool %>% filter(TF %in% highlight_tfs)
lab <- bind_rows(lab_sig, lab_hi) %>% distinct(TF, .keep_all = TRUE)

# Matching tables for the figure
write_csv(counts, file.path(out_dir, "TF_activity_volcano_PDO_vs_tumors_counts.csv"))
write_csv(lab %>% select(TF, diff, p, FDR, status),
          file.path(out_dir, "TF_activity_volcano_PDO_vs_tumors_labeled_TFs.csv"))
sig_tbl <- df %>%
  filter(status != "Not significant") %>%
  transmute(TF, direction = as.character(status), diff, p, FDR) %>%
  arrange(direction, desc(abs(diff)))
write_csv(sig_tbl, file.path(out_dir, "TF_activity_volcano_PDO_vs_tumors_significant_TFs.csv"))
cat("Significant TFs saved:", nrow(sig_tbl), "rows\n")

# ---- 3. Volcano plot -------------------------------------------------
cols <- c("Not significant" = "grey70", "Lower in PDOs" = "#2166ac", "Higher in PDOs" = "#b2182b")

p <- ggplot(df, aes(x = diff, y = neg_log_fdr, colour = status)) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey40") +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = panel_title,
       x = "Mean TF activity difference (PDOs minus tumors)",
       y = "-log10 FDR (Wilcoxon, BH)",
       caption = paste0("Exploratory (4 vs 31); differences also reflect culture and cell composition. Dashed line: FDR = ", fdr_cut,
                        ". Labels: up to ", n_label, " TFs per direction with the largest |difference|, plus ",
                        paste(highlight_tfs, collapse = ", "), ".")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12))

p <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  p + ggrepel::geom_text_repel(data = lab, aes(label = TF), size = 2.8,
                               max.overlaps = Inf, show.legend = FALSE,
                               colour = "grey20", segment.size = 0.2)
} else {
  p + geom_text(data = lab, aes(label = TF), size = 2.8, vjust = -0.6,
                check_overlap = TRUE, show.legend = FALSE, colour = "grey20")
}

ggsave(file.path(out_dir, "TF_activity_volcano_PDO_vs_tumors.pdf"), p, width = 7, height = 6.5)
ggsave(file.path(out_dir, "TF_activity_volcano_PDO_vs_tumors.png"), p, width = 7, height = 6.5, dpi = 300)

cat("Done. Results in:", out_dir, "\n")
