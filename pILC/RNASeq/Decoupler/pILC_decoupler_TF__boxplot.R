# =====================================================================
# pILC: boxplot of selected TF activities by RNA cluster, with FDR on top
# Reads tables saved by pILC_decoupler_TF_activity.R (no need to rerun it)
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# ---- EDIT HERE ------------------------------------------------------
tfs       <- c("AR", "ESR1", "FOXA1", "HDAC1", "CIITA", "BARX1")  # TFs to plot
n_col     <- 3                                                     # panels per row
plot_name <- "TF_activity_selected_TFs_boxplot_FDR"                # output name (PNG + PDF)
plot_w    <- 10                                                    # width, inches
plot_h    <- 7                                                     # height, inches
# ---------------------------------------------------------------------

root    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
out_dir <- file.path(root, "Decoupler")
cluster_file <- file.path(
  root, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv"
)

# ---- 1. Load saved tables -------------------------------------------
tf_df <- read_csv(file.path(out_dir, "TF_activity_per_sample.csv"),
                  show_col_types = FALSE, name_repair = "minimal")
stats <- read_csv(file.path(out_dir, "TF_cluster_stats.csv"),
                  show_col_types = FALSE) %>%
  distinct(TF, KW_p, KW_FDR)   # global test across clusters, one row per TF

missing_tfs <- setdiff(tfs, tf_df$TF)
if (length(missing_tfs) > 0) {
  warning("TFs not scored (not in network or too few targets): ",
          paste(missing_tfs, collapse = ", "))
}
tfs <- tfs[tfs %in% tf_df$TF]
if (length(tfs) == 0) stop("None of the requested TFs were scored.")

# ---- 2. Cluster labels (same harmonization as the main script) ------
norm_id <- function(x) gsub("-", ".", x, fixed = TRUE)
clu <- read_csv(cluster_file, show_col_types = FALSE)
cl_raw <- as.character(clu$cluster)
clu$cluster_lab <- ifelse(grepl("^[0-9]+$", cl_raw),
                          paste0("Cluster", cl_raw), cl_raw)
clu$sample <- norm_id(clu$sample)
cl_levels  <- sort(unique(clu$cluster_lab))
cl_map     <- setNames(clu$cluster_lab, clu$sample)

pal     <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
cl_cols <- setNames(pal[seq_along(cl_levels)], cl_levels)

# ---- 3. Long table with FDR in the panel header ---------------------
fdr_txt <- function(x) ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))

lab_df <- stats %>%
  filter(TF %in% tfs) %>%
  mutate(panel = paste0(TF, "\nKruskal-Wallis FDR = ", fdr_txt(KW_FDR)))
panel_levels <- lab_df$panel[match(tfs, lab_df$TF)]

box_df <- tf_df %>%
  filter(TF %in% tfs) %>%
  pivot_longer(-TF, names_to = "sample", values_to = "score") %>%
  mutate(Cluster = factor(cl_map[norm_id(sample)], levels = cl_levels)) %>%
  left_join(select(lab_df, TF, panel), by = "TF") %>%
  mutate(panel = factor(panel, levels = panel_levels))

if (any(is.na(box_df$Cluster))) {
  stop("Some samples have no cluster label; check sample IDs.")
}

# ---- 4. Plot ---------------------------------------------------------
set.seed(42)  # fixed jitter positions
p <- ggplot(box_df, aes(x = Cluster, y = score, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2) +
  facet_wrap(~panel, scales = "free_y", ncol = n_col) +
  scale_fill_manual(values = cl_cols, guide = "none") +
  labs(y = "TF activity (ulm score)", x = NULL,
       caption = "FDR: Kruskal-Wallis across clusters, BH-adjusted over all scored TFs") +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, paste0(plot_name, ".pdf")), p,
       width = plot_w, height = plot_h)
ggsave(file.path(out_dir, paste0(plot_name, ".png")), p,
       width = plot_w, height = plot_h, dpi = 300)

# Matching stats table for the plotted TFs
write_csv(lab_df %>% select(TF, KW_p, KW_FDR),
          file.path(out_dir, paste0(plot_name, "_stats.csv")))

cat("Saved:", file.path(out_dir, plot_name), "(.png/.pdf/_stats.csv)\n")