# =====================================================================
# PDO: boxplot of selected TF activities, tumors by RNA cluster with the 4 PDOs overlaid
# Reads tables saved by PDO_decoupler_TF_activity.R (no need to rerun it)
# The FDR in each panel header is the exploratory PDO-vs-tumors test (4 vs 31).
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
plot_name <- "TF_activity_selected_TFs_boxplot_PDO_vs_tumors"      # output name (PNG + PDF)
plot_w    <- 10                                                    # width, inches
plot_h    <- 7                                                     # height, inches
# ---------------------------------------------------------------------

root    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
out_dir <- file.path(root, "PDO", "Decoupler")
cluster_file <- file.path(root, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")

# ---- 1. Load saved tables -------------------------------------------
tf_p <- read_csv(file.path(out_dir, "TF_activity_per_sample_PDO.csv"),
                 show_col_types = FALSE, name_repair = "minimal")
tf_t <- read_csv(file.path(out_dir, "TF_activity_per_sample_tumors_common_genes.csv"),
                 show_col_types = FALSE, name_repair = "minimal")
stats <- read_csv(file.path(out_dir, "TF_PDO_vs_tumors_stats.csv"), show_col_types = FALSE) %>%
  distinct(TF, p, FDR)

missing_tfs <- setdiff(tfs, tf_p$TF)
if (length(missing_tfs) > 0)
  warning("TFs not scored (not in network or too few targets): ", paste(missing_tfs, collapse = ", "))
tfs <- tfs[tfs %in% tf_p$TF]
if (length(tfs) == 0) stop("None of the requested TFs were scored.")

# ---- 2. Cluster labels (same harmonization as the main script) ------
norm_id <- function(x) gsub("-", ".", x, fixed = TRUE)
clu <- read_csv(cluster_file, show_col_types = FALSE)
cl_raw <- as.character(clu$cluster)
clu$cluster_lab <- ifelse(grepl("^[0-9]+$", cl_raw), paste0("Cluster", cl_raw), cl_raw)
clu$sample <- norm_id(clu$sample)
cl_levels  <- sort(unique(clu$cluster_lab))
cl_map     <- setNames(clu$cluster_lab, clu$sample)

pal     <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")
cl_cols <- setNames(pal[seq_along(cl_levels)], cl_levels)
grp_levels <- c(cl_levels, "PDO")

# ---- 3. Long tables with FDR in the panel header --------------------
fdr_txt <- function(x) ifelse(x < 0.001, "< 0.001", sprintf("%.3f", x))

lab_df <- stats %>%
  filter(TF %in% tfs) %>%
  mutate(panel = paste0(TF, "\nPDO vs tumors FDR = ", fdr_txt(FDR)))
panel_levels <- lab_df$panel[match(tfs, lab_df$TF)]

tum_df <- tf_t %>%
  filter(TF %in% tfs) %>%
  pivot_longer(-TF, names_to = "sample", values_to = "score") %>%
  mutate(Group = cl_map[norm_id(sample)]) %>%
  left_join(select(lab_df, TF, panel), by = "TF") %>%
  mutate(panel = factor(panel, levels = panel_levels))
if (any(is.na(tum_df$Group))) stop("Some tumor samples have no cluster label; check sample IDs.")

pdo_df <- tf_p %>%
  filter(TF %in% tfs) %>%
  pivot_longer(-TF, names_to = "sample", values_to = "score") %>%
  mutate(Group = "PDO") %>%
  left_join(select(lab_df, TF, panel), by = "TF") %>%
  mutate(panel = factor(panel, levels = panel_levels))
pdo_cols <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[seq_along(unique(pdo_df$sample))],
                     unique(pdo_df$sample))

# ---- 4. Plot ---------------------------------------------------------
n_t <- table(factor(unique(tum_df[, c("sample", "Group")])$Group, levels = cl_levels))
x_labs <- c(setNames(paste0(sub("Cluster", "Cluster ", names(n_t)), "\n(n = ", as.integer(n_t), ")"), names(n_t)),
            PDO = paste0("PDO\n(n = ", length(pdo_cols), ")"))

p <- ggplot() +
  geom_boxplot(data = tum_df, aes(x = Group, y = score, fill = Group),
               outlier.shape = NA, alpha = 0.7) +
  geom_point(data = tum_df, aes(x = Group, y = score),
             position = position_jitter(width = 0.15, height = 0, seed = 42), size = 1.2, colour = "grey30") +
  geom_point(data = pdo_df, aes(x = Group, y = score, colour = sample),
             position = position_jitter(width = 0.1, height = 0, seed = 42), size = 3) +
  facet_wrap(~panel, scales = "free_y", ncol = n_col) +
  scale_fill_manual(values = cl_cols, guide = "none") +
  scale_colour_manual(values = pdo_cols, name = "PDO") +
  scale_x_discrete(limits = grp_levels, labels = x_labs[grp_levels]) +
  labs(y = "TF activity (ulm score)", x = NULL,
       caption = paste0("Boxes: tumors by RNA cluster; grey points: tumors; coloured points: PDOs. ",
                        "FDR: PDOs vs all tumors (Wilcoxon, BH over all scored TFs), exploratory.")) +
  theme_bw(base_size = 11) +
  theme(plot.caption = element_text(size = 9, hjust = 0))

ggsave(file.path(out_dir, paste0(plot_name, ".pdf")), p, width = plot_w, height = plot_h)
ggsave(file.path(out_dir, paste0(plot_name, ".png")), p, width = plot_w, height = plot_h, dpi = 300)

# Matching stats table for the plotted TFs
write_csv(lab_df %>% select(TF, p, FDR),
          file.path(out_dir, paste0(plot_name, "_stats.csv")))

cat("Saved:", file.path(out_dir, plot_name), "(.png/.pdf/_stats.csv)\n")
