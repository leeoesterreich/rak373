# ============================================================
# Downstream summary for the PDO InstaPrism (PBCA reference) results
#
# READS saved files only. It does NOT rerun InstaPrism and does NOT overwrite
# any existing file; everything new goes to PDO/InstaPrism_PBCA_reference/Downstream/.
#
# With n = 4 PDOs and no RNA clusters the tumor script's inferential parts
# (Kruskal-Wallis across clusters, Grade 2 vs 3, Spearman with clinical variables,
# plasma-cell sensitivity) cannot be done. This script is DESCRIPTIVE:
#   1. Key fractions + ESTIMATE scores per PDO, with the position of each PDO in the tumor distribution
#   2. Figure: the 9 key fractions, tumors by cluster with the PDOs overlaid
#   3. Cell-state fractions per PDO versus the tumor mean
# ============================================================

suppressPackageStartupMessages({
  library(InstaPrism)   # only so readRDS() can rebuild the saved S4 InstaPrism object
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths and settings ----------------
RNA_DIR    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR    <- file.path(RNA_DIR, "PDO")
INSTA_P    <- file.path(PDO_DIR, "InstaPrism_PBCA_reference")                       # PDO results
INSTA_T    <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/InstaPrism_PBCA_reference/Stranded")  # tumor results
OUT_DIR    <- file.path(INSTA_P, "Downstream")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# PDO inputs (required)
FRAC_P   <- file.path(INSTA_P, "PDO_InstaPrism_PBCA_fraction_with_summary.csv")
DECONV_P <- file.path(INSTA_P, "PDO_InstaPrism_PBCA_deconv_result.rds")
EST_P    <- file.path(PDO_DIR, "immunedeconv/TPM/ESTIMATE_scores_PDO.csv")          # optional

# Tumor inputs (optional: used for the comparison)
FRAC_T    <- file.path(INSTA_T, "pILC_InstaPrism_PBCA_fraction_with_summary.csv")
DECONV_T  <- file.path(INSTA_T, "pILC_InstaPrism_PBCA_deconv_result.rds")
CLUSTER_F <- file.path(RNA_DIR, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
EST_T     <- file.path(RNA_DIR, "immunedeconv/Stranded/TPM/ESTIMATE_scores.csv")
PBCA_REF_FILE <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/PBCA_InstaPrism_reference.Rds")  # optional: state -> cell type

CLUSTER_LEVELS <- c("1", "2", "3")
cluster_cols   <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")
pdo_cols       <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
TOP_STATES     <- 25

# Same 9 features as the tumor boxplots (name in the CSV = label on the plot)
PLOT_FEATURES <- c(
  "Fibroblasts"                   = "Fibroblasts (CAF)",
  "Perivascular-like (PVL) Cells" = "PVL Cells",
  "Endothelial Cells"             = "Endothelial Cells",
  "Stromal_Fraction"              = "Stromal fraction",
  "Cancer Epithelial Cells"       = "Cancer Epithelial Cells",
  "Epithelial_Fraction"           = "Epithelial fraction",
  "Macrophages"                   = "Macrophages (TAM)",
  "Plasma Cells"                  = "Plasma Cells",
  "Immune_Fraction"               = "Immune fraction"
)

COMPARTMENTS <- list(
  CAF              = c("Fibroblasts"),
  Stromal_vascular = c("Endothelial Cells", "Perivascular-like (PVL) Cells"),
  Epithelial       = c("Cancer Epithelial Cells", "Epithelial Cells", "Myoepithelial Cells"),
  TAM              = c("Macrophages"),
  Other_immune     = c("CD8+ T Cells", "CD4+ T Cells", "Regulatory T Cells", "NK Cells", "B Cells",
                       "Plasma Cells", "Monocytes", "Dendritic Cells", "MDSCs", "Mast Cells", "Neutrophils")
)
TUMOUR_STATE_PATTERN <- "TNBC|HER2|ER\\+|Luminal|_BC_|CID|Patient"   # patient-specific cancer-cell profiles

# ---------------- 2. Helpers ----------------
need <- function(path, what) if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)

# Same sample-key harmonisation as the other scripts: BO-046 -> BO_46, TP09-09T -> TP09_9
make_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- sub("([0-9])T$", "\\1", x)
  gsub("_0+([0-9])", "_\\1", x)
}

comp_of <- function(ct) {
  hit <- names(COMPARTMENTS)[vapply(COMPARTMENTS, function(v) ct %in% v, logical(1))]
  if (length(hit) > 0) hit[1] else NA_character_
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
}

# Cell-state matrix (samples x states) from a saved InstaPrism object, rows named by sample key
get_states <- function(rds, keys) {
  obj <- readRDS(rds)
  if (!"Post.ini.cs" %in% slotNames(obj)) stop("no 'Post.ini.cs' element in ", rds)
  cs    <- slot(obj, "Post.ini.cs")
  theta <- if (isS4(cs) && "theta" %in% slotNames(cs)) slot(cs, "theta") else cs
  theta <- as.matrix(theta)
  if (is.null(rownames(theta)) || is.null(colnames(theta))) stop("state matrix has no row/column names")
  hit_rows <- sum(make_key(rownames(theta)) %in% keys)
  hit_cols <- sum(make_key(colnames(theta)) %in% keys)
  m <- if (hit_cols >= hit_rows) t(theta) else theta       # samples must be rows
  rownames(m) <- make_key(rownames(m))
  if (!all(keys %in% rownames(m))) stop("not all samples are present in the state matrix of ", rds)
  m[keys, , drop = FALSE]
}

# ---------------- 3. Load the data ----------------
need(FRAC_P, "PDO fraction file"); need(DECONV_P, "PDO InstaPrism result")

frac_p <- fread(FRAC_P, data.table = FALSE, check.names = FALSE)
frac_p$key <- make_key(frac_p$Sample_ID)
if (nrow(frac_p) != 4) stop("Expected 4 PDO samples, found ", nrow(frac_p))
miss <- setdiff(names(PLOT_FEATURES), colnames(frac_p))
if (length(miss) > 0) stop("Expected columns not in the PDO fraction file: ", paste(miss, collapse = ", "))
pdo_col_map <- setNames(pdo_cols[seq_len(nrow(frac_p))], frac_p$Sample_ID)

# ESTIMATE scores of the PDOs (optional)
if (file.exists(EST_P)) {
  e <- fread(EST_P, data.table = FALSE)
  e$key <- make_key(e$Sample)
  e <- e[, c("key", "StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")]
  frac_p <- left_join(frac_p, e, by = "key")
} else message("PDO ESTIMATE file not found; ESTIMATE columns skipped.")

# Tumors (optional)
have_t <- all(file.exists(c(FRAC_T, CLUSTER_F)))
if (have_t) {
  frac_t <- fread(FRAC_T, data.table = FALSE, check.names = FALSE)
  frac_t$key <- make_key(frac_t$Sample_ID)
  cl_raw <- fread(CLUSTER_F, data.table = FALSE)
  s_hit <- grep("^sample$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
  s_col <- if (length(s_hit) > 0) s_hit[1] else colnames(cl_raw)[1]
  c_hit <- grep("^cluster$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
  c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(cl_raw)[ncol(cl_raw)]
  cl <- data.frame(key = make_key(cl_raw[[s_col]]),
                   num = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cl_raw[[c_col]])))))
  if (anyNA(cl$num) || !all(cl$num %in% 1:3)) stop("Cluster values must be 1, 2 or 3 in ", CLUSTER_F)
  frac_t$Cluster <- factor(as.character(cl$num[match(frac_t$key, cl$key)]), levels = CLUSTER_LEVELS)
  if (anyNA(frac_t$Cluster)) stop("No cluster found for: ", paste(frac_t$Sample_ID[is.na(frac_t$Cluster)], collapse = ", "))
  message("Tumors loaded: ", nrow(frac_t)); print(table(frac_t$Cluster))
} else message("Tumor fraction or cluster file not found; tumor comparison skipped.")

# ---------------- 4. Key fractions per PDO and position in the tumor distribution ----------------
feat <- names(PLOT_FEATURES)
pdo_wide <- frac_p[, c("Sample_ID", feat, intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"), colnames(frac_p)))]
write.csv(pdo_wide, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_key_fractions_per_PDO.csv"), row.names = FALSE)
print(pdo_wide)

summ <- do.call(rbind, lapply(feat, function(f) {
  do.call(rbind, lapply(seq_len(nrow(frac_p)), function(i) {
    v  <- frac_p[[f]][i]
    tv <- if (have_t) frac_t[[f]] else NA_real_
    data.frame(feature = f, Sample = frac_p$Sample_ID[i], PDO_value = v,
               tumor_median = median(tv), tumor_min = min(tv), tumor_max = max(tv),
               PDO_percentile_in_tumors = if (have_t) 100 * mean(tv <= v) else NA_real_,
               within_tumor_range = if (have_t) (v >= min(tv) & v <= max(tv)) else NA)
  }))
}))
write.csv(summ, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_key_fractions_vs_tumors.csv"), row.names = FALSE)

# ---------------- 5. Figure: tumors by cluster with the PDOs overlaid ----------------
grp_levels <- if (have_t) c(CLUSTER_LEVELS, "PDO") else "PDO"
x_labs <- if (have_t) {
  n_cl <- table(frac_t$Cluster)
  c(setNames(paste0("Cluster ", names(n_cl), "\n(n = ", as.integer(n_cl), ")"), names(n_cl)),
    PDO = paste0("PDO\n(n = ", nrow(frac_p), ")"))
} else c(PDO = paste0("PDO\n(n = ", nrow(frac_p), ")"))

mk_panel <- function(f, lab) {
  d_p <- data.frame(Group = "PDO", Value = frac_p[[f]], Sample = frac_p$Sample_ID)
  p <- ggplot()
  if (have_t) {
    d_t <- data.frame(Group = as.character(frac_t$Cluster), Value = frac_t[[f]])
    p <- p +
      geom_boxplot(data = d_t, aes(x = Group, y = Value, fill = Group),
                   width = 0.6, outlier.shape = NA, colour = "grey20") +
      geom_point(data = d_t, aes(x = Group, y = Value),
                 position = position_jitter(width = 0.12, height = 0, seed = SEED),
                 size = 1.5, alpha = 0.6, colour = "grey30")
  }
  p +
    geom_point(data = d_p, aes(x = Group, y = Value, colour = Sample),
               position = position_jitter(width = 0.1, height = 0, seed = SEED), size = 3.3) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_colour_manual(values = pdo_col_map, name = "PDO") +
    scale_x_discrete(limits = grp_levels, labels = x_labs[grp_levels]) +
    labs(title = lab, x = NULL, y = "Estimated fraction") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14))
}
fig_box <- wrap_plots(Map(mk_panel, feat, unname(PLOT_FEATURES)), ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = paste0("Boxes: tumors by RNA cluster (median and interquartile range); grey points: individual tumors. ",
                     "Coloured points: the 4 PDOs. Descriptive only; no test is applied (n = 4)."),
    theme = theme(plot.caption = element_text(size = 11, hjust = 0)))
save_fig(fig_box, "PDO_InstaPrism_PBCA_key_fractions_PDO_vs_tumors", w = 13, h = 11)

# ---------------- 6. Cell-state fractions (Post.ini.cs): PDO versus tumor mean ----------------
tryCatch({
  st_p <- get_states(DECONV_P, frac_p$key)
  message("PDO cell states: ", ncol(st_p), "; row sums range ",
          paste(round(range(rowSums(st_p)), 3), collapse = " to "), " (expected about 1).")

  # Optional: parent cell type of each state from the reference map
  parent <- rep(NA_character_, ncol(st_p)); names(parent) <- colnames(st_p)
  if (file.exists(PBCA_REF_FILE)) {
    try({
      mp    <- readRDS(PBCA_REF_FILE)@map
      st2ct <- unlist(lapply(names(mp), function(ct) setNames(rep(ct, length(mp[[ct]])), mp[[ct]])))
      parent[] <- unname(st2ct[colnames(st_p)])
    }, silent = TRUE)
  }
  # Compartment: named cell types as such; all patient-specific cancer-cell profiles count as Epithelial
  state_comp <- vapply(colnames(st_p), function(s) {
    cp <- comp_of(s)
    if (is.na(cp) && !is.na(parent[s])) cp <- comp_of(parent[s])
    if (is.na(cp) && grepl(TUMOUR_STATE_PATTERN, s)) cp <- "Epithelial"
    cp
  }, character(1))

  # Tumor cell states (same helper), if available
  tum_mean <- tum_median <- rep(NA_real_, ncol(st_p)); names(tum_mean) <- names(tum_median) <- colnames(st_p)
  if (have_t && file.exists(DECONV_T)) {
    st_t <- get_states(DECONV_T, frac_t$key)
    common <- intersect(colnames(st_p), colnames(st_t))
    tum_mean[common]   <- colMeans(st_t[, common, drop = FALSE])
    tum_median[common] <- apply(st_t[, common, drop = FALSE], 2, median)
  } else message("Tumor InstaPrism object not found; tumor columns left empty.")

  state_tab <- data.frame(state = colnames(st_p), parent_cell_type = unname(parent),
                          compartment = unname(state_comp), t(st_p),
                          PDO_mean = colMeans(st_p), PDO_max = apply(st_p, 2, max),
                          tumor_mean = unname(tum_mean), tumor_median = unname(tum_median),
                          check.names = FALSE, row.names = NULL)
  colnames(state_tab)[4:7] <- frac_p$Sample_ID                     # PDO sample names for the 4 fraction columns
  state_tab <- state_tab[order(-state_tab$PDO_mean), ]
  write.csv(state_tab, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_cell_states.csv"), row.names = FALSE)

  # Share of each compartment in each PDO (sum of the state fractions)
  comp_tab <- as.data.frame(t(rowsum(t(st_p), group = ifelse(is.na(state_comp), "Unassigned", state_comp))))
  comp_tab <- cbind(Sample = rownames(comp_tab), comp_tab)
  write.csv(comp_tab, file.path(OUT_DIR, "PDO_InstaPrism_PBCA_compartment_fractions.csv"), row.names = FALSE)
  message("Compartment fractions per PDO:"); print(comp_tab, row.names = FALSE)

  # Figure: top states by mean PDO fraction; columns = 4 PDOs (+ tumor mean)
  top <- head(state_tab, TOP_STATES)
  cols_show <- c(frac_p$Sample_ID, if (have_t && any(!is.na(top$tumor_mean))) "tumor_mean")
  hm <- do.call(rbind, lapply(cols_show, function(cn)
    data.frame(state = top$state, Column = cn, pct = 100 * top[[cn]], stringsAsFactors = FALSE)))
  hm$state  <- factor(hm$state, levels = rev(top$state))
  hm$Column <- factor(hm$Column, levels = cols_show,
                      labels = ifelse(cols_show == "tumor_mean", paste0("Tumor mean\n(n = ", if (have_t) nrow(frac_t) else NA, ")"), cols_show))
  fig_states <- ggplot(hm, aes(x = Column, y = state, fill = pct)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 3.2) +
    scale_fill_gradient(low = "white", high = "#B2182B", name = "Fraction (%)") +
    labs(title = paste0("Top ", nrow(top), " cell states in the PDOs"),
         caption = "Numbers: estimated fraction (%). Tumor mean = mean over all tumors.", x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
          plot.caption = element_text(hjust = 0, size = 9))
  save_fig(fig_states, "PDO_InstaPrism_PBCA_cell_states_top", w = 8, h = 0.3 * nrow(top) + 2.5)
}, error = function(e) {
  message("SKIPPED cell-state analysis: ", conditionMessage(e))
})

message("Done. Outputs in: ", OUT_DIR)
