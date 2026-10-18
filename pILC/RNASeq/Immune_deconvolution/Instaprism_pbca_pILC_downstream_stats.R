# ============================================================
# Downstream statistics for the pILC InstaPrism (PBCA reference) results
#
# READS saved files only. It does NOT rerun InstaPrism and does NOT overwrite
# any existing InstaPrism / ESTIMATE output; everything new goes to Downstream_stats/.
#
#   1. Cluster (Kruskal-Wallis + pairwise Wilcoxon) and Grade 2 vs 3 (Wilcoxon) statistics
#   2. Boxplots by cluster for the key cell types and summary fractions
#   3. Spearman correlations with clinical variables and ESTIMATE scores
#   4. Cell-state fractions (deconv_res@Post.ini.cs): list + comparison across clusters
#   5. Sensitivity check of plasma-cell-rich samples
# ============================================================

suppressPackageStartupMessages({
  library(InstaPrism)   # only so readRDS() can rebuild the saved S4 InstaPrism object
  library(data.table)   # fread()
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)    # combine panels
})

SEED <- 42
set.seed(SEED)

# ---------------- 1. Paths and settings ----------------
RNA_DIR   <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
INSTA_DIR <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/InstaPrism_PBCA_reference/Stranded")
OUT_DIR   <- file.path(INSTA_DIR, "Downstream_stats")      # new folder, nothing existing is touched
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Saved InstaPrism outputs (read only)
FRAC_FILE   <- file.path(INSTA_DIR, "pILC_InstaPrism_PBCA_fraction_with_summary.csv")
DECONV_FILE <- file.path(INSTA_DIR, "pILC_InstaPrism_PBCA_deconv_result.rds")
ORDER_FILE  <- file.path(INSTA_DIR, "pILC_InstaPrism_PBCA_sample_cluster_order.csv")   # used only as a label check

# Other inputs
CLUSTER_FILE    <- file.path(RNA_DIR, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")  # cluster labels used for ALL analyses (GP numbering)
CLUSTER_K3_FILE <- file.path(RNA_DIR, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3.csv")           # used only as a label check
CLINICAL_FILE   <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/Clinical_correlation/pILC_ClinicoPathological_details.csv"  # note: sits in pILC/, not pILC/RNASeq/
ESTIMATE_FILE   <- file.path(RNA_DIR, "immunedeconv/Stranded/TPM/ESTIMATE_scores.csv")
TPM_FILE        <- file.path(RNA_DIR, "4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt")
PBCA_REF_FILE   <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/PBCA_InstaPrism_reference.Rds")  # optional: state -> cell type names

# Analysis settings
CLUSTER_LEVELS <- c("1", "2", "3")
cluster_cols   <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")   # Dark2, as in the other cluster figures
FDR_ALPHA      <- 0.05
PLASMA_Z       <- 3.5     # plasma-rich = modified z-score above this ...
PLASMA_ABS     <- 0.20    # ... AND plasma fraction at least this (avoids flagging modest values when the MAD is small)
MIN_STATE_MEAN <- 0.005   # cell states with an overall mean fraction below this are not tested
TOP_STATES     <- 25      # states shown in the cell-state heatmap

# Clinical column names are matched by pattern. If the wrong column is picked,
# put the exact header here, e.g. CLIN_COL["grade"] <- "Grade"
CLIN_COL <- c(sample = NA_character_, grade = NA_character_, cellularity = NA_character_,
              tils = NA_character_, pleomorphic = NA_character_)

# Cell types / fractions shown in the boxplots (name in the CSV = label on the plot), grouped by compartment:
# CAF / stromal (row 1-2), epithelial (row 2), TAM and other immune (row 3).
# These 9 also form the "primary" FDR family in the tables.
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
SUMMARY_COLS <- c("Immune_Fraction", "Stromal_Fraction", "Epithelial_Fraction")

# Compartments: which InstaPrism cell types belong to which biological compartment (edit if you define them differently)
COMPARTMENTS <- list(
  CAF              = c("Fibroblasts"),
  Stromal_vascular = c("Endothelial Cells", "Perivascular-like (PVL) Cells"),
  Epithelial       = c("Cancer Epithelial Cells", "Epithelial Cells", "Myoepithelial Cells"),
  TAM              = c("Macrophages"),
  Other_immune     = c("CD8+ T Cells", "CD4+ T Cells", "Regulatory T Cells", "NK Cells", "B Cells",
                       "Plasma Cells", "Monocytes", "Dendritic Cells", "MDSCs", "Mast Cells", "Neutrophils")
)
FOCUS_COMPARTMENTS    <- c("CAF", "Stromal_vascular", "Epithelial", "TAM")   # get cell-state boxplots
STATES_PER_COMPARTMENT <- 6                                                  # top states (by Kruskal-Wallis p) per figure
# Fallback only if the reference map cannot give the parent cell type of each state: match on the state name
COMPARTMENT_PATTERNS <- c(CAF = "fibro|caf", Stromal_vascular = "endothel|perivasc|pvl|pericyte",
                          Epithelial = "epithel|cancer", TAM = "macroph")

# Pre-specified InstaPrism vs ESTIMATE agreement pairs (all expected to be positive)
AGREE_PAIRS <- data.frame(
  feature  = c("Immune_Fraction", "Stromal_Fraction", "Epithelial_Fraction", "Cancer Epithelial Cells"),
  variable = c("ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "ESTIMATE_TumorPurity", "ESTIMATE_TumorPurity"),
  stringsAsFactors = FALSE
)

# Genes used to judge whether plasma-rich samples look like a real plasma-cell infiltrate
PLASMA_GENES <- c("MZB1", "JCHAIN", "SDC1", "TNFRSF17", "XBP1", "PRDM1", "CD38")  # plasma-cell program
HEAT_GENES   <- c(PLASMA_GENES, "IGKC", "IGLC2", "IGHG1", "IGHA1", "IGHM", "MS4A1", "CD79A")  # rows of the marker heatmap

# ---------------- 2. Helper functions ----------------

need <- function(path, what) if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)

# Make sample names comparable across files:
#   TP03-117 / TP03.117_S25 -> TP03_117 ;  TP09-9 / TP09-09 / TP09.9_S34 -> TP09_9
make_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_S[0-9]+$", "", x)          # sequencing suffix (_S25, _S34, ...)
  x <- gsub("[^A-Z0-9]+", "_", x)       # '-' and '.' -> '_'
  x <- gsub("^_|_$", "", x)
  x <- sub("([0-9])T$", "\\1", x)       # trailing tumour letter, e.g. TP09-09T
  gsub("_0+([0-9])", "_\\1", x)         # drop leading zeros: TP09_09 -> TP09_9
}

# p-value text for figures: "= 0.012" / "< 0.001" / "n/a"
fmt_p <- function(p) ifelse(is.na(p), "n/a",
                            ifelse(p < 0.001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3))))

# Compartment of a cell type (NA if not listed) and a helper that adds a 'compartment' column
comp_of <- function(ct) {
  hit <- names(COMPARTMENTS)[vapply(COMPARTMENTS, function(v) ct %in% v, logical(1))]
  if (length(hit) > 0) hit[1] else NA_character_
}
add_comp <- function(df) {
  df$compartment <- ifelse(df$feature %in% SUMMARY_COLS, "Summary fraction",
                           unname(vapply(as.character(df$feature), comp_of, character(1))))
  df
}

# Find a clinical column by pattern (or use the manual override)
find_col <- function(df, pattern, override, label) {
  if (!is.na(override)) {
    if (!override %in% colnames(df)) stop("CLIN_COL value '", override, "' is not a column of the clinical file.")
    return(override)
  }
  hit <- grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) {
    message("  [", label, "] no column matched. Available columns: ", paste(colnames(df), collapse = " | "))
    return(NA_character_)
  }
  if (length(hit) > 1)
    message("  [", label, "] several columns matched (", paste(hit, collapse = ", "),
            "); using the first. Set CLIN_COL to choose another.")
  hit[1]
}

# Numeric parsing: keeps digits and '.', so "30%" -> 30. Warns about values such as "<5" or "10-20".
to_num <- function(x, label) {
  x <- as.character(x)
  if (any(grepl("[<>]|[0-9]\\s*-\\s*[0-9]", x)))
    warning(label, ": some values look like '<5' or ranges; they are parsed as the bare number. Please check.")
  suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
}

# Grade: first digit of the value; only 2 and 3 are kept
to_grade <- function(x) {
  g <- suppressWarnings(as.integer(sub("^\\D*(\\d).*$", "\\1", as.character(x))))
  ifelse(g %in% c(2L, 3L), g, NA_integer_)
}

# Kruskal-Wallis across clusters for one feature, plus medians / means / eta-squared (H-based effect size)
kw_one <- function(v, g) {
  ok <- is.finite(v) & !is.na(g)
  v <- v[ok]; g <- droplevels(g[ok])
  out <- data.frame(n = length(v))
  for (lv in CLUSTER_LEVELS) {
    out[[paste0("n_c", lv)]]      <- sum(g == lv)
    out[[paste0("median_c", lv)]] <- if (any(g == lv)) median(v[g == lv]) else NA_real_
    out[[paste0("mean_c", lv)]]   <- if (any(g == lv)) mean(v[g == lv])   else NA_real_
  }
  testable <- length(unique(v)) > 1 && nlevels(g) >= 2 && all(table(g) >= 3)
  if (testable) {
    kt <- kruskal.test(v ~ g)
    out$KW_H   <- unname(kt$statistic)
    out$KW_p   <- kt$p.value
    out$eta2_H <- max(0, (unname(kt$statistic) - nlevels(g) + 1) / (length(v) - nlevels(g)))
  } else {
    out$KW_H <- NA_real_; out$KW_p <- NA_real_; out$eta2_H <- NA_real_
  }
  out
}

# Pairwise Wilcoxon rank-sum tests between clusters for one feature.
# rank_biserial > 0 means group1 tends to be higher than group2 (range -1 to 1).
pair_one <- function(v, g) {
  ok <- is.finite(v) & !is.na(g)
  v <- v[ok]; g <- as.character(g[ok])
  do.call(rbind, lapply(combn(CLUSTER_LEVELS, 2, simplify = FALSE), function(pr) {
    x <- v[g == pr[1]]; y <- v[g == pr[2]]
    res <- data.frame(group1 = pr[1], group2 = pr[2], n1 = length(x), n2 = length(y),
                      median1 = if (length(x)) median(x) else NA_real_,
                      median2 = if (length(y)) median(y) else NA_real_,
                      W = NA_real_, p = NA_real_, rank_biserial = NA_real_)
    if (length(x) >= 3 && length(y) >= 3 && length(unique(c(x, y))) > 1) {
      wt <- wilcox.test(x, y, exact = FALSE)        # normal approximation (ties are common with fractions)
      res$W <- unname(wt$statistic); res$p <- wt$p.value
      res$rank_biserial <- 2 * unname(wt$statistic) / (length(x) * length(y)) - 1
    }
    res
  }))
}

# Spearman correlation (rank-based, robust to the extreme plasma samples)
spearman_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2)
    return(data.frame(n = sum(ok), rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.frame(n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
}

# Save a figure as PNG + PDF
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
}

# ---------------- 3. Load and harmonize the data ----------------
for (f in c(FRAC_FILE, DECONV_FILE, CLUSTER_FILE, CLINICAL_FILE, ESTIMATE_FILE, TPM_FILE)) need(f, "Input file")

# 3a. InstaPrism fractions (cell types + Immune/Stromal/Epithelial summaries)
frac <- fread(FRAC_FILE, data.table = FALSE, check.names = FALSE)   # keep names such as "CD8+ T Cells"
if (!"Sample_ID" %in% colnames(frac)) stop("Column 'Sample_ID' not found in ", FRAC_FILE)
frac$key <- make_key(frac$Sample_ID)
if (anyDuplicated(frac$key)) stop("Duplicated sample keys in the fraction file.")

cell_cols <- setdiff(colnames(frac), c("Sample_ID", "key", SUMMARY_COLS))
miss_feat <- setdiff(c(names(PLOT_FEATURES), "Cancer Epithelial Cells"), colnames(frac))
if (length(miss_feat) > 0)
  stop("Expected columns not found in the fraction file: ", paste(miss_feat, collapse = ", "),
       "\nColumns present: ", paste(colnames(frac), collapse = " | "))
message("Samples in fraction file: ", nrow(frac), "; cell types: ", length(cell_cols))

# 3b. Cluster membership (GP numbering, from Cluster_membership_k3_GPlabels.csv)
cl_raw <- fread(CLUSTER_FILE, data.table = FALSE)
message("Cluster file columns: ", paste(colnames(cl_raw), collapse = " | "))
s_hit <- grep("^sample$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
s_col <- if (length(s_hit) > 0) s_hit[1] else colnames(cl_raw)[1]
c_hit <- grep("^cluster$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(cl_raw)[ncol(cl_raw)]
message("Using sample column '", s_col, "' and cluster column '", c_col, "'")

cl <- data.frame(key = make_key(cl_raw[[s_col]]),
                 num = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cl_raw[[c_col]])))))
if (anyNA(cl$num) || !all(cl$num %in% 1:3)) stop("Cluster values must be 1, 2 or 3 in ", CLUSTER_FILE)
if (anyDuplicated(cl$key)) stop("Duplicated sample keys in the cluster file.")
cl$Cluster <- factor(as.character(cl$num), levels = CLUSTER_LEVELS)

dat <- left_join(frac, cl[, c("key", "Cluster")], by = "key")
if (anyNA(dat$Cluster))
  stop("No cluster found for: ", paste(dat$Sample_ID[is.na(dat$Cluster)], collapse = ", "))
message("Samples per cluster:"); print(table(dat$Cluster))

# Label checks: Cluster_membership_k3.csv (used by the earlier InstaPrism script) should equal the GP-labels file, and the earlier InstaPrism plots
# should use the same numbering (the plot-order file may predate the relabelling).
label_check <- function(other_file, what) {
  if (!file.exists(other_file)) { message("  (", what, " not found, check skipped)"); return(invisible(NULL)) }
  o <- fread(other_file, data.table = FALSE)
  o_s <- grep("^sample", colnames(o), ignore.case = TRUE, value = TRUE)[1]
  o_c <- grep("^cluster$", colnames(o), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(o_s) || is.na(o_c)) { message("  (", what, ": sample/cluster columns not recognised, check skipped)"); return(invisible(NULL)) }
  chk <- merge(data.frame(key = dat$key, now = as.character(dat$Cluster)),
               data.frame(key = make_key(o[[o_s]]), other = as.character(o[[o_c]])), by = "key")
  message("Cluster numbers: current file (rows) vs ", what, " (columns), ", nrow(chk), " shared samples:")
  print(table(current = chk$now, other = chk$other))
  if (any(chk$now != chk$other))
    warning("Cluster numbering differs from ", what, ". Figures made from that file use different cluster numbers.")
}
label_check(CLUSTER_K3_FILE, "Cluster_membership_k3.csv")
label_check(ORDER_FILE, "the earlier InstaPrism sample_cluster_order.csv")

# 3c. Clinical table
clin <- fread(CLINICAL_FILE, data.table = FALSE, check.names = FALSE)
message("Clinical file columns: ", paste(colnames(clin), collapse = " | "))
col_sample <- find_col(clin, "^sample", CLIN_COL["sample"], "sample ID")
if (is.na(col_sample)) col_sample <- colnames(clin)[1]
col_grade  <- find_col(clin, "grade",       CLIN_COL["grade"],       "grade")
col_cell   <- find_col(clin, "cellularity", CLIN_COL["cellularity"], "tumour cellularity")
col_tils   <- find_col(clin, "^tils?$|^tils?[ _.(%]", CLIN_COL["tils"], "TILs")
col_pleo   <- find_col(clin, "pleomorph",   CLIN_COL["pleomorphic"], "pleomorphic component")
message("Clinical columns used -> sample: ", col_sample, " | grade: ", col_grade, " | cellularity: ", col_cell,
        " | TILs: ", col_tils, " | pleomorphic: ", col_pleo)

clin_sub <- data.frame(
  key                   = make_key(clin[[col_sample]]),
  Grade                 = if (!is.na(col_grade)) to_grade(clin[[col_grade]]) else NA_integer_,
  Tumor_cellularity     = if (!is.na(col_cell))  to_num(clin[[col_cell]], "Tumor cellularity") else NA_real_,
  TILs                  = if (!is.na(col_tils))  to_num(clin[[col_tils]], "TILs") else NA_real_,
  Pleomorphic_component = if (!is.na(col_pleo))  to_num(clin[[col_pleo]], "Pleomorphic component") else NA_real_
)
if (anyDuplicated(clin_sub$key))
  stop("Duplicated sample keys in the clinical file: ", paste(clin_sub$key[duplicated(clin_sub$key)], collapse = ", "))
dat <- left_join(dat, clin_sub, by = "key")
no_clin <- dat$Sample_ID[!dat$key %in% clin_sub$key]
if (length(no_clin) > 0) message("Samples not found in the clinical file: ", paste(no_clin, collapse = ", "))
message("Clinical values available (non-missing out of ", nrow(dat), "):")
print(colSums(!is.na(dat[, c("Grade", "Tumor_cellularity", "TILs", "Pleomorphic_component")])))

# 3d. ESTIMATE scores (TPM-based, from the corrected counts)
est <- fread(ESTIMATE_FILE, data.table = FALSE)
if (!all(c("Sample", "ImmuneScore", "StromalScore", "TumorPurity") %in% colnames(est)))
  stop("ESTIMATE file needs the columns Sample, ImmuneScore, StromalScore, TumorPurity. Found: ",
       paste(colnames(est), collapse = " | "))
est$key <- make_key(est$Sample)
est_sub <- est[, c("key", "ImmuneScore", "StromalScore", "TumorPurity")]
colnames(est_sub)[2:4] <- paste0("ESTIMATE_", colnames(est_sub)[2:4])
dat <- left_join(dat, est_sub, by = "key")
if (anyNA(dat$ESTIMATE_ImmuneScore))
  warning("No ESTIMATE score for: ", paste(dat$Sample_ID[is.na(dat$ESTIMATE_ImmuneScore)], collapse = ", "))
if ("Cluster" %in% colnames(est)) {
  e_cl <- as.character(est$Cluster[match(dat$key, est$key)])
  if (any(e_cl != as.character(dat$Cluster), na.rm = TRUE))
    warning("Cluster numbers in the ESTIMATE file differ from the current cluster file for ",
            sum(e_cl != as.character(dat$Cluster), na.rm = TRUE), " samples.")
}

# ---------------- 4. Flag plasma-rich samples (used in several sections) ----------------
pl <- dat[["Plasma Cells"]]
dat$Plasma_modified_z <- if (mad(pl) > 0) (pl - median(pl)) / mad(pl) else NA_real_   # robust z: (x - median) / MAD
dat$Plasma_rich <- pl >= PLASMA_ABS & (is.na(dat$Plasma_modified_z) | dat$Plasma_modified_z > PLASMA_Z)
message("Plasma-rich samples (fraction >= ", PLASMA_ABS, " and modified z > ", PLASMA_Z, "): ",
        if (any(dat$Plasma_rich)) paste(dat$Sample_ID[dat$Plasma_rich], collapse = ", ") else "none")

features <- c(SUMMARY_COLS, cell_cols)   # every feature tested
feature_type <- function(f) ifelse(f %in% SUMMARY_COLS, "Summary fraction", "Cell type")

# ---------------- 5. Statistics: clusters (KW + pairwise Wilcoxon) and Grade 2 vs 3 ----------------
# FDR families (BH):
#   * "all features"  = every cell type + the 3 summary fractions (conservative, shown on the plots)
#   * "primary"       = the 9 pre-specified boxplot features (PLOT_FEATURES)
# Pairwise tests are only to be interpreted when the Kruskal-Wallis FDR (all features) is < FDR_ALPHA.

kw_tab <- do.call(rbind, lapply(features, function(f) cbind(feature = f, kw_one(dat[[f]], dat$Cluster))))
kw_tab$feature_type <- feature_type(kw_tab$feature)
kw_tab <- add_comp(kw_tab)
kw_tab$FDR_BH_all_features <- p.adjust(kw_tab$KW_p, method = "BH")
prim <- kw_tab$feature %in% names(PLOT_FEATURES)
kw_tab$FDR_BH_primary_features <- NA_real_
kw_tab$FDR_BH_primary_features[prim] <- p.adjust(kw_tab$KW_p[prim], method = "BH")
kw_tab <- kw_tab[order(kw_tab$KW_p), ]
write.csv(kw_tab, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cluster_stats.csv"), row.names = FALSE)

pw_tab <- do.call(rbind, lapply(features, function(f) cbind(feature = f, pair_one(dat[[f]], dat$Cluster))))
pw_tab$feature_type <- feature_type(pw_tab$feature)
pw_tab <- add_comp(pw_tab)
pw_tab <- pw_tab %>%
  group_by(feature) %>%
  mutate(p_adj_BH_within_feature = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  as.data.frame()
pw_tab$FDR_BH_all_pairwise_tests <- p.adjust(pw_tab$p, method = "BH")
pw_tab$KW_FDR_all_features <- kw_tab$FDR_BH_all_features[match(pw_tab$feature, kw_tab$feature)]
pw_tab$pairwise_interpretable <- !is.na(pw_tab$KW_FDR_all_features) & pw_tab$KW_FDR_all_features < FDR_ALPHA
pw_tab <- pw_tab[order(pw_tab$p), ]
write.csv(pw_tab, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cluster_pairwise_wilcoxon.csv"), row.names = FALSE)

message("Top cluster differences (Kruskal-Wallis):")
print(head(kw_tab[, c("feature", "KW_p", "FDR_BH_all_features", "FDR_BH_primary_features", "eta2_H")], 10), row.names = FALSE)

# Grade 2 vs Grade 3 (Wilcoxon rank-sum). rank_biserial > 0 means higher in Grade 3.
has_grade <- sum(!is.na(dat$Grade)) >= 6 && all(c(2L, 3L) %in% dat$Grade)
if (has_grade) {
  message("Grade counts:"); print(table(Grade = dat$Grade, useNA = "ifany"))
  message("Grade by cluster (grade may be confounded with cluster):")
  print(table(Cluster = dat$Cluster, Grade = dat$Grade))
  gr_tab <- do.call(rbind, lapply(features, function(f) {
    v <- dat[[f]]; g <- dat$Grade
    x <- v[!is.na(g) & g == 3L & is.finite(v)]; y <- v[!is.na(g) & g == 2L & is.finite(v)]
    res <- data.frame(feature = f, feature_type = feature_type(f), n_grade2 = length(y), n_grade3 = length(x),
                      median_grade2 = if (length(y)) median(y) else NA_real_,
                      median_grade3 = if (length(x)) median(x) else NA_real_,
                      W = NA_real_, p = NA_real_, rank_biserial = NA_real_)
    if (length(x) >= 3 && length(y) >= 3 && length(unique(c(x, y))) > 1) {
      wt <- wilcox.test(x, y, exact = FALSE)
      res$W <- unname(wt$statistic); res$p <- wt$p.value
      res$rank_biserial <- 2 * unname(wt$statistic) / (length(x) * length(y)) - 1
    }
    res
  }))
  gr_tab <- add_comp(gr_tab)
  gr_tab$FDR_BH_all_features <- p.adjust(gr_tab$p, method = "BH")
  prim_g <- gr_tab$feature %in% names(PLOT_FEATURES)
  gr_tab$FDR_BH_primary_features <- NA_real_
  gr_tab$FDR_BH_primary_features[prim_g] <- p.adjust(gr_tab$p[prim_g], method = "BH")
  gr_tab <- gr_tab[order(gr_tab$p), ]
  write.csv(gr_tab, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_grade_stats.csv"), row.names = FALSE)
  message("Top Grade 2 vs 3 differences:")
  print(head(gr_tab[, c("feature", "p", "FDR_BH_all_features", "rank_biserial")], 10), row.names = FALSE)
} else {
  message("SKIPPED grade statistics: no usable Grade 2 / Grade 3 values were parsed (see the clinical columns above).")
}

# ---------------- 6. Boxplots by cluster ----------------
n_cl   <- table(dat$Cluster)
x_labs <- setNames(paste0("Cluster ", names(n_cl), "\n(n = ", as.integer(n_cl), ")"), names(n_cl))

# One cluster boxplot for any vector of values (used for cell types and for cell states)
box_core <- function(v, title, subtitle) {
  d <- data.frame(Cluster = dat$Cluster, Value = v)
  ggplot(d, aes(x = Cluster, y = Value, fill = Cluster)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, colour = "grey20") +
    geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 1.8, alpha = 0.7, colour = "black") +
    scale_fill_manual(values = cluster_cols) +
    scale_x_discrete(labels = x_labs) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Estimated fraction") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11))
}
make_box <- function(feat, lab) {
  fdr <- kw_tab$FDR_BH_all_features[kw_tab$feature == feat]
  box_core(dat[[feat]], lab, paste("Kruskal-Wallis FDR", fmt_p(fdr)))
}
box_plots <- Map(make_box, names(PLOT_FEATURES), unname(PLOT_FEATURES))
fig_box <- wrap_plots(box_plots, ncol = 3) +
  plot_annotation(
    caption = paste0("Boxes show the median and interquartile range; points are individual samples.\n",
                     "Kruskal-Wallis test across the 3 RNA clusters; FDR = Benjamini-Hochberg over all ",
                     nrow(kw_tab), " cell types and summary fractions."),
    theme = theme(plot.caption = element_text(size = 11, hjust = 0)))
save_fig(fig_box, "pILC_InstaPrism_PBCA_key_fractions_by_cluster", w = 13, h = 11)

# ---------------- 7. Spearman correlations: clinical variables and ESTIMATE ----------------
corr_vars <- c("Tumor_cellularity", "TILs", "Pleomorphic_component",
               "ESTIMATE_ImmuneScore", "ESTIMATE_StromalScore", "ESTIMATE_TumorPurity")
corr_vars <- corr_vars[sapply(corr_vars, function(v) sum(is.finite(dat[[v]])) >= 8)]
message("Variables used for correlation: ", paste(corr_vars, collapse = ", "))

cor_tab <- do.call(rbind, lapply(features, function(f)
  do.call(rbind, lapply(corr_vars, function(v) cbind(feature = f, variable = v, spearman_one(dat[[f]], dat[[v]]))))))
cor_tab$variable_type <- ifelse(grepl("^ESTIMATE_", cor_tab$variable), "ESTIMATE", "Clinical")
cor_tab$FDR_BH_all_tests <- p.adjust(cor_tab$p, method = "BH")      # exploratory family: every feature x variable
pair_key <- paste(AGREE_PAIRS$feature, AGREE_PAIRS$variable)
cor_tab$prespecified_agreement_pair <- paste(cor_tab$feature, cor_tab$variable) %in% pair_key
cor_tab$FDR_BH_prespecified_pairs <- NA_real_
idx <- cor_tab$prespecified_agreement_pair
cor_tab$FDR_BH_prespecified_pairs[idx] <- p.adjust(cor_tab$p[idx], method = "BH")
cor_tab <- cor_tab[order(cor_tab$p), ]
write.csv(cor_tab, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_spearman_correlations.csv"), row.names = FALSE)
message("InstaPrism vs ESTIMATE (pre-specified pairs, expected positive):")
print(cor_tab[cor_tab$prespecified_agreement_pair, c("feature", "variable", "n", "rho", "p", "FDR_BH_prespecified_pairs")],
      row.names = FALSE)

# 7a. Agreement scatterplots (plasma-rich samples shown as triangles)
agree_plots <- lapply(seq_len(nrow(AGREE_PAIRS)), function(i) {
  f <- AGREE_PAIRS$feature[i]; v <- AGREE_PAIRS$variable[i]
  if (!v %in% corr_vars) return(NULL)
  d  <- data.frame(x = dat[[v]], y = dat[[f]], Cluster = dat$Cluster,
                   Group = ifelse(dat$Plasma_rich, "plasma-rich", "other"))
  st <- spearman_one(d$x, d$y)
  ggplot(d, aes(x = x, y = y, colour = Cluster, shape = Group)) +
    geom_point(size = 2.6, alpha = 0.85) +
    scale_colour_manual(values = cluster_cols, name = "RNA cluster", labels = paste("Cluster", CLUSTER_LEVELS)) +
    scale_shape_manual(values = c("other" = 16, "plasma-rich" = 17), name = NULL) +
    labs(title = paste(gsub("_", " ", f), "vs", sub("ESTIMATE_", "ESTIMATE ", v)),
         subtitle = paste0("Spearman rho = ", formatC(st$rho, format = "f", digits = 2), ", p ", fmt_p(st$p),
                           " (n = ", st$n, ")"),
         x = sub("ESTIMATE_", "ESTIMATE ", v), y = paste(gsub("_", " ", f), "(InstaPrism)")) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12.5), plot.subtitle = element_text(size = 11))
})
agree_plots <- Filter(Negate(is.null), agree_plots)
if (length(agree_plots) > 0) {
  fig_agree <- wrap_plots(agree_plots, ncol = 2) + plot_layout(guides = "collect") +
    plot_annotation(caption = "Triangles: plasma-rich samples (see the plasma sensitivity figure).",
                    theme = theme(plot.caption = element_text(size = 10.5, hjust = 0)))
  save_fig(fig_agree, "pILC_InstaPrism_PBCA_InstaPrism_vs_ESTIMATE_agreement", w = 11, h = 9)
}

# 7b. Heatmap of the full Spearman matrix (stars: FDR over all tests < FDR_ALPHA)
feat_order <- c(SUMMARY_COLS, cell_cols[order(-colMeans(dat[, cell_cols, drop = FALSE]))])
cor_plot <- cor_tab %>%
  mutate(feature  = factor(feature, levels = rev(feat_order)),
         variable = factor(variable, levels = corr_vars),
         star     = ifelse(!is.na(FDR_BH_all_tests) & FDR_BH_all_tests < FDR_ALPHA, "*", ""))
fig_cor <- ggplot(cor_plot, aes(x = variable, y = feature, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = star), size = 6, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1),
                       na.value = "grey90", name = "Spearman rho") +
  labs(title = "InstaPrism fractions vs clinical variables and ESTIMATE scores",
       caption = paste0("* FDR < ", FDR_ALPHA, " (Benjamini-Hochberg over all ", nrow(cor_tab), " correlations)."),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))
save_fig(fig_cor, "pILC_InstaPrism_PBCA_spearman_heatmap", w = 8, h = 8.5)

# ---------------- 8. Cell-state fractions (deconv_res@Post.ini.cs) ----------------
# Wrapped in tryCatch: if the object does not have the expected structure the script reports why and continues.
tryCatch({
  deconv_res <- readRDS(DECONV_FILE)
  slots <- if (isS4(deconv_res)) slotNames(deconv_res) else names(deconv_res)
  message("Slots / elements of the saved InstaPrism object: ", paste(slots, collapse = ", "))
  if (!"Post.ini.cs" %in% slots) stop("no 'Post.ini.cs' element in the saved object")

  cs_obj <- if (isS4(deconv_res)) slot(deconv_res, "Post.ini.cs") else deconv_res[["Post.ini.cs"]]
  theta  <- if (isS4(cs_obj) && "theta" %in% slotNames(cs_obj)) slot(cs_obj, "theta") else cs_obj
  theta  <- as.matrix(theta)
  if (is.null(rownames(theta)) || is.null(colnames(theta))) stop("state matrix has no row/column names")

  # Orientation: samples must be rows. Decide by which dimension matches our sample keys.
  hit_rows <- sum(make_key(rownames(theta)) %in% dat$key)
  hit_cols <- sum(make_key(colnames(theta)) %in% dat$key)
  st_mat <- if (hit_cols >= hit_rows) t(theta) else theta
  message("State matrix (as stored): ", nrow(theta), " x ", ncol(theta), "; samples matched on ",
          if (hit_cols >= hit_rows) "columns (transposed)" else "rows", ".")
  rownames(st_mat) <- make_key(rownames(st_mat))
  if (!all(dat$key %in% rownames(st_mat))) stop("not all samples are present in the state matrix")
  st_mat <- st_mat[dat$key, , drop = FALSE]
  message("Cell states: ", ncol(st_mat), "; row sums range ",
          paste(round(range(rowSums(st_mat)), 3), collapse = " to "), " (expected about 1).")

  # Optional: parent cell type of each state, taken from the PBCA reference map
  parent <- rep(NA_character_, ncol(st_mat)); names(parent) <- colnames(st_mat)
  if (file.exists(PBCA_REF_FILE)) {
    try({
      ref <- readRDS(PBCA_REF_FILE)
      mp  <- ref@map
      st2ct <- unlist(lapply(names(mp), function(ct) setNames(rep(ct, length(mp[[ct]])), mp[[ct]])))
      parent[] <- unname(st2ct[colnames(st_mat)])
    }, silent = TRUE)
  }
  if (all(is.na(parent))) message("Parent cell types could not be mapped from the reference (column left empty).")

  # Compartment of each state. In the PBCA reference the 17 named cell types are the only non-tumour states;
  # all other columns are patient-specific tumour-cell profiles (e.g. "TNBC_Patient 1", "ER+_CID3941"), i.e. epithelial.
  TUMOUR_STATE_PATTERN <- "TNBC|HER2|ER\\+|Luminal|_BC_|CID|Patient"
  state_comp <- vapply(colnames(st_mat), function(st_) {
    cp <- comp_of(st_)                                         # state name is itself a cell type
    if (is.na(cp) && !is.na(parent[st_])) cp <- comp_of(parent[st_])   # else use the reference map, if it worked
    if (is.na(cp) && grepl(TUMOUR_STATE_PATTERN, st_)) cp <- "Epithelial"
    cp
  }, character(1))
  names(state_comp) <- colnames(st_mat)
  message("States per compartment:"); print(table(state_comp, useNA = "ifany"))

  # List of states with overall abundance
  state_list <- data.frame(
    state = colnames(st_mat), parent_cell_type = unname(parent), compartment = unname(state_comp),
    mean_fraction   = colMeans(st_mat),
    median_fraction = apply(st_mat, 2, median),
    max_fraction    = apply(st_mat, 2, max),
    n_samples_ge_1pct = colSums(st_mat >= 0.01),
    tested_across_clusters = colMeans(st_mat) >= MIN_STATE_MEAN,
    row.names = NULL
  )
  state_list <- state_list[order(-state_list$mean_fraction), ]
  write.csv(state_list, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cell_states_list.csv"), row.names = FALSE)
  message(nrow(state_list), " states listed; ", sum(state_list$tested_across_clusters),
          " have mean fraction >= ", MIN_STATE_MEAN, " and are tested across clusters.")

  # Kruskal-Wallis across clusters for the tested states (filter does not use cluster labels)
  test_states <- state_list$state[state_list$tested_across_clusters]
  if (length(test_states) >= 2) {
    st_tab <- do.call(rbind, lapply(test_states, function(s) cbind(state = s, kw_one(st_mat[, s], dat$Cluster))))
    st_tab$parent_cell_type <- unname(parent[st_tab$state])
    st_tab$FDR_BH_states <- p.adjust(st_tab$KW_p, method = "BH")                 # over all tested states
    st_tab$compartment <- unname(state_comp[st_tab$state])
    st_tab$FDR_BH_within_compartment <- NA_real_                                   # over the tested states of the same compartment
    for (cp in unique(na.omit(st_tab$compartment))) {
      ii <- which(st_tab$compartment == cp)
      st_tab$FDR_BH_within_compartment[ii] <- p.adjust(st_tab$KW_p[ii], method = "BH")
    }
    st_tab <- st_tab[order(st_tab$KW_p), ]
    write.csv(st_tab, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_cell_states_cluster_stats.csv"), row.names = FALSE)

    # Boxplots of the top states (smallest Kruskal-Wallis p) in each focus compartment
    for (cp in FOCUS_COMPARTMENTS) {
      cs_all <- st_tab[!is.na(st_tab$compartment) & st_tab$compartment == cp, ]
      if (nrow(cs_all) < 2) { message("  Compartment ", cp, " has ", nrow(cs_all),
                                      " state(s) (= the cell type itself, already in the main boxplots); figure skipped."); next }
      cs <- head(cs_all[order(cs_all$KW_p), ], STATES_PER_COMPARTMENT)
      state_plots <- lapply(seq_len(nrow(cs)), function(k)
        box_core(st_mat[, cs$state[k]], cs$state[k],
                 paste0("Kruskal-Wallis FDR ", fmt_p(cs$FDR_BH_within_compartment[k]),
                        if (!is.na(cs$parent_cell_type[k])) paste0("  |  ", cs$parent_cell_type[k]) else "")))
      fig_cp <- wrap_plots(state_plots, ncol = 3) +
        plot_annotation(
          title = paste0(cp, " cell states"),
          caption = paste0("Top ", nrow(cs), " of ", nrow(cs_all), " tested states in this compartment (smallest Kruskal-Wallis p). ",
                           "FDR = Benjamini-Hochberg within the compartment."),
          theme = theme(plot.title = element_text(face = "bold", size = 16),
                        plot.caption = element_text(size = 10.5, hjust = 0)))
      save_fig(fig_cp, paste0("pILC_InstaPrism_PBCA_cell_states_", cp, "_by_cluster"),
               w = 13, h = 4.3 * ceiling(nrow(cs) / 3) + 1)
    }
    message("States with Kruskal-Wallis FDR < ", FDR_ALPHA, ": ", sum(st_tab$FDR_BH_states < FDR_ALPHA, na.rm = TRUE))

    # Heatmap: cluster means of the top states, scaled per state (z-score across the 3 clusters)
    top <- head(st_tab, TOP_STATES)
    hm <- do.call(rbind, lapply(seq_len(nrow(top)), function(i) {
      mu <- unname(sapply(CLUSTER_LEVELS, function(lv) mean(st_mat[dat$Cluster == lv, top$state[i]])))
      z  <- if (sd(mu) > 0) (mu - mean(mu)) / sd(mu) else rep(0, length(mu))
      data.frame(state = top$state[i], Cluster = CLUSTER_LEVELS, mean_frac = mu, z = z,
                 label = paste0(top$state[i], ifelse(!is.na(top$FDR_BH_states[i]) & top$FDR_BH_states[i] < FDR_ALPHA, "  *", "")),
                 stringsAsFactors = FALSE)
    }))
    hm$label <- factor(hm$label, levels = rev(unique(hm$label)))
    fig_states <- ggplot(hm, aes(x = Cluster, y = label, fill = z)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.1f", 100 * mean_frac)), size = 3.2) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "Scaled mean\n(z across\nclusters)") +
      scale_x_discrete(labels = x_labs) +
      labs(title = paste0("Top ", nrow(top), " cell states by Kruskal-Wallis p-value"),
           caption = "Numbers: mean fraction (%) in the cluster; colour: scaled per state. * = Kruskal-Wallis FDR < 0.05 (BH over tested states).",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0, size = 9))
    save_fig(fig_states, "pILC_InstaPrism_PBCA_cell_states_top_by_cluster", w = 9, h = 0.3 * nrow(top) + 2.5)
  } else {
    message("Fewer than 2 states pass the abundance filter; cluster comparison of states skipped.")
  }
}, error = function(e) {
  message("SKIPPED cell-state analysis: ", conditionMessage(e))
})

# ---------------- 9. Sensitivity check of plasma-rich samples ----------------
# Question: are the plasma-rich samples a real plasma-cell infiltrate or an artefact?
# Evidence from the TPM matrix: plasma marker genes, immunoglobulin share, kappa/lambda balance,
# agreement with ESTIMATE / TILs; then repeat the cluster tests and agreement without the flagged samples.
if (!any(dat$Plasma_rich)) {
  message("SKIPPED plasma sensitivity: no sample passes the plasma-rich rule.")
} else {
  tpm <- fread(TPM_FILE, data.table = FALSE)
  tpm_m <- as.matrix(tpm[, -1]); rownames(tpm_m) <- tpm[[1]]
  colnames(tpm_m) <- make_key(colnames(tpm_m))
  if (!all(dat$key %in% colnames(tpm_m))) stop("Not all samples are present in the TPM matrix.")
  tpm_m <- tpm_m[, dat$key, drop = FALSE]
  storage.mode(tpm_m) <- "numeric"

  have_pl <- intersect(PLASMA_GENES, rownames(tpm_m))
  message("Plasma marker genes found: ", paste(have_pl, collapse = ", "),
          if (length(setdiff(PLASMA_GENES, have_pl)) > 0) paste0(" | missing: ", paste(setdiff(PLASMA_GENES, have_pl), collapse = ", ")) else "")
  if (length(have_pl) < 3) stop("Fewer than 3 plasma marker genes found in the TPM matrix; check the gene symbols.")

  # Plasma marker score = mean of per-gene z-scores of log2(TPM + 1) across the 31 samples
  lg <- log2(tpm_m[have_pl, , drop = FALSE] + 1)
  zs <- t(scale(t(lg))); zs[!is.finite(zs)] <- 0
  dat$Plasma_marker_score <- colMeans(zs)

  # Immunoglobulin share of the TPM in the matrix and kappa fraction of light chains
  ig_genes  <- grep("^(IG[HKL][VDJC]|IGLL[0-9])", rownames(tpm_m), value = TRUE)
  dat$IG_TPM_share <- colSums(tpm_m[ig_genes, , drop = FALSE]) / colSums(tpm_m)
  lam_genes <- intersect(c("IGLC1", "IGLC2", "IGLC3", "IGLC7"), rownames(tpm_m))
  kap <- if ("IGKC" %in% rownames(tpm_m)) tpm_m["IGKC", ] else rep(NA_real_, ncol(tpm_m))
  lam <- if (length(lam_genes) > 0) colSums(tpm_m[lam_genes, , drop = FALSE]) else rep(NA_real_, ncol(tpm_m))
  dat$Kappa_fraction_of_light_chain <- kap / (kap + lam)
  message("Immunoglobulin genes in the TPM matrix: ", length(ig_genes))

  dat$Immune_fraction_without_plasma <- ifelse(pl < 1, (dat$Immune_Fraction - pl) / (1 - pl), NA_real_)
  dat$Rank_plasma_fraction <- rank(-pl, ties.method = "min")
  dat$Rank_plasma_marker_score <- rank(-dat$Plasma_marker_score, ties.method = "min")
  dat$Rank_IG_TPM_share <- rank(-dat$IG_TPM_share, ties.method = "min")

  ev_cols <- c("Sample_ID", "Cluster", "Plasma Cells", "Plasma_modified_z", "Plasma_rich",
               "Plasma_marker_score", "IG_TPM_share", "Kappa_fraction_of_light_chain",
               "Rank_plasma_fraction", "Rank_plasma_marker_score", "Rank_IG_TPM_share",
               "Immune_Fraction", "Immune_fraction_without_plasma",
               "ESTIMATE_ImmuneScore", "TILs", "Tumor_cellularity")
  ev_cols <- ev_cols[ev_cols %in% colnames(dat)]
  evidence <- dat[order(-pl), ev_cols]
  write.csv(evidence, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_plasma_samples_evidence.csv"), row.names = FALSE)
  message("Evidence for the flagged samples (rank 1 = highest of ", nrow(dat), "):")
  print(evidence[evidence$Plasma_rich, c("Sample_ID", "Plasma Cells", "Plasma_marker_score",
                                         "Rank_plasma_marker_score", "IG_TPM_share", "Rank_IG_TPM_share")],
        row.names = FALSE)

  # Leave-flagged-out: cluster Kruskal-Wallis for the primary (boxplot) features, and the 4 agreement pairs
  keep <- !dat$Plasma_rich
  sens_kw <- do.call(rbind, lapply(names(PLOT_FEATURES), function(f) {
    a <- kw_one(dat[[f]], dat$Cluster); b <- kw_one(dat[[f]][keep], dat$Cluster[keep])
    data.frame(block = "cluster_KW", feature = f, variable = NA_character_,
               n_all = a$n, n_excluding_flagged = b$n,
               stat_all = a$KW_H, p_all = a$KW_p, stat_excluding_flagged = b$KW_H, p_excluding_flagged = b$KW_p)
  }))
  sens_ag <- do.call(rbind, lapply(seq_len(nrow(AGREE_PAIRS)), function(i) {
    f <- AGREE_PAIRS$feature[i]; v <- AGREE_PAIRS$variable[i]
    if (!v %in% corr_vars) return(NULL)
    a <- spearman_one(dat[[f]], dat[[v]]); b <- spearman_one(dat[[f]][keep], dat[[v]][keep])
    data.frame(block = "agreement_spearman", feature = f, variable = v,
               n_all = a$n, n_excluding_flagged = b$n,
               stat_all = a$rho, p_all = a$p, stat_excluding_flagged = b$rho, p_excluding_flagged = b$p)
  }))
  sens <- rbind(sens_kw, sens_ag)
  sens$FDR_all <- NA_real_; sens$FDR_excluding_flagged <- NA_real_
  for (b in unique(sens$block)) {            # BH within each block (cluster tests; agreement pairs)
    i <- sens$block == b
    sens$FDR_all[i] <- p.adjust(sens$p_all[i], method = "BH")
    sens$FDR_excluding_flagged[i] <- p.adjust(sens$p_excluding_flagged[i], method = "BH")
  }
  write.csv(sens, file.path(OUT_DIR, "pILC_InstaPrism_PBCA_plasma_leave_flagged_out.csv"), row.names = FALSE)
  message("Leave-flagged-out comparison (stat = H for cluster_KW, rho for agreement):")
  print(sens[, c("block", "feature", "variable", "n_all", "n_excluding_flagged",
                 "stat_all", "stat_excluding_flagged", "FDR_all", "FDR_excluding_flagged")], row.names = FALSE)

  # Figure: (A) plasma fraction vs plasma marker score, (B) vs IG share, (C) marker heatmap ordered by plasma fraction
  sc_d <- data.frame(Sample_ID = dat$Sample_ID, plasma = pl, marker = dat$Plasma_marker_score,
                     ig = 100 * dat$IG_TPM_share, Cluster = dat$Cluster,
                     Group = ifelse(dat$Plasma_rich, "plasma-rich", "other"))
  mk_scatter <- function(yvar, ylab, title) {
    st <- spearman_one(sc_d$plasma, sc_d[[yvar]])
    ggplot(sc_d, aes(x = plasma, y = .data[[yvar]], colour = Cluster, shape = Group)) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_text(data = sc_d[sc_d$Group == "plasma-rich", ], aes(label = Sample_ID),
                size = 3, vjust = -0.9, show.legend = FALSE) +
      scale_colour_manual(values = cluster_cols, name = "RNA cluster", labels = paste("Cluster", CLUSTER_LEVELS)) +
      scale_shape_manual(values = c("other" = 16, "plasma-rich" = 17), name = NULL) +
      labs(title = title,
           subtitle = paste0("Spearman rho = ", formatC(st$rho, format = "f", digits = 2), ", p ", fmt_p(st$p)),
           x = "InstaPrism plasma-cell fraction", y = ylab) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 12.5))
  }
  pA <- mk_scatter("marker", "Plasma marker score (mean z of log2 TPM+1)", "A. Plasma fraction vs plasma marker genes")
  pB <- mk_scatter("ig", "Immunoglobulin share of TPM (%)", "B. Plasma fraction vs immunoglobulin share")

  heat_genes <- intersect(HEAT_GENES, rownames(tpm_m))
  hz <- t(scale(t(log2(tpm_m[heat_genes, , drop = FALSE] + 1)))); hz[!is.finite(hz)] <- 0
  hz <- as.data.frame(hz); hz$Gene <- rownames(hz)
  hl <- pivot_longer(hz, cols = -Gene, names_to = "key", values_to = "z")
  ord_keys <- dat$key[order(-pl)]
  lab_map  <- setNames(paste0(dat$Sample_ID, ifelse(dat$Plasma_rich, " *", "")), dat$key)
  hl$Sample <- factor(lab_map[hl$key], levels = lab_map[ord_keys])
  hl$Gene   <- factor(hl$Gene, levels = rev(heat_genes))
  pH <- ggplot(hl, aes(x = Sample, y = Gene, fill = z)) +
    geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "z-score of\nlog2(TPM+1)") +
    labs(title = "C. Marker genes, samples ordered by plasma fraction (high to low); * = plasma-rich",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 12.5))

  fig_plasma <- ((pA | pB) + plot_layout(guides = "collect")) / pH + plot_layout(heights = c(1, 1.1))
  save_fig(fig_plasma, "pILC_InstaPrism_PBCA_plasma_sensitivity", w = 15, h = 10)
}

message("Done. Outputs in: ", OUT_DIR)