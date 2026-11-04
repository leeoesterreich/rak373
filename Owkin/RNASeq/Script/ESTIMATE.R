# ==============================================================================
# Script 3: ESTIMATE (immunedeconv) tissue-content check, 8 Primary vs 8 Recurrence
#   A. ESTIMATE scores (Stromal, Immune, ESTIMATE, Tumor purity):
#      TPM = main input, TMM log-CPM = sensitivity check (as in the pILC analysis)
#   B. Paired tests of the ESTIMATE scores
#   C. Is the drop in the ER activity index explained by the change in tissue content?
# Needs the output of script 2 (geneset_composite_scores_per_sample.csv)
# Output folder: .../Okwin/RNASeq/Estimate
# ==============================================================================

suppressPackageStartupMessages({
  library(immunedeconv)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  # tidyverse packages are loaded last so their functions are not masked
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})
select <- dplyr::select   # AnnotationDbi masks dplyr::select
rename <- dplyr::rename   # S4Vectors masks dplyr::rename

SEED <- 42
set.seed(SEED)

# ---- 0. Settings (CHANGE THESE) --------------------------------------------------
BASE        <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/RNASeq"
TPM_FILE    <- file.path(BASE, "counts/matrix/gene_tpm.tsv")
COUNTS_FILE <- file.path(BASE, "counts/matrix/gene_counts.tsv")
META_FILE   <- file.path(BASE, "sample_metadata.csv")            # block, patient, sample_type
OUT_DIR     <- file.path(BASE, "Estimate")                    # results of script 3
SCORES_FILE <- file.path(BASE, "ER_signaling", "geneset_composite_scores_per_sample.csv")   # input from script 2
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EST_NAMES <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
ER_FEATS  <- c("ER_activity_index", "Resistance_index")   # scores adjusted for tissue content

# ---- Helpers ---------------------------------------------------------------------
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# Exact paired Wilcoxon + paired t-test on (Recurrence - Primary); df: patient, timepoint, value
paired_test <- function(df) {
  w <- df %>%
    select(patient, timepoint, value) %>%
    pivot_wider(names_from = timepoint, values_from = value) %>%
    drop_na()
  tibble(
    n_pairs      = nrow(w),
    median_delta = median(w$Recurrence - w$Primary),
    mean_delta   = mean(w$Recurrence - w$Primary),
    n_up         = sum(w$Recurrence > w$Primary),
    p_wilcoxon   = suppressWarnings(wilcox.test(w$Recurrence, w$Primary, paired = TRUE, exact = TRUE)$p.value),
    p_ttest      = t.test(w$Recurrence, w$Primary, paired = TRUE)$p.value
  )
}

# Ensembl IDs (versioned) -> HGNC symbols; duplicated symbols collapsed by maximum (as in pILC)
to_symbols <- function(m) {
  m   <- as.matrix(m)
  storage.mode(m) <- "numeric"
  ens <- sub("\\..*$", "", rownames(m))
  sym <- unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL",
                                      column = "SYMBOL", multiVals = "first"))
  ok  <- !is.na(sym)
  m   <- m[ok, , drop = FALSE]
  sym <- sym[ok]
  dup <- duplicated(sym) | duplicated(sym, fromLast = TRUE)
  m_u <- m[!dup, , drop = FALSE]
  rownames(m_u) <- sym[!dup]
  if (any(dup)) {
    idx <- split(which(dup), sym[dup])
    m_d <- do.call(rbind, lapply(idx, function(i) apply(m[i, , drop = FALSE], 2, max)))
    m_u <- rbind(m_u, m_d)
  }
  m_u
}

# ESTIMATE through immunedeconv (same call as in the pILC analysis); one row per sample
run_estimate <- function(expr) {
  res <- as.data.frame(immunedeconv::deconvolute_estimate(expr))   # scores x samples
  if ("cell_type" %in% colnames(res)) { rownames(res) <- res$cell_type; res$cell_type <- NULL }
  message("  Score labels returned: ", paste(rownames(res), collapse = " | "))

  # Standardize labels (e.g. "Stromal score" -> "StromalScore")
  std <- c(stromalscore = "StromalScore", immunescore = "ImmuneScore",
           estimatescore = "ESTIMATEScore", tumorpurity = "TumorPurity")
  key <- tolower(gsub("[^A-Za-z]", "", rownames(res)))
  rownames(res) <- ifelse(key %in% names(std), std[key], rownames(res))
  missing <- setdiff(EST_NAMES, rownames(res))
  if (length(missing) > 0) stop("Scores not found: ", paste(missing, collapse = ", "))

  res[EST_NAMES, , drop = FALSE] %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("Sample")
}

# ---- 1. Metadata (same 8 pairs as scripts 1 and 2) -------------------------------
meta <- read.csv(META_FILE, stringsAsFactors = FALSE) %>%
  mutate(
    sample    = paste0(block, "_mRNA"),
    timepoint = factor(ifelse(sample_type == "Primary", "Primary", "Recurrence"),
                       levels = c("Primary", "Recurrence"))
  ) %>%
  group_by(patient) %>%
  filter(any(timepoint == "Primary"), any(timepoint == "Recurrence")) %>%
  ungroup() %>%
  mutate(patient = factor(patient)) %>%
  arrange(patient, timepoint)
stopifnot(all(table(meta$patient) == 2))
cat("Pairs used:", nlevels(meta$patient), "\n")

# ---- 2. ESTIMATE: TPM (main) and TMM log-CPM (sensitivity) -----------------------
tpm <- read.delim(TPM_FILE, row.names = 1, check.names = FALSE)
stopifnot(all(meta$sample %in% colnames(tpm)))
message("ESTIMATE on TPM")
est_tpm <- run_estimate(to_symbols(tpm[, meta$sample]))

counts <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE)
stopifnot(all(meta$sample %in% colnames(counts)))
y      <- calcNormFactors(DGEList(counts = as.matrix(counts[, meta$sample])))
logcpm <- cpm(y, log = TRUE, prior.count = 2)
message("ESTIMATE on TMM log-CPM (sensitivity)")
est_lcpm <- run_estimate(to_symbols(logcpm))

shared <- inner_join(est_tpm, est_lcpm, by = "Sample", suffix = c("_TPM", "_logCPM"))
agreement <- tibble(
  Score = EST_NAMES,
  spearman_rho_TPM_vs_logCPM = sapply(EST_NAMES, function(s)
    cor(shared[[paste0(s, "_TPM")]], shared[[paste0(s, "_logCPM")]], method = "spearman"))
)
write.csv(agreement, file.path(OUT_DIR, "ESTIMATE_TPM_vs_logCPM_correlation.csv"), row.names = FALSE)
print(as.data.frame(agreement))

# ---- 3. Join with metadata and the ER scores from script 2 -----------------------
if (!file.exists(SCORES_FILE)) stop("Run script 2 first; not found: ", SCORES_FILE)
er_scores <- read.csv(SCORES_FILE, stringsAsFactors = FALSE) %>%
  select(sample, all_of(ER_FEATS))

scores <- est_tpm %>%
  left_join(meta %>% select(sample, block, patient, timepoint), by = c("Sample" = "sample")) %>%
  left_join(er_scores, by = c("Sample" = "sample")) %>%
  arrange(patient, timepoint)
stopifnot(!anyNA(scores$ER_activity_index), !anyNA(scores$patient))
write.csv(scores, file.path(OUT_DIR, "ESTIMATE_scores_TPM.csv"), row.names = FALSE)
write.csv(est_lcpm, file.path(OUT_DIR, "ESTIMATE_scores_logCPM.csv"), row.names = FALSE)

# ---- 4. Paired tests of the ESTIMATE scores --------------------------------------
long_est <- scores %>%
  select(patient, timepoint, all_of(EST_NAMES)) %>%
  pivot_longer(all_of(EST_NAMES), names_to = "feature", values_to = "value") %>%
  mutate(feature = factor(feature, levels = EST_NAMES))

est_stats <- long_est %>%
  group_by(feature) %>% group_modify(~ paired_test(.x)) %>% ungroup() %>%
  mutate(FDR_wilcoxon = p.adjust(p_wilcoxon, "BH"),
         FDR_ttest    = p.adjust(p_ttest, "BH"))
# ESTIMATEScore = Stromal + Immune and TumorPurity is a function of ESTIMATEScore,
# so these four tests are not independent.
write.csv(est_stats, file.path(OUT_DIR, "ESTIMATE_paired_stats.csv"), row.names = FALSE)
print(as.data.frame(est_stats))

# Result shown on each panel: exact paired Wilcoxon p, number of pairs that go up, median change
est_labels <- est_stats %>%
  mutate(panel = sprintf("%s\nWilcoxon p = %.3f | %d of %d up",
                         feature, p_wilcoxon, n_up, n_pairs)) %>%
  select(feature, panel)
long_est_plot <- long_est %>%
  left_join(est_labels, by = "feature") %>%
  mutate(panel = factor(panel, levels = est_labels$panel))

p_pair <- ggplot(long_est_plot, aes(timepoint, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_line(aes(group = patient, colour = patient), alpha = 0.8) +
  geom_point(aes(colour = patient), size = 2) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Score", colour = "Patient",
       caption = "Exact paired Wilcoxon test, 8 pairs. Exploratory; the four scores are not independent.") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))
save_fig(p_pair, "ESTIMATE_scores_paired", w = 13, h = 4.6)

# ---- 5. Is the ER change explained by the change in tissue content? ---------------
# Per patient: delta = Recurrence - Primary
dl <- scores %>%
  select(patient, timepoint, all_of(c(EST_NAMES, ER_FEATS))) %>%
  pivot_longer(-c(patient, timepoint), names_to = "feature", values_to = "value") %>%
  pivot_wider(names_from = timepoint, values_from = value) %>%
  mutate(delta = Recurrence - Primary) %>%
  select(patient, feature, delta) %>%
  pivot_wider(names_from = feature, values_from = delta)

covs <- c("TumorPurity", "StromalScore", "ImmuneScore")   # one covariate at a time (8 pairs)

# 5a. Spearman correlation of the changes (8 pairs)
delta_cor <- bind_rows(lapply(ER_FEATS, function(o) {
  bind_rows(lapply(covs, function(x) {
    ct <- suppressWarnings(cor.test(dl[[o]], dl[[x]], method = "spearman", exact = FALSE))
    tibble(outcome = o, covariate = x, spearman_rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
write.csv(delta_cor, file.path(OUT_DIR, "ESTIMATE_delta_correlations.csv"), row.names = FALSE)

# 5b. Adjusted change: delta_outcome ~ 1 + delta_covariate. The intercept is the change
#     in the outcome expected if tissue composition did not change (same model as
#     patient + timepoint + covariate in a paired design).
adj <- bind_rows(lapply(ER_FEATS, function(o) {
  bind_rows(lapply(covs, function(x) {
    fit <- lm(dl[[o]] ~ dl[[x]])
    co  <- summary(fit)$coefficients
    ci  <- confint(fit)
    tibble(outcome = o, covariate = x, n_pairs = nrow(dl),
           unadjusted_mean_delta = mean(dl[[o]]),
           adjusted_delta = co[1, 1], adjusted_se = co[1, 2],
           adjusted_ci_low = ci[1, 1], adjusted_ci_high = ci[1, 2],
           p_adjusted = co[1, 4], slope = co[2, 1], p_slope = co[2, 4],
           r_squared = summary(fit)$r.squared)
  }))
}))
write.csv(adj, file.path(OUT_DIR, "ER_change_adjusted_for_ESTIMATE.csv"), row.names = FALSE)
print(as.data.frame(adj))

# 5c. Figure: change in ER index vs change in each ESTIMATE score
sc_long <- dl %>%
  select(patient, ER_activity_index, all_of(covs)) %>%
  pivot_longer(all_of(covs), names_to = "covariate", values_to = "delta_x")
p_sc <- ggplot(sc_long, aes(delta_x, ER_activity_index)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40", formula = y ~ x) +
  geom_point(size = 2.5) +
  geom_text(aes(label = patient), size = 3, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~covariate, scales = "free_x", nrow = 1) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(x = "Change in ESTIMATE score (Recurrence - Primary)",
       y = "Change in ER activity index (Recurrence - Primary)") +
  theme_bw(base_size = 11)
save_fig(p_sc, "ESTIMATE_vs_ER_index_change", w = 12, h = 4.5)

sessionInfo()