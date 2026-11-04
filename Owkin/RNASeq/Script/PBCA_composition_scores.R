# ==============================================================================
# Script 4: tissue composition from PBCA-derived gene lists (epithelial / stromal / immune)
#   1. Rank genes by compartment specificity in the PBCA InstaPrism reference
#      (same logic as the pILC cancer-epithelial script)
#   2. Remove ER-axis genes so the lists are independent of the ER scores
#   3. Composite score per compartment = mean z-score of log2(TPM+1) of its top-N genes
#   4. Tumor content index = z(Epithelial) - mean(z(Stromal), z(Immune))
#   5. Paired tests (8 pairs) and ER-score change adjusted for composition
# Needs the output of script 2 (and, optionally, script 3 for the ESTIMATE check)
# Output folder: .../Okwin/RNASeq/Deconvolution
# ==============================================================================

suppressPackageStartupMessages({
  library(InstaPrism)      # only so readRDS() can rebuild the S4 reference object
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readxl)
  library(msigdbr)
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
BASE          <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/RNASeq"
TPM_FILE      <- file.path(BASE, "counts/matrix/gene_tpm.tsv")
META_FILE     <- file.path(BASE, "sample_metadata.csv")            # block, patient, sample_type
PBCA_REF_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Immune_deconvolution_BC_Atlas/PBCA_InstaPrism_reference.Rds"
GENESET_FILE  <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/ER_scoring_gene_sets.xlsx"
SCORES_FILE   <- file.path(BASE, "ER_signaling", "geneset_composite_scores_per_sample.csv")   # script 2
ESTIMATE_FILE <- file.path(BASE, "Estimate", "ESTIMATE_scores_TPM.csv")                       # script 3 (optional)
OUT_DIR       <- file.path(BASE, "Deconvolution")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_TOP              <- 100                        # genes per compartment in the main analysis
N_GRID             <- c(25, 50, 100, 200, 300)   # gene-number sensitivity check
MIN_TPM_MEAN       <- 5                          # gene must have mean TPM >= this in the 16 samples
MIN_VALUE_QUANTILE <- 0.75                       # gene must be in the top 25% of the compartment profile
PSEUDO             <- 1e-6                       # pseudocount for the log2 ratio

ER_FEATS <- c("ER_activity_index", "Resistance_index")   # scores adjusted for composition

# Genes removed from ALL lists so that composition scores are independent of the ER scores
EXCLUDE_SETS <- c("er_upreg", "er_pos_corr", "er_neg_corr", "TamR_up", "TamR_down")   # sets in the workbook
ER_PANEL     <- c("ESR1", "FOXA1", "GATA3", "TFAP2C", "GRHL2", "XBP1", "MYB", "SPDEF", "PBX1",
                  "NR2F2", "FOXM1", "E2F1", "MYC", "JUN", "FOS", "RELA", "NFKB1", "STAT3", "SOX9", "HIF1A")

# The 17 named (non-cancer) cell types of the PBCA reference; every other column is a cancer-cell profile
IMMUNE_TYPES     <- c("CD8+ T Cells", "CD4+ T Cells", "Regulatory T Cells", "NK Cells", "B Cells",
                      "Plasma Cells", "Macrophages", "Monocytes", "Dendritic Cells", "MDSCs",
                      "Mast Cells", "Neutrophils")
STROMAL_TYPES    <- c("Fibroblasts", "Perivascular-like (PVL) Cells", "Endothelial Cells")
NORMAL_EPI_TYPES <- c("Epithelial Cells", "Myoepithelial Cells")
NAMED_TYPES      <- c(IMMUNE_TYPES, STROMAL_TYPES, NORMAL_EPI_TYPES)

COMPARTMENTS <- c("Epithelial", "Stromal", "Immune")
SCORE_NAMES  <- c("PBCA_Epithelial", "PBCA_Stromal", "PBCA_Immune", "PBCA_TumorContent")

# ---- Helpers ---------------------------------------------------------------------
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
zs <- function(x) as.numeric(scale(x))

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

# Per-patient change (Recurrence - Primary) for a set of score columns
delta_table <- function(df, feats) {
  df %>%
    select(patient, timepoint, all_of(feats)) %>%
    pivot_longer(all_of(feats), names_to = "feature", values_to = "value") %>%
    pivot_wider(names_from = timepoint, values_from = value) %>%
    mutate(delta = Recurrence - Primary) %>%
    select(patient, feature, delta) %>%
    pivot_wider(names_from = feature, values_from = delta)
}

# delta_outcome ~ 1 + delta_covariate: the intercept is the outcome change expected if
# the composition score did not change (same as patient + timepoint + covariate)
adjust_fit <- function(dl, outcome, covariate) {
  fit <- lm(dl[[outcome]] ~ dl[[covariate]])
  co  <- summary(fit)$coefficients
  ci  <- confint(fit)
  tibble(outcome = outcome, covariate = covariate, n_pairs = nrow(dl),
         unadjusted_mean_delta = mean(dl[[outcome]]),
         adjusted_delta = co[1, 1], adjusted_se = co[1, 2],
         adjusted_ci_low = ci[1, 1], adjusted_ci_high = ci[1, 2],
         p_adjusted = co[1, 4], slope = co[2, 1], p_slope = co[2, 4],
         r_squared = summary(fit)$r.squared)
}

# ---- 1. Metadata (same 8 pairs as the other scripts) -----------------------------
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

# ---- 2. Bulk TPM (gene symbols) ---------------------------------------------------
tpm <- read.delim(TPM_FILE, row.names = 1, check.names = FALSE)
stopifnot(all(meta$sample %in% colnames(tpm)))
tpm_sym  <- to_symbols(tpm[, meta$sample])
mean_tpm <- rowMeans(tpm_sym)
logt     <- log2(tpm_sym + 1)

# ---- 3. Genes to remove (ER axis) --------------------------------------------------
gs_long <- read_excel(GENESET_FILE, sheet = "Genes (long)")
colnames(gs_long) <- c("set", "category", "gene")
hallmark <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                     error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
exclude_genes <- unique(c(
  gs_long$gene[gs_long$set %in% EXCLUDE_SETS],
  hallmark$gene_symbol[hallmark$gs_name %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY",
                                               "HALLMARK_ESTROGEN_RESPONSE_LATE")],
  ER_PANEL))
cat("ER-axis genes excluded from all lists:", length(exclude_genes), "\n")

# ---- 4. Compartment-specific gene lists from the PBCA reference ---------------------
ref <- readRDS(PBCA_REF_FILE)
phi <- as.matrix(ref@phi.cs)                      # genes x profiles (cell types + cancer-cell profiles)
message("Reference: ", nrow(phi), " genes x ", ncol(phi), " profiles")
if (mean(grepl("^ENSG", rownames(phi))) > 0.5)
  stop("Reference row names look like Ensembl IDs; gene symbols are needed.")
missing_named <- setdiff(NAMED_TYPES, colnames(phi))
if (length(missing_named) > 0) stop("Cell types not in the reference: ", paste(missing_named, collapse = ", "))

cancer_cols <- setdiff(colnames(phi), NAMED_TYPES)
message("Cancer-cell profiles used: ", length(cancer_cols))
if (length(cancer_cols) < 3) stop("Too few cancer-cell profiles found.")
cancer_mean <- rowMeans(phi[, cancer_cols, drop = FALSE])
cancer_max  <- apply(phi[, cancer_cols, drop = FALSE], 1, max)
rowmax      <- function(cols) apply(phi[, cols, drop = FALSE], 1, max)

# Value of each gene in the compartment, and its highest value outside the compartment
comp_value <- list(
  Epithelial = cancer_mean,                     # as in the pILC cancer-epithelial script
  Stromal    = rowmax(STROMAL_TYPES),
  Immune     = rowmax(IMMUNE_TYPES))
other_value <- list(
  Epithelial = rowmax(c(STROMAL_TYPES, IMMUNE_TYPES)),                                  # non-epithelial types
  Stromal    = pmax(rowmax(c(IMMUNE_TYPES, NORMAL_EPI_TYPES)), cancer_max),
  Immune     = pmax(rowmax(c(STROMAL_TYPES, NORMAL_EPI_TYPES)), cancer_max))

rank_compartment <- function(comp) {
  df <- tibble(gene = rownames(phi), compartment = comp,
               compartment_value = comp_value[[comp]], other_value = other_value[[comp]]) %>%
    mutate(log2_specificity = log2((compartment_value + PSEUDO) / (other_value + PSEUDO)),
           mean_TPM = unname(mean_tpm[gene]))
  cutoff <- quantile(df$compartment_value, MIN_VALUE_QUANTILE)
  df <- df %>%
    filter(!is.na(mean_TPM), mean_TPM >= MIN_TPM_MEAN, compartment_value >= cutoff) %>%
    arrange(-log2_specificity)
  n_removed <- sum(head(df$gene %in% exclude_genes, N_TOP))   # ER-axis genes among the top N
  df <- df %>% filter(!gene %in% exclude_genes) %>% mutate(rank = row_number())
  list(table = df, n_removed_top = n_removed)
}
ranked     <- setNames(lapply(COMPARTMENTS, rank_compartment), COMPARTMENTS)
gene_lists <- lapply(ranked, `[[`, "table")
for (k in COMPARTMENTS) {
  if (nrow(gene_lists[[k]]) < max(N_GRID)) stop("Fewer genes than max(N_GRID) pass the filters for ", k)
}

list_summary <- tibble(
  compartment = COMPARTMENTS,
  n_genes_passing_filters = sapply(gene_lists, nrow),
  n_ER_axis_genes_removed_from_top_N = sapply(ranked, `[[`, "n_removed_top"),
  top10_genes = sapply(gene_lists, function(d) paste(head(d$gene, 10), collapse = ", ")))
write.csv(list_summary, file.path(OUT_DIR, "PBCA_gene_list_summary.csv"), row.names = FALSE)
print(as.data.frame(list_summary))
write.csv(bind_rows(lapply(gene_lists, function(d) head(d, max(N_GRID)))),
          file.path(OUT_DIR, "PBCA_compartment_gene_lists.csv"), row.names = FALSE)

# ---- 5. Composite scores ------------------------------------------------------------
# Compartment score = mean z-score (across the 16 samples) of log2(TPM+1) of the top-N genes
comp_score <- function(genes) {
  g <- intersect(genes, rownames(logt))
  Z <- t(scale(t(logt[g, , drop = FALSE])))
  Z[is.na(Z)] <- 0
  colMeans(Z)
}
make_scores <- function(n) {
  s <- sapply(COMPARTMENTS, function(k) comp_score(head(gene_lists[[k]]$gene, n)))
  s <- as.data.frame(s)
  colnames(s) <- c("PBCA_Epithelial", "PBCA_Stromal", "PBCA_Immune")
  # Tumor content index: high = epithelial-rich, low = stromal / immune-rich
  s$PBCA_TumorContent <- zs(s$PBCA_Epithelial) - (zs(s$PBCA_Stromal) + zs(s$PBCA_Immune)) / 2
  rownames(s) <- colnames(logt)
  rownames_to_column(s, "Sample")
}

if (!file.exists(SCORES_FILE)) stop("Run script 2 first; not found: ", SCORES_FILE)
er_scores <- read.csv(SCORES_FILE, stringsAsFactors = FALSE) %>% select(sample, all_of(ER_FEATS))

attach_info <- function(s) {
  s %>%
    left_join(meta %>% select(sample, block, patient, timepoint), by = c("Sample" = "sample")) %>%
    left_join(er_scores, by = c("Sample" = "sample")) %>%
    arrange(patient, timepoint)
}
scores <- attach_info(make_scores(N_TOP))
stopifnot(!anyNA(scores$ER_activity_index), !anyNA(scores$patient))
write.csv(scores, file.path(OUT_DIR, "PBCA_composition_scores_per_sample.csv"), row.names = FALSE)

# ---- 6. Paired tests ----------------------------------------------------------------
long_sc <- scores %>%
  select(patient, timepoint, all_of(SCORE_NAMES)) %>%
  pivot_longer(all_of(SCORE_NAMES), names_to = "feature", values_to = "value") %>%
  mutate(feature = factor(feature, levels = SCORE_NAMES))

sc_stats <- long_sc %>%
  group_by(feature) %>% group_modify(~ paired_test(.x)) %>% ungroup() %>%
  mutate(FDR_wilcoxon = p.adjust(p_wilcoxon, "BH"), FDR_ttest = p.adjust(p_ttest, "BH"))
write.csv(sc_stats, file.path(OUT_DIR, "PBCA_composition_paired_stats.csv"), row.names = FALSE)
print(as.data.frame(sc_stats))

p_pair <- ggplot(long_sc, aes(timepoint, value, group = patient, colour = patient)) +
  geom_line(alpha = 0.8) + geom_point(size = 2) +
  facet_wrap(~feature, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Composite score", colour = "Patient") +
  theme_bw(base_size = 11)
save_fig(p_pair, "PBCA_composition_scores_paired", w = 12, h = 4)

# ---- 7. ER-score change adjusted for composition ------------------------------------
# Primary covariate (pre-specified): PBCA_TumorContent. The three compartments are sensitivity checks.
dl <- delta_table(scores, c(SCORE_NAMES, ER_FEATS))

delta_cor <- bind_rows(lapply(ER_FEATS, function(o) bind_rows(lapply(SCORE_NAMES, function(x) {
  ct <- suppressWarnings(cor.test(dl[[o]], dl[[x]], method = "spearman", exact = FALSE))
  tibble(outcome = o, covariate = x, spearman_rho = unname(ct$estimate), p = ct$p.value)
}))))
write.csv(delta_cor, file.path(OUT_DIR, "PBCA_delta_correlations.csv"), row.names = FALSE)

adj <- bind_rows(lapply(ER_FEATS, function(o)
  bind_rows(lapply(SCORE_NAMES, function(x) adjust_fit(dl, o, x)))))
write.csv(adj, file.path(OUT_DIR, "ER_change_adjusted_for_PBCA_composition.csv"), row.names = FALSE)
print(as.data.frame(adj))

sc_long <- dl %>%
  select(patient, ER_activity_index, all_of(SCORE_NAMES)) %>%
  pivot_longer(all_of(SCORE_NAMES), names_to = "covariate", values_to = "delta_x") %>%
  mutate(covariate = factor(covariate, levels = SCORE_NAMES))
p_sc <- ggplot(sc_long, aes(delta_x, ER_activity_index)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40", formula = y ~ x) +
  geom_point(size = 2.5) +
  geom_text(aes(label = patient), size = 3, vjust = -0.8, check_overlap = TRUE) +
  facet_wrap(~covariate, scales = "free_x", nrow = 1) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(x = "Change in composite score (Recurrence - Primary)",
       y = "Change in ER activity index (Recurrence - Primary)") +
  theme_bw(base_size = 11)
save_fig(p_sc, "PBCA_composition_vs_ER_index_change", w = 14, h = 4.5)

# ---- 8. Does the gene-list length change the answer? --------------------------------
sens <- bind_rows(lapply(N_GRID, function(n) {
  sn  <- attach_info(make_scores(n))
  dln <- delta_table(sn, c(SCORE_NAMES, ER_FEATS))
  bind_rows(lapply(SCORE_NAMES, function(x) adjust_fit(dln, "ER_activity_index", x))) %>%
    mutate(n_genes_per_compartment = n, .before = 1)
}))
write.csv(sens, file.path(OUT_DIR, "ER_change_adjusted_by_gene_number.csv"), row.names = FALSE)
print(as.data.frame(sens %>% select(n_genes_per_compartment, covariate, unadjusted_mean_delta,
                                    adjusted_delta, adjusted_ci_low, adjusted_ci_high, p_adjusted)))

# ---- 9. Check against the ESTIMATE scores of script 3 (optional) ---------------------
if (file.exists(ESTIMATE_FILE)) {
  es  <- read.csv(ESTIMATE_FILE, stringsAsFactors = FALSE)
  chk <- scores %>% left_join(es %>% select(Sample, StromalScore, ImmuneScore, TumorPurity), by = "Sample")
  pairs <- list(c("PBCA_Epithelial", "TumorPurity"), c("PBCA_Stromal", "StromalScore"),
                c("PBCA_Immune", "ImmuneScore"), c("PBCA_TumorContent", "TumorPurity"))
  est_cor <- bind_rows(lapply(pairs, function(pr) {
    a  <- chk[[pr[1]]]
    b  <- chk[[pr[2]]]
    ct <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE))
    tibble(PBCA_score = pr[1], ESTIMATE_score = pr[2], spearman_rho = unname(ct$estimate),
           p = ct$p.value, n_samples = sum(is.finite(a) & is.finite(b)))
  }))
  write.csv(est_cor, file.path(OUT_DIR, "PBCA_vs_ESTIMATE_correlation.csv"), row.names = FALSE)
  print(as.data.frame(est_cor))
} else message("ESTIMATE file not found; check against ESTIMATE skipped.")

sessionInfo()