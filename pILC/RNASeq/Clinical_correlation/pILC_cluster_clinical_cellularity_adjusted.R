# pILC: association of RNA clusters (NMF k=3) with clinical features,
# BEFORE and AFTER adjustment for tumor cellularity.
# Outputs: one table (unadjusted vs adjusted results) + 3 figures (PNG + PDF).

suppressPackageStartupMessages(library(ggplot2))
SEED <- 42
set.seed(SEED)

# ---------------------------------------------------------------- 1. Paths
root   <- Sys.getenv("PILC_ROOT", "/ix1/alee/LO_LAB/Personal/Rahul/pILC")
clus_f <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
clin_f <- file.path(root, "Clinical_correlation/pILC_ClinicoPathological_details.csv")
out_dir <- file.path(root, "Clinical_correlation/Cluster_vs_clinical_cellularity_adjusted")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")

# ---------------------------------------------------------------- 2. Load and merge
norm_id <- function(x) {                       # TP09.9 / TP09-9 -> TP09-09
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}
clus <- read.csv(clus_f); clus$Key <- norm_id(clus$sample)
clin <- read.csv(clin_f, check.names = FALSE, strip.white = TRUE); clin$Key <- norm_id(clin$Sample)
miss <- setdiff(clus$Key, clin$Key)
if (length(miss)) stop("No clinical data for: ", paste(miss, collapse = ", "))
d <- merge(clus[, c("Key", "cluster")], clin, by = "Key")
d$Cluster <- factor(d$cluster, levels = 1:3)

# ---------------------------------------------------------------- 3. Clean variables
# "<5" style values are censored: use half the limit (2.5), a simple and stated choice
pct <- function(x) {
  x <- trimws(as.character(x))
  out <- suppressWarnings(as.numeric(sub("^<", "", x)))
  out[grepl("^<", x) & !is.na(out)] <- out[grepl("^<", x) & !is.na(out)] / 2
  out
}
tw  <- function(col) trimws(as.character(d[[col]]))
d$cell      <- pct(d[["tumor_cellularity (%)"]])       # covariate to adjust for
d$cell10    <- d$cell / 10                              # effect per 10% cellularity
d$Age       <- suppressWarnings(as.numeric(d[["Age at Diagnosis"]]))
d$TILs      <- pct(d[["TILs (%)"]])
d$Pleo_cell <- pct(d[["pleomorphic_invasive_compoent_cellularity (%)"]])
d$Grade_III <- ifelse(tw("Grade/Differentiation-Desc") == "Grade III", 1,
                      ifelse(tw("Grade/Differentiation-Desc") == "Grade II", 0, NA))
d$ER_pos    <- ifelse(tw("Tumor Marker #1 (ER)") == "Positive", 1,
                      ifelse(tw("Tumor Marker #1 (ER)") == "Negative", 0, NA))
d$PR_pos    <- ifelse(tw("Tumor Marker #2 (PR)") == "Positive", 1,
                      ifelse(tw("Tumor Marker #2 (PR)") == "Negative", 0, NA))
d$HER2_pos  <- ifelse(tw("BREAST: Her-2/Neu-Desc") == "Positive", 1,
                      ifelse(tw("BREAST: Her-2/Neu-Desc") == "Negative", 0, NA))
d$HRpos_HER2neg <- ifelse(tw("HR+ Her2-") == "Yes", 1, 0)

cont_vars <- c(TILs = "TILs (%)", Pleo_cell = "Pleomorphic component cellularity (%)")
bin_vars  <- c(Grade_III = "Grade III (vs II)", ER_pos = "ER positive", PR_pos = "PR positive",
               HER2_pos = "HER2 positive", HRpos_HER2neg = "HR+/HER2- subtype")

# ---------------------------------------------------------------- 4. Confounder check
kw_cell <- kruskal.test(cell ~ Cluster, d)$p.value
message("Tumor cellularity by cluster: Kruskal-Wallis p = ", signif(kw_cell, 3))
print(tapply(d$cell, d$Cluster, median, na.rm = TRUE))

# ---------------------------------------------------------------- 5. Tests
# Continuous: rank-based (robust for small n). Unadjusted = Kruskal-Wallis / lm(rank ~ cluster);
# adjusted = lm(rank ~ cellularity + cluster) vs lm(rank ~ cellularity) (F test for cluster).
summ_cont <- function(x, cl) paste(sapply(levels(cl), function(k) {
  v <- x[cl == k & !is.na(x)]
  sprintf("c%s %.0f [%.0f-%.0f] (n=%d)", k, median(v), quantile(v, .25), quantile(v, .75), length(v))
}), collapse = "; ")
test_cont <- function(v, label) {
  df <- d[!is.na(d[[v]]) & !is.na(d$cell), c("Cluster", "cell10", v)]; names(df)[3] <- "y"
  df$r <- rank(df$y)
  m_un0 <- lm(r ~ 1, df);     m_un1 <- lm(r ~ Cluster, df)
  m_ad0 <- lm(r ~ cell10, df); m_ad1 <- lm(r ~ cell10 + Cluster, df)
  data.frame(Feature = label, Type = "continuous", N = nrow(df),
             p_unadjusted = kruskal.test(y ~ Cluster, df)$p.value,
             p_adjusted   = anova(m_ad0, m_ad1)[2, "Pr(>F)"],
             effect_unadjusted = 1 - deviance(m_un1) / deviance(m_un0),   # R2 of cluster on ranks
             effect_adjusted   = 1 - deviance(m_ad1) / deviance(m_ad0),   # partial R2 of cluster
             by_cluster = summ_cont(df$y, df$Cluster), separation = FALSE)
}
# Binary: Fisher exact (unadjusted) and logistic LRT for cluster, with/without cellularity.
test_bin <- function(v, label) {
  df <- d[!is.na(d[[v]]) & !is.na(d$cell), c("Cluster", "cell10", v)]; names(df)[3] <- "y"
  tab <- table(df$Cluster, df$y)
  f0 <- suppressWarnings(glm(y ~ cell10, binomial, df))
  f1 <- suppressWarnings(glm(y ~ cell10 + Cluster, binomial, df))
  sep <- !f1$converged || any(fitted(f1) < 1e-8 | fitted(f1) > 1 - 1e-8)
  by_cl <- paste(sapply(levels(df$Cluster), function(k)
    sprintf("c%s %d/%d", k, sum(df$y[df$Cluster == k]), sum(df$Cluster == k))), collapse = "; ")
  data.frame(Feature = label, Type = "binary", N = nrow(df),
             p_unadjusted = tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_),
             p_adjusted   = suppressWarnings(anova(f0, f1, test = "Chisq")[2, "Pr(>Chi)"]),
             effect_unadjusted = NA_real_, effect_adjusted = NA_real_,
             by_cluster = by_cl, separation = sep)
}
res <- rbind(
  do.call(rbind, Map(test_cont, names(cont_vars), cont_vars)),
  do.call(rbind, Map(test_bin,  names(bin_vars),  bin_vars)))
res$FDR_unadjusted <- p.adjust(res$p_unadjusted, "BH")
res$FDR_adjusted   <- p.adjust(res$p_adjusted,   "BH")
rownames(res) <- NULL
# Tumor cellularity itself (the confounder): unadjusted only, kept out of the FDR family
res <- rbind(res, data.frame(
  Feature = "Tumor cellularity (%) [confounder; not in FDR]", Type = "continuous",
  N = sum(!is.na(d$cell)), p_unadjusted = kw_cell, p_adjusted = NA_real_,
  effect_unadjusted = NA_real_, effect_adjusted = NA_real_,
  by_cluster = summ_cont(d$cell, d$Cluster), separation = FALSE,
  FDR_unadjusted = NA_real_, FDR_adjusted = NA_real_))
write.csv(res, file.path(out_dir, "Cluster_vs_clinical_unadjusted_vs_cellularity_adjusted.csv"), row.names = FALSE)

# ---- Pairwise cluster comparisons (post hoc) for the continuous features
# Unadjusted: Wilcoxon rank-sum. Adjusted: F test of cluster in lm(rank ~ cellularity + cluster).
# BH correction across the 3 cluster pairs within each feature.
pair_test <- function(v, label, adjust = TRUE) {
  prs <- combn(levels(d$Cluster), 2, simplify = FALSE)
  out <- do.call(rbind, lapply(prs, function(pr) {
    df <- d[d$Cluster %in% pr & !is.na(d[[v]]) & !is.na(d$cell), c("Cluster", "cell10", v)]
    names(df)[3] <- "y"; df$Cluster <- droplevels(df$Cluster); df$r <- rank(df$y)
    med <- tapply(df$y, df$Cluster, median)
    data.frame(Feature = label, Comparison = paste0("cluster ", pr[1], " vs ", pr[2]), N = nrow(df),
               median_first = med[[1]], median_second = med[[2]],
               p_unadjusted = wilcox.test(y ~ Cluster, df, exact = FALSE)$p.value,
               p_adjusted = if (adjust) anova(lm(r ~ cell10, df), lm(r ~ cell10 + Cluster, df))[2, "Pr(>F)"] else NA_real_)
  }))
  out$FDR_unadjusted <- p.adjust(out$p_unadjusted, "BH")
  out$FDR_adjusted   <- p.adjust(out$p_adjusted,   "BH")
  out
}
pair_res <- do.call(rbind, c(
  list(pair_test("cell", "Tumor cellularity (%)", adjust = FALSE)),
  Map(function(v, lab) pair_test(v, lab), names(cont_vars), cont_vars)))
rownames(pair_res) <- NULL
write.csv(pair_res, file.path(out_dir, "Cluster_pairwise_comparisons_unadjusted_vs_cellularity_adjusted.csv"), row.names = FALSE)
print(pair_res)
print(res[, c("Feature", "N", "p_unadjusted", "p_adjusted", "FDR_unadjusted", "FDR_adjusted", "separation")])

# ---------------------------------------------------------------- 6. Figures
th <- theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))
save_both <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
}

# (a) Continuous features by cluster: violin + box + points, with pairwise brackets.
# Brackets use pairwise Wilcoxon (unadjusted figure) or the cellularity-adjusted pairwise tests
# (adjusted figure), BH-corrected across the 3 cluster pairs within each feature.
LABEL_TYPE <- "stars"          # "stars" (ns, *, **, ***) or "p" to print numbers
lab_fun <- function(p) {
  if (LABEL_TYPE == "p") {
    ifelse(p < 1e-4, "p < 0.0001", sprintf("p = %.3f", p))
  } else {
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  }
}
fmt_p  <- function(p) ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))
short  <- function(x) sub(" \\(%\\)$", "", x)

draw_violin <- function(long, pair_df, p_col, facet_lab, title, ylab, caption) {
  long$Facet <- factor(facet_lab[as.character(long$Feature)], levels = unname(facet_lab))
  # brackets: 1 vs 2 lowest, 1 vs 3 middle, 2 vs 3 top
  br <- do.call(rbind, lapply(split(pair_df, pair_df$Feature), function(pr) {
    f   <- pr$Feature[1]
    v   <- long$Value[long$Feature == f]
    rng <- diff(range(v)); step <- 0.12 * rng
    pts <- strsplit(sub("cluster ", "", pr$Comparison), " vs ")
    data.frame(Feature = f, x1 = as.numeric(sapply(pts, `[`, 1)), x2 = as.numeric(sapply(pts, `[`, 2)),
               y = max(v) + step * seq_len(nrow(pr)), tick = 0.025 * rng,
               label = lab_fun(pr[[p_col]]), stringsAsFactors = FALSE)
  }))
  br$Facet <- factor(facet_lab[br$Feature], levels = unname(facet_lab))
  ggplot(long, aes(Cluster, Value)) +
    geom_violin(aes(fill = Cluster), color = "black", scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
    geom_jitter(width = 0.06, size = 1.8, alpha = 0.8) +
    geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = br, aes(x = x1, xend = x1, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_segment(data = br, aes(x = x2, xend = x2, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_text(data = br, aes(x = (x1 + x2) / 2, y = y, label = label), inherit.aes = FALSE,
              vjust = -0.2, size = 4.5) +
    facet_wrap(~Facet, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = cluster_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +   # room for the top bracket label
    labs(title = title, y = ylab, caption = paste(strwrap(caption, 90), collapse = "\n")) +   # wrap long caption
    th + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5),
               plot.caption = element_text(hjust = 0, size = 10), strip.text = element_text(face = "bold"),
               panel.border = element_rect(color = "grey30", fill = NA), axis.line = element_blank())
}

feat_levels <- c("Tumor cellularity (%)", unname(cont_vars))

# (a1) Before adjustment: raw values
long <- do.call(rbind, c(
  list(data.frame(Cluster = d$Cluster, Feature = "Tumor cellularity (%)", Value = d$cell)),
  Map(function(v, lab) data.frame(Cluster = d$Cluster, Feature = lab, Value = d[[v]]),
      names(cont_vars), cont_vars)))
long <- long[!is.na(long$Value), ]
kw_tab <- c(kw_cell, res$p_unadjusted[match(unname(cont_vars), res$Feature)])
facet_un <- setNames(paste0(short(feat_levels), "\nKruskal-Wallis p = ", fmt_p(kw_tab)), feat_levels)
p_a <- draw_violin(long, pair_res, "FDR_unadjusted", facet_un,
                   "Comparison of Clinical Features Across Clusters", "Value (in %)",
                   "Pairwise Wilcoxon rank-sum tests, BH-adjusted within each feature; ns = not significant")
save_both(p_a, "Clinical_continuous_by_cluster_unadjusted", 12, 5)

# (a2) After adjustment for tumor cellularity (tumor cellularity itself cannot be adjusted)
# Value = rank percentile (0-100) minus the tumor-cellularity effect, taken from the model
# percentile ~ cellularity + cluster, so the cluster differences shown match the adjusted tests.
adj_long <- do.call(rbind, Map(function(v, lab) {
  df <- d[!is.na(d[[v]]) & !is.na(d$cell), c("Cluster", "cell10", v)]; names(df)[3] <- "y"
  df$pct <- 100 * (rank(df$y) - 0.5) / nrow(df)
  b <- coef(lm(pct ~ cell10 + Cluster, df))[["cell10"]]
  data.frame(Cluster = df$Cluster, Feature = lab, Value = df$pct - b * (df$cell10 - mean(df$cell10)))
}, names(cont_vars), cont_vars))
adj_p   <- res$p_adjusted[match(unname(cont_vars), res$Feature)]
facet_ad <- setNames(paste0(short(unname(cont_vars)), "\nAdjusted p = ", fmt_p(adj_p)), unname(cont_vars))
p_adj <- draw_violin(adj_long, pair_res[pair_res$Feature %in% cont_vars, ], "FDR_adjusted", facet_ad,
                     "Clinical features across clusters, adjusted for tumor cellularity",
                     "Cellularity-adjusted rank percentile",
                     "Rank percentile minus the tumor-cellularity effect; brackets: cellularity-adjusted pairwise tests, BH-adjusted within each feature")
save_both(p_adj, "Clinical_continuous_by_cluster_cellularity_adjusted", 9.5, 5.5)

# (b) Binary features: proportion by cluster
bl <- do.call(rbind, Map(function(v, lab) {
  df <- d[!is.na(d[[v]]), ]
  agg <- aggregate(df[[v]], list(Cluster = df$Cluster), mean)
  data.frame(Cluster = agg$Cluster, Feature = ifelse(lab == "Grade III (vs II)", "Grade III", lab), Proportion = agg$x)  # panel title only; the results table keeps the full name
}, names(bin_vars), bin_vars))
p_b <- ggplot(bl, aes(Cluster, Proportion, fill = Cluster)) +
  geom_col(width = 0.7) + facet_wrap(~Feature, nrow = 1) +
  scale_fill_manual(values = cluster_cols) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Binary clinical features by RNA cluster", y = "Proportion of tumors") +
  th + theme(legend.position = "none")
save_both(p_b, "Clinical_binary_by_cluster", 13, 4.5)

# (c) Before vs after adjustment: -log10(p) per feature
rp <- res[!grepl("confounder", res$Feature), ]   # cellularity row has no adjusted result
pl <- rbind(data.frame(Feature = rp$Feature, Model = "Unadjusted", p = rp$p_unadjusted),
            data.frame(Feature = rp$Feature, Model = "Adjusted for cellularity", p = rp$p_adjusted))
pl$Feature <- factor(pl$Feature, levels = rev(rp$Feature))
p_c <- ggplot(pl, aes(-log10(p), Feature, color = Model)) +
  geom_line(aes(group = Feature), color = "grey60") +
  geom_point(size = 3.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("Unadjusted" = "#4E79A7", "Adjusted for cellularity" = "#E15759")) +
  labs(title = "Cluster association with clinical features", x = expression(-log[10](p)), y = NULL) + th
save_both(p_c, "Clinical_p_unadjusted_vs_adjusted", 9, 5)

message("Done. Output in ", out_dir)