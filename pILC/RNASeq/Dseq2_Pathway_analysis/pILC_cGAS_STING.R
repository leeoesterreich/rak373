###############################################################################
# pILC RNA-seq | script: cGAS-STING signature across the three NMF clusters (focus: Cluster 2),
#                           including the check that the Cluster 2 difference is more than a TILs effect.
# Reads the outputs of pILC_01a_DESeq2_clusters.R (run it first). This single script replaces
# pILC_07_cGAS_STING.R and pILC_07b_cGAS_STING_TILs_adjusted.R.
#
# Signature : the 7-gene cGAS-STING panel of Chen et al., npj Breast Cancer 2024 (C6orf150, CCL5, CXCL10,
#             IRF3, TBK1, TMEM173, STAT1). Current HGNC symbols are used: C6orf150 = CGAS, TMEM173 = STING1.
#             Score = mean of the per-gene z-scores (VST, cluster-blind). The paper's exact formula was not
#             available to me, so this follows the same idea and may not be identical.
# Extra     : the panel is also split into a sensing/signalling core (CGAS, STING1, TBK1, IRF3) and a
#             downstream output (CCL5, CXCL10, STAT1); the output genes are interferon / T-cell genes that
#             rise with immune infiltration, the core genes much less.
# Clinical  : tumor cellularity, TILs and pleomorphic component are read here directly from the clinical file
#             with a parser that handles "5", "5%", "5-10", "<5", ">50" and reports every conversion
#             (the as.numeric() in script 01a turns any non-numeric entry into NA without a warning).
# Tests     : Kruskal-Wallis across clusters, pairwise Wilcoxon, Cluster 2 vs the rest (Wilcoxon, all samples),
#             Cluster 2 vs the rest adjusted for TILs (linear models, three versions), Grade 2 vs 3,
#             Spearman with the clinical variables.
# TILs model: value ~ Cluster2 (unadjusted) and value ~ Cluster2 + TILs (adjusted), on the samples with a TILs
#             value. "Cluster 2 effect" = difference in mean value, Cluster 2 minus the other samples. If it stays
#             similar and significant after adding TILs, TILs do not explain the difference. Versions: linear TILs;
#             log2(TILs + 1) (reduces the pull of very high TILs); linear without the highest-TILs sample.
# Reference : Chen M, Yu S, van der Sluis T, Zwager MC, Schroder CP, van der Vegt B, et al. (van Vugt MATM, corresponding).
#             cGAS-STING pathway expression correlates with genomic instability and immune cell infiltration in
#             breast cancer. npj Breast Cancer. 2024;10:1. doi:10.1038/s41523-023-00609-z
#             (the author list is shortened: please complete it from the journal page before citing)
# Note      : clusters were derived from the same expression data, so p-values describe how the groups differ and
#             are not independent validation. With ~31 samples the results are exploratory. TILs adjusts for one
#             variable only; other cell types (e.g. plasma cells) can still contribute.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
PILC_ROOT     <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"
ROOT          <- file.path(PILC_ROOT, "RNASeq")
SHARED        <- file.path(ROOT, "DESeq2_stranded", "Shared_inputs_for_pathway_analysis")
META_FILE     <- file.path(SHARED, "Sample_metadata.tsv")
VST_FILE      <- file.path(SHARED, "VST_blind_ssGSEA_input.tsv")        # gene symbols x samples (cluster-blind VST)
CLINICAL_FILE <- file.path(PILC_ROOT, "Clinical_correlation", "pILC_ClinicoPathological_details.csv")

OUT_ROOT <- file.path(ROOT, "Pathway_analysis", "cGAS_STING")            # no dates in names
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

FOCUS   <- "Cluster 2"      # cluster compared against the rest
FDR_CUT <- 0.05
LT_RULE <- "half"           # how an entry like "<5" is converted: "half" -> 2.5, "value" -> 5

# gene -> accepted symbols (first one found in the matrix is used)
PANEL <- list(CGAS    = c("CGAS", "MB21D1", "C6orf150"),
              STING1  = c("STING1", "TMEM173"),
              TBK1    = "TBK1",
              IRF3    = "IRF3",
              STAT1   = "STAT1",
              CCL5    = "CCL5",
              CXCL10  = "CXCL10")
CORE   <- c("CGAS", "STING1", "TBK1", "IRF3")
OUTPUT <- c("CCL5", "CXCL10", "STAT1")
REFERENCE <- c("Signature source: Chen M, Yu S, van der Sluis T, et al. cGAS-STING pathway expression correlates with genomic instability and immune cell infiltration in breast cancer. npj Breast Cancer. 2024;10:1. doi:10.1038/s41523-023-00609-z",
               "Panel genes (as published): C6orf150, CCL5, CXCL10, IRF3, TBK1, TMEM173, STAT1; current HGNC symbols: C6orf150 = CGAS, TMEM173 = STING1.",
               "Score used here: mean of per-gene z-scores of the cluster-blind VST matrix (the paper's exact formula was not checked).")

# clinical variables read from the clinical file (name = output column, value = column-name pattern)
CLIN_PATTERNS <- c(Tumor_cellularity = "^tumor_cellularity", TILs = "^TILs", Pleomorphic_component = "pleomorph")

for (f in c(META_FILE, VST_FILE, CLINICAL_FILE)) if (!file.exists(f)) stop("Missing input (run pILC_01a first): ", f)

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
grade_colors   <- c("Grade 2" = "#66C2A5", "Grade 3" = "#E78AC3")
GRADE_PAIR     <- names(grade_colors)
clin_vars      <- names(CLIN_PATTERNS)

# ------------------------------- 2. HELPERS ----------------------------------
norm_id <- function(x) {                       # same sample-ID style as script 01a
  x <- gsub("-", ".", sub("_S[0-9]+$", "", trimws(as.character(x))))
  sub("\\.0+([0-9])", ".\\1", x)
}
save_plot <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
}
save_pheatmap <- function(draw_fun, stem, w, h) {
  ph <- draw_fun()                                   # pheatmap(..., silent = TRUE) returns the figure without drawing it
  png(paste0(stem, ".png"), width = w * 300, height = h * 300, res = 300)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = w, height = h)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
  sz <- file.size(paste0(stem, c(".png", ".pdf")))
  if (any(is.na(sz) | sz < 10000)) warning("Heatmap file is missing or very small: ", stem)
}
fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), formatC(x, format = "f", digits = 3))
pre <- function(x) file.path(OUT_ROOT, paste0("cGAS_STING_", x))

# One clinical entry -> number. Returns the number and the rule used, so every conversion can be checked.
parse_one <- function(raw) {
  s <- tolower(trimws(as.character(raw)))
  s <- gsub("%", "", s)
  s <- gsub("[–—]", "-", s)          # en dash / em dash -> hyphen
  s <- gsub("\\s+", "", s)
  num <- "[0-9]+(\\.[0-9]+)?"
  if (is.na(raw) || s %in% c("", "na", "n/a", "nd", "unknown", "-", "."))
    return(list(value = NA_real_, rule = "missing"))
  if (grepl(paste0("^", num, "$"), s))
    return(list(value = as.numeric(s), rule = "number"))
  if (grepl(paste0("^", num, "(-|to)", num, "$"), s)) {
    p <- as.numeric(strsplit(s, "-|to")[[1]])
    return(list(value = mean(p), rule = "range: midpoint"))
  }
  if (grepl(paste0("^<=?", num, "$"), s)) {
    v <- as.numeric(sub("^<=?", "", s))
    return(list(value = if (LT_RULE == "half") v / 2 else v, rule = paste0("less than: ", LT_RULE)))
  }
  if (grepl(paste0("^>=?", num, "$"), s))
    return(list(value = as.numeric(sub("^>=?", "", s)), rule = "greater than: limit value"))
  list(value = NA_real_, rule = "NOT PARSED")
}

# ------------------------------- 3. LOAD DATA --------------------------------
meta_dt <- fread(META_FILE)
need <- c("Sample", "Cluster", "Grade")
if (length(setdiff(need, names(meta_dt)))) stop("Sample_metadata.tsv lacks: ", paste(setdiff(need, names(meta_dt)), collapse = ", "))
meta_dt <- meta_dt[, need, with = FALSE]                        # clinical variables come from the clinical file below
meta_dt[, Sample_key := norm_id(Sample)]
meta_dt[!(Grade %in% GRADE_PAIR), Grade := NA_character_]        # blank / other grades -> NA
if (!(FOCUS %in% cluster_levels)) stop("FOCUS must be one of: ", paste(cluster_levels, collapse = ", "))

# 3a. Clinical variables from the clinical file (robust numeric parsing, every conversion is listed)
clin <- fread(CLINICAL_FILE)
if (!("Sample" %in% names(clin))) stop("The clinical file needs a column named 'Sample'.")
clin[, Sample_key := norm_id(Sample)]
if (anyDuplicated(clin$Sample_key)) stop("Duplicated samples in the clinical file.")
parse_tabs <- lapply(clin_vars, function(v) {
  col <- grep(CLIN_PATTERNS[[v]], names(clin), ignore.case = TRUE, value = TRUE)[1]
  message("Clinical column for ", v, ": ", ifelse(is.na(col), "NOT FOUND (variable skipped)", col))
  if (is.na(col)) return(NULL)
  pr <- lapply(clin[[col]], parse_one)
  data.table(Variable = v, Sample_key = clin$Sample_key, Raw = as.character(clin[[col]]),
             Value = vapply(pr, `[[`, numeric(1), "value"), Rule = vapply(pr, `[[`, character(1), "rule"))
})
parse_all <- rbindlist(parse_tabs)
clin_vars <- intersect(clin_vars, unique(parse_all$Variable))
fwrite(parse_all, pre("clinical_parsing_check.tsv"), sep = "\t", quote = FALSE, na = "")
message("How the clinical entries were converted (all entries: clinical_parsing_check.tsv):")
print(dcast(parse_all[, .N, by = .(Variable, Rule)], Rule ~ Variable, value.var = "N", fill = 0))
if (any(parse_all$Rule == "NOT PARSED")) {
  message("These entries could NOT be converted; check them:")
  print(parse_all[Rule == "NOT PARSED", .(Variable, Sample_key, Raw)])
}
clin_wide <- dcast(parse_all, Sample_key ~ Variable, value.var = "Value")
meta_dt <- merge(meta_dt, clin_wide, by = "Sample_key", all.x = TRUE, sort = FALSE)
if (anyNA(match(meta_dt$Sample_key, clin$Sample_key)))
  message("No clinical row for: ", paste(meta_dt$Sample[is.na(match(meta_dt$Sample_key, clin$Sample_key))], collapse = ", "))

metadata <- as.data.frame(meta_dt)
rownames(metadata) <- metadata$Sample
metadata$Cluster <- factor(metadata$Cluster, levels = cluster_levels)
message("Cluster sizes:"); print(table(metadata$Cluster))
message("Samples with a numeric value per clinical variable:")
print(colSums(!is.na(metadata[, clin_vars, drop = FALSE])))
HAS_TILS <- "TILs" %in% clin_vars && sum(is.finite(metadata$TILs)) >= 10
if (HAS_TILS && sum(is.finite(metadata$TILs)) < 25)
  warning("Fewer than 25 samples have a TILs value: check clinical_parsing_check.tsv before trusting the TILs models.")

vst_dt <- fread(VST_FILE)
expr <- as.matrix(vst_dt[, names(vst_dt)[-1], with = FALSE])
rownames(expr) <- vst_dt[[1]]
if (!all(rownames(metadata) %in% colnames(expr))) stop("Samples of Sample_metadata.tsv are missing in the VST matrix.")
expr <- expr[, rownames(metadata)]

# ------------------------------- 4. SIGNATURE --------------------------------
# 4a. Which symbol was found for each gene (check this table: a missing gene lowers the score's meaning)
gene_map <- rbindlist(lapply(names(PANEL), function(g) {
  hit <- PANEL[[g]][PANEL[[g]] %in% rownames(expr)][1]
  data.table(Gene = g, Accepted_symbols = paste(PANEL[[g]], collapse = " / "),
             Symbol_used = ifelse(is.na(hit), NA_character_, hit), Found = !is.na(hit))
}))
fwrite(gene_map, pre("gene_map.tsv"), sep = "\t", quote = FALSE, na = "")
print(gene_map)
writeLines(REFERENCE, pre("signature_reference.txt"))
if (sum(gene_map$Found) < 5) stop("Fewer than 5 of the 7 panel genes were found: check gene_map.tsv.")
if (any(!gene_map$Found)) warning("Not found in the matrix: ", paste(gene_map$Gene[!gene_map$Found], collapse = ", "))
found <- gene_map[Found == TRUE]

gene_expr <- expr[found$Symbol_used, , drop = FALSE]
rownames(gene_expr) <- found$Gene

# 4b. Scores = mean of per-gene z-scores across the samples
zg <- t(scale(t(gene_expr)))
mean_z <- function(genes) {
  g <- intersect(genes, rownames(zg))
  if (length(g) == 0) rep(NA_real_, ncol(zg)) else colMeans(zg[g, , drop = FALSE])
}
scores <- rbind(cGAS_STING_score        = mean_z(rownames(zg)),
                Sensing_core_score      = mean_z(CORE),
                Downstream_output_score = mean_z(OUTPUT))
scores <- scores[!apply(scores, 1, function(x) all(is.na(x))), , drop = FALSE]
feat <- rbind(scores, gene_expr)                       # scores (z-score units) + genes (VST units)
feature_names <- rownames(feat)
main_scores   <- intersect(c("cGAS_STING_score", "Sensing_core_score", "Downstream_output_score"), feature_names)
main_score    <- "cGAS_STING_score"

sample_tab <- data.table(Sample = rownames(metadata), Cluster = metadata$Cluster, Grade = metadata$Grade,
                         metadata[, clin_vars, drop = FALSE])
sample_tab <- merge(sample_tab, data.table(Sample = colnames(feat), t(feat)), by = "Sample", sort = FALSE)
sample_tab[, Focus := Cluster == FOCUS]
fwrite(sample_tab[, !"Focus"], pre("scores_and_gene_expression_per_sample.tsv"), sep = "\t", quote = FALSE, na = "")

fl <- melt(as.data.table(t(feat), keep.rownames = "Sample"), id.vars = "Sample",
           variable.name = "Feature", value.name = "Value")
fl[, Feature := as.character(Feature)]
fl[, Cluster := metadata[Sample, "Cluster"]]
fl[, Grade   := metadata[Sample, "Grade"]]
fl[, Focus   := Cluster == FOCUS]

# ------------------------------- 5. STATISTICS -------------------------------
# 5a. Kruskal-Wallis across the three clusters (BH over all features) + median per cluster
kw <- fl[, .(Kruskal_p = kruskal.test(Value ~ Cluster)$p.value), by = Feature]
kw[, Kruskal_FDR := p.adjust(Kruskal_p, method = "BH")]
med <- dcast(fl[, .(Median = median(Value)), by = .(Feature, Cluster)], Feature ~ Cluster, value.var = "Median")
med_cols <- paste0("Median_", sub(" ", "_", cluster_levels))
setnames(med, cluster_levels, med_cols)
kw <- merge(kw, med, by = "Feature")
kw[, Highest_cluster := cluster_levels[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = med_cols]
setorder(kw, Kruskal_FDR, Kruskal_p)
fwrite(kw, pre("clusters_KruskalWallis.tsv"), sep = "\t", quote = FALSE)
print(kw)

# 5b. Pairwise Wilcoxon (difference = first minus second; BH over all feature x pair tests)
pw <- rbindlist(lapply(combn(cluster_levels, 2, simplify = FALSE), function(pr) {
  fl[Cluster %in% pr, {
    a <- Value[Cluster == pr[1]]; b <- Value[Cluster == pr[2]]
    .(Comparison = paste0(sub(" ", "", pr[1]), "_vs_", sub(" ", "", pr[2])),
      Median_diff = median(a) - median(b),
      Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
  }, by = Feature]
}))
pw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
setorder(pw, Wilcoxon_FDR, Wilcoxon_p)
fwrite(pw, pre("clusters_pairwise_Wilcoxon.tsv"), sep = "\t", quote = FALSE)

# 5c. FOCUS cluster vs all other samples (Wilcoxon, all samples; BH over features)
fv <- fl[, {
  a <- Value[Focus]; b <- Value[!Focus]
  .(n_focus = length(a), n_rest = length(b), Median_focus = median(a), Median_rest = median(b),
    Median_diff = median(a) - median(b), Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
}, by = Feature]
fv[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
setorder(fv, Wilcoxon_FDR, Wilcoxon_p)
fwrite(fv, pre("focus_cluster_vs_rest_Wilcoxon.tsv"), sep = "\t", quote = FALSE)
print(fv)

# 5d. TILs-adjusted models: is the FOCUS-cluster difference more than a TILs effect?
res <- NULL
if (HAS_TILS) {
  have <- sample_tab[is.finite(TILs)]
  message(sprintf("TILs models: %d samples with a TILs value (%d in %s, %d in the other clusters)",
                  nrow(have), sum(have$Focus), FOCUS, sum(!have$Focus)))
  if (nrow(sample_tab) > nrow(have)) {
    message("Samples WITHOUT a TILs value:"); print(sample_tab[!is.finite(TILs), .(Sample, Cluster)])
  }
  # adjustment is only reliable if the clusters overlap in TILs
  til_by_cluster <- have[, .(n = .N, TILs_median = median(TILs), TILs_min = min(TILs), TILs_max = max(TILs)), by = Cluster][order(Cluster)]
  fwrite(til_by_cluster, pre("TILs_by_cluster.tsv"), sep = "\t", quote = FALSE)
  print(til_by_cluster)

  fit_feature <- function(dd, f, tils_term) {
    m <- data.frame(y = dd[[f]], Focus = dd$Focus, TILs = dd$TILs)
    m <- m[complete.cases(m), ]
    if (nrow(m) < 8 || length(unique(m$Focus)) < 2 || sum(m$Focus) < 3) return(NULL)
    cu <- summary(lm(y ~ Focus, m))$coefficients
    ca <- summary(lm(reformulate(c("Focus", tils_term), response = "y"), m))$coefficients
    data.table(Feature = f, n = nrow(m), n_focus = sum(m$Focus), n_rest = sum(!m$Focus),
               Effect_unadjusted = cu["FocusTRUE", "Estimate"], Unadjusted_p = cu["FocusTRUE", "Pr(>|t|)"],
               Effect_adjusted = ca["FocusTRUE", "Estimate"],  Adjusted_p = ca["FocusTRUE", "Pr(>|t|)"],
               TILs_coefficient = ca[3, "Estimate"], TILs_p = ca[3, "Pr(>|t|)"])
  }
  model_specs <- list(linear             = list(term = "TILs",           drop_top = FALSE),
                      log2               = list(term = "log2(TILs + 1)", drop_top = FALSE),
                      linear_no_top_TILs = list(term = "TILs",           drop_top = TRUE))
  res <- rbindlist(lapply(names(model_specs), function(nm) {
    sp <- model_specs[[nm]]
    dd <- if (sp$drop_top) have[TILs < max(TILs)] else have
    r  <- rbindlist(lapply(feature_names, function(f) fit_feature(dd, f, sp$term)))
    r[, Model := nm]
    r[, Unadjusted_FDR := p.adjust(Unadjusted_p, method = "BH")]
    r[, Adjusted_FDR   := p.adjust(Adjusted_p,   method = "BH")]
    r[, Change_in_effect_pct := 100 * (Effect_adjusted - Effect_unadjusted) / abs(Effect_unadjusted)]
    r[, Holds_after_adjustment := Adjusted_FDR < FDR_CUT & sign(Effect_adjusted) == sign(Effect_unadjusted)]
    r
  }))
  setcolorder(res, c("Model", "Feature"))
  res[, Feature := factor(Feature, levels = feature_names)]
  setorder(res, Model, Feature)
  res[, Feature := as.character(Feature)]
  fwrite(res, pre("TILs_adjusted_models.tsv"), sep = "\t", quote = FALSE)
  message(FOCUS, " effect before and after adjusting for TILs (linear model):")
  print(res[Model == "linear", .(Feature, n, Effect_unadjusted = round(Effect_unadjusted, 2),
                                 Effect_adjusted = round(Effect_adjusted, 2),
                                 Adjusted_FDR = signif(Adjusted_FDR, 2), Holds_after_adjustment)])
} else message("TILs missing or too few values: TILs-adjusted models skipped.")

# 5e. Grade 2 vs Grade 3 (Wilcoxon; difference = first minus second)
sl_g <- fl[Grade %in% GRADE_PAIR]
if (all(table(factor(sl_g[Feature == feature_names[1], Grade], levels = GRADE_PAIR)) >= 3)) {
  gw <- sl_g[, {
    a <- Value[Grade == GRADE_PAIR[1]]; b <- Value[Grade == GRADE_PAIR[2]]
    .(Median_diff = median(a) - median(b), Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
  }, by = Feature]
  gw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
  setorder(gw, Wilcoxon_FDR, Wilcoxon_p)
  fwrite(gw, pre("grade_Wilcoxon.tsv"), sep = "\t", quote = FALSE)
} else message("Fewer than 3 samples in a grade group: grade test skipped.")

# 5f. Spearman correlation with the clinical variables (BH over all tests)
cor_tab <- rbindlist(lapply(feature_names, function(s) rbindlist(lapply(clin_vars, function(v) {
  x <- feat[s, rownames(metadata)]; y <- as.numeric(metadata[[v]])
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.table(Feature = s, Variable = v, n = sum(ok), Spearman_rho = unname(ct$estimate), p = ct$p.value)
}))))
if (nrow(cor_tab) > 0) {
  cor_tab[, FDR := p.adjust(p, method = "BH")]
  setorder(cor_tab, FDR, p)
  fwrite(cor_tab, pre("clinical_correlation_Spearman.tsv"), sep = "\t", quote = FALSE)
}

# ------------------------------- 6. FIGURES ----------------------------------
# 6a. Heatmaps (row z-score clipped at +/-2.5): (i) rows and samples clustered, (ii) samples grouped by cluster
ann <- data.frame(Cluster = metadata$Cluster, Grade = factor(metadata$Grade, levels = names(grade_colors)),   # blank/other grades -> NA
                  metadata[, intersect("TILs", clin_vars), drop = FALSE], row.names = rownames(metadata))
ann <- ann[, colSums(!is.na(ann)) > 0, drop = FALSE]
ann_cols <- list(Cluster = cluster_colors, Grade = grade_colors)
ann_cols <- ann_cols[names(ann_cols) %in% names(ann)]
z_all <- t(scale(t(feat)))
z_all[z_all >  2.5] <-  2.5
z_all[z_all < -2.5] <- -2.5
z_all <- z_all[complete.cases(z_all), , drop = FALSE]

plot_hm <- function(z, main, stem, by_cluster = FALSE) {
  gaps <- NULL
  if (by_cluster) {
    grp <- as.character(ann[colnames(z), "Cluster"])
    ord <- unlist(lapply(cluster_levels, function(g) {
      cols <- colnames(z)[grp == g]
      if (length(cols) > 2) cols[hclust(dist(t(z[, cols, drop = FALSE])), "ward.D2")$order] else cols
    }))
    z    <- z[, ord, drop = FALSE]
    gaps <- as.integer(head(cumsum(table(factor(as.character(ann[ord, "Cluster"]), levels = cluster_levels))), -1))
  }
  draw <- function() {
    pheatmap::pheatmap(z, silent = TRUE, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
                       annotation_col = ann, annotation_colors = ann_cols,
                       cluster_rows = TRUE, cluster_cols = !by_cluster, gaps_col = gaps,
                       clustering_method = "ward.D2", fontsize_row = 9, fontsize_col = 7,
                       border_color = NA, main = main)
  }
  save_pheatmap(draw, stem, 10, 5.5)
}
plot_hm(z_all, "cGAS-STING panel and scores (row z-score; rows and samples clustered)", pre("heatmap_clustered"))
plot_hm(z_all, "cGAS-STING panel and scores (row z-score; samples grouped by NMF cluster)",
        pre("heatmap_ordered_by_cluster"), by_cluster = TRUE)

# 6b. Boxplots by cluster: the scores and the individual genes (Kruskal-Wallis FDR in each title)
lab <- kw[, .(Feature, Label = paste0(Feature, "\nKruskal-Wallis FDR = ", fmt_p(Kruskal_FDR)))]
sl  <- merge(fl, lab, by = "Feature")
sl[, Label := factor(Label, levels = lab[match(feature_names, Feature), Label])]
box_layer <- list(geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6),
                  geom_jitter(width = 0.12, size = 1.8, shape = 21, fill = "white"),
                  scale_fill_manual(values = cluster_colors),
                  theme_bw(base_size = 11),
                  theme(legend.position = "none", panel.grid.minor = element_blank()))
p_facet <- ggplot(sl, aes(Cluster, Value, fill = Cluster)) + box_layer +
  facet_wrap(~ Label, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Score (z-score units) or VST expression") +
  theme(strip.text = element_text(size = 8), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_facet, pre("boxplots_by_cluster_facet"), 11, 3.2 * ceiling(length(feature_names) / 3) + 0.6)

p_main <- ggplot(fl[Feature == main_score], aes(Cluster, Value, fill = Cluster)) + box_layer +
  labs(title = "cGAS-STING score",
       subtitle = paste0("Kruskal-Wallis FDR = ", fmt_p(kw[Feature == main_score, Kruskal_FDR]),
                         "\n", FOCUS, " vs rest Wilcoxon FDR = ", fmt_p(fv[Feature == main_score, Wilcoxon_FDR])),
       x = NULL, y = "Mean z-score of the panel genes") +
  theme(plot.title = element_text(face = "bold"))
save_plot(p_main, pre("score_by_cluster"), 5, 4.5)

# 6c. Scores vs TILs with one fitted line for the FOCUS cluster and one for the other samples
if (HAS_TILS) {
  long <- melt(have, id.vars = c("Sample", "Cluster", "TILs", "Focus"), measure.vars = main_scores,
               variable.name = "Score", value.name = "Value")
  long[, Score := factor(Score, levels = main_scores)]
  ct_main <- suppressWarnings(cor.test(have[[main_score]], have$TILs, method = "spearman", exact = FALSE))
  p_sc <- ggplot(long, aes(TILs, Value)) +
    geom_smooth(aes(group = Focus, linetype = Focus), method = "lm", se = FALSE, colour = "grey30", linewidth = 0.6) +
    geom_point(aes(colour = Cluster), size = 2.8) +
    scale_colour_manual(values = cluster_colors) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                          labels = c("TRUE" = FOCUS, "FALSE" = "Other clusters"), name = "Fitted line") +
    facet_wrap(~ Score, scales = "free_y", ncol = 3) +
    labs(x = "TILs (%)", y = "Score (mean z-score)",
         title = paste0("cGAS-STING scores vs TILs (", nrow(have), " samples with TILs)"),
         subtitle = paste0("cGAS-STING score vs TILs: Spearman rho = ", round(unname(ct_main$estimate), 2),
                           ", p = ", fmt_p(ct_main$p.value),
                           ". Lines one above the other = the ", FOCUS, " difference persists at the same TILs")) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  save_plot(p_sc, pre("score_vs_TILs_by_group"), 11, 4.3)

  # 6d. Cluster effect before vs after adjustment (linear model); * = adjusted FDR < FDR_CUT
  pl <- melt(res[Model == "linear"], id.vars = c("Feature", "Adjusted_FDR"),
             measure.vars = c("Effect_unadjusted", "Effect_adjusted"), variable.name = "Type", value.name = "Effect")
  pl[, Type := factor(fifelse(Type == "Effect_unadjusted", "Unadjusted", "Adjusted for TILs"), levels = c("Unadjusted", "Adjusted for TILs"))]
  pl[, Feature := factor(Feature, levels = rev(feature_names))]
  p_eff <- ggplot(pl, aes(Effect, Feature, fill = Type)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(data = unique(pl[Type == "Adjusted for TILs", .(Feature, Effect, Adjusted_FDR)]),
              aes(label = ifelse(Adjusted_FDR < FDR_CUT, "*", ""), fill = NULL), hjust = -0.3, size = 6, show.legend = FALSE) +
    scale_fill_manual(values = c("Unadjusted" = "grey65", "Adjusted for TILs" = "#D95F02")) +
    labs(x = paste0(FOCUS, " minus other clusters (z-score units for scores, VST units for genes)"), y = NULL, fill = NULL,
         title = paste0(FOCUS, " effect before and after adjusting for TILs"),
         subtitle = paste0("* adjusted FDR < ", FDR_CUT, " (BH over features); linear model, n = ", nrow(have))) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  save_plot(p_eff, pre("effect_before_after_adjustment"), 8, 5)
}

message("Done. Outputs in: ", OUT_ROOT)