###############################################################################
# pILC RNA-seq | script: GSVA (Hallmark) across the NMF clusters and Grade 2 vs 3
# Reads the outputs of pILC_01a (Sample_metadata.tsv, VST matrix). Run it after pILC_02 if you want the
# comparison with ssGSEA (that part is skipped when the ssGSEA files are missing).
#   Scores : GSVA (gsvaParam, Gaussian kernel for the continuous VST values), Hallmark sets
#   Tests  : Kruskal-Wallis + pairwise Wilcoxon (BH-FDR) as in the ssGSEA analysis, so the two methods
#            can be compared directly; limma on the scores (moderated F and pairwise contrasts) as the
#            second, model-based test; Grade 2 vs 3 (Wilcoxon and limma)
#   Compare: per-pathway correlation of GSVA and ssGSEA scores; which pathways are significant in both
# Note   : clusters were derived from the same expression data, so p-values describe how the groups
#          differ and are not independent validation.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)         # needs >= 1.50 (gsvaParam)
  library(ggplot2)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
ROOT    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
DE_ROOT <- file.path(ROOT, "DESeq2_stranded")
SHARED  <- file.path(DE_ROOT, "Shared_inputs_for_pathway_analysis")
META_FILE    <- file.path(SHARED, "Sample_metadata.tsv")
VST_FILE     <- file.path(SHARED, "VST_blind_ssGSEA_input.tsv")
HALLMARK_GMT <- file.path(ROOT, "Pathway_analysis", "Dataset", "h.all.v2025.1.Hs.symbols.gmt")

# ssGSEA results from pILC_02 (only used for the comparison)
SS_SCORES <- file.path(ROOT, "Pathway_analysis", "ssGSEA", "ssGSEA_Hallmark_scores.csv")
SS_KW     <- file.path(ROOT, "Pathway_analysis", "ssGSEA", "ssGSEA_Hallmark_clusters_KruskalWallis.csv")

OUT_ROOT <- file.path(ROOT, "Pathway_analysis", "GSVA")              # no dates in names
OUT_IND  <- file.path(OUT_ROOT, "GSVA_boxplots_clusters_individual")
dir.create(OUT_IND, recursive = TRUE, showWarnings = FALSE)

MIN_SIZE <- 15          # same gene-set size limits as the fgsea analysis in pILC_02
MAX_SIZE <- 500
KCDF     <- "Gaussian"  # VST values are continuous (use "Poisson" only for integer counts)
FDR_CUT  <- 0.05

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
grade_colors   <- c("Grade 2" = "#8BC34A", "Grade 3" = "#6A1B9A")
GRADE_PAIR     <- names(grade_colors)

for (f in c(META_FILE, VST_FILE, HALLMARK_GMT))
  if (!file.exists(f)) stop("Missing input: ", f, if (f != HALLMARK_GMT) "  (run pILC_01a first)" else "")
if (packageVersion("GSVA") < "1.50.0") stop("GSVA >= 1.50 is needed (gsvaParam).")

# ------------------------------- 2. HELPERS ----------------------------------
save_plot <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
}
save_pheatmap <- function(draw_fun, stem, w, h) {
  png(paste0(stem, ".png"), width = w * 300, height = h * 300, res = 300)
  draw_fun(); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = w, height = h)
  draw_fun(); invisible(dev.off())
}
fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), formatC(x, format = "f", digits = 3))

# limma on the pathway scores: pairwise contrasts (BH over all pathway x pair tests) and, for 3 groups, the moderated F test
limma_pairs <- function(scores, group, levs) {
  g  <- factor(group, levels = levs)
  lv <- paste0("G", seq_along(levs))
  design <- model.matrix(~ 0 + g)
  colnames(design) <- lv
  prs <- combn(seq_along(levs), 2, simplify = FALSE)
  cm  <- limma::makeContrasts(contrasts = vapply(prs, function(p) paste0(lv[p[1]], "-", lv[p[2]]), character(1)), levels = design)
  nm  <- vapply(prs, function(p) paste0(gsub(" ", "", levs[p[1]]), "_vs_", gsub(" ", "", levs[p[2]])), character(1))
  colnames(cm) <- nm
  fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(scores, design), cm))
  overall <- NULL
  if (length(levs) > 2) {
    tt <- limma::topTable(fit, number = Inf, sort.by = "none")
    overall <- data.table(Pathway = rownames(tt), limma_F = tt$F, limma_F_p = tt$P.Value)
    overall[, limma_F_FDR := p.adjust(limma_F_p, method = "BH")]
    setorder(overall, limma_F_FDR, limma_F_p)
  }
  pw <- rbindlist(lapply(nm, function(n) {
    tt <- limma::topTable(fit, coef = n, number = Inf, sort.by = "none")
    data.table(Pathway = rownames(tt), Comparison = n, Mean_score_diff = tt$logFC, t = tt$t, limma_p = tt$P.Value)
  }))
  pw[, limma_FDR := p.adjust(limma_p, method = "BH")]
  setorder(pw, limma_FDR, limma_p)
  list(pairs = pw, overall = overall)
}

# ------------------------------- 3. LOAD DATA --------------------------------
meta_dt <- fread(META_FILE)
metadata <- as.data.frame(meta_dt)
rownames(metadata) <- metadata$Sample
metadata$Cluster <- factor(metadata$Cluster, levels = cluster_levels)
metadata$Grade   <- factor(metadata$Grade, levels = GRADE_PAIR)
message("Cluster sizes:"); print(table(metadata$Cluster))

vst_dt <- fread(VST_FILE)
expr <- as.matrix(vst_dt[, names(vst_dt)[-1], with = FALSE])
rownames(expr) <- vst_dt[[1]]
if (!all(rownames(metadata) %in% colnames(expr))) stop("Samples of Sample_metadata.tsv are missing in the VST matrix.")
expr <- expr[, rownames(metadata)]
expr <- expr[apply(expr, 1, sd) > 0, , drop = FALSE]                     # GSVA cannot use constant genes

gmt <- strsplit(readLines(HALLMARK_GMT), "\t")
hallmark <- setNames(lapply(gmt, function(x) unique(x[-(1:2)])), sub("^HALLMARK_", "", vapply(gmt, `[`, character(1), 1)))
coverage <- mean(unique(unlist(hallmark)) %in% rownames(expr))
message(sprintf("GSVA input: %d genes x %d samples; %.1f%% of Hallmark genes present", nrow(expr), ncol(expr), 100 * coverage))
if (coverage < 0.8) warning("Less than 80% of Hallmark genes are in the input: check the gene filter in pILC_01a.")

# ------------------------------- 4. GSVA -------------------------------------
gv <- gsva(gsvaParam(exprData = expr, geneSets = hallmark, minSize = MIN_SIZE, maxSize = MAX_SIZE, kcdf = KCDF))
message("GSVA scores: ", nrow(gv), " pathways x ", ncol(gv), " samples (scores range ", round(min(gv), 2), " to ", round(max(gv), 2), ")")
fwrite(as.data.table(gv, keep.rownames = "Pathway"), file.path(OUT_ROOT, "GSVA_Hallmark_scores.csv"))

score_long <- melt(as.data.table(t(gv), keep.rownames = "Sample"), id.vars = "Sample", variable.name = "Pathway", value.name = "Score")
score_long[, Pathway := as.character(Pathway)]
score_long[, Cluster := metadata[Sample, "Cluster"]]
score_long[, Grade   := metadata[Sample, "Grade"]]

# ------------------------------- 5. STATISTICS -------------------------------
# 5a. Kruskal-Wallis across the clusters (BH over all pathways) + median per cluster
kw <- score_long[, .(Kruskal_p = kruskal.test(Score ~ Cluster)$p.value), by = Pathway]
kw[, Kruskal_FDR := p.adjust(Kruskal_p, method = "BH")]
med <- dcast(score_long[, .(Median = median(Score)), by = .(Pathway, Cluster)], Pathway ~ Cluster, value.var = "Median")
med_cols <- paste0("Median_", sub(" ", "_", cluster_levels))
setnames(med, cluster_levels, med_cols)
kw <- merge(kw, med, by = "Pathway")
kw[, Highest_cluster := cluster_levels[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = med_cols]
setorder(kw, Kruskal_FDR, Kruskal_p)
fwrite(kw, file.path(OUT_ROOT, "GSVA_Hallmark_clusters_KruskalWallis.csv"))

# 5b. Pairwise Wilcoxon (difference = first minus second; BH over all pathway x pair tests)
pw <- rbindlist(lapply(combn(cluster_levels, 2, simplify = FALSE), function(pr) {
  score_long[Cluster %in% pr, {
    a <- Score[Cluster == pr[1]]; b <- Score[Cluster == pr[2]]
    .(Comparison = paste0(sub(" ", "", pr[1]), "_vs_", sub(" ", "", pr[2])),
      Median_diff = median(a) - median(b), Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
  }, by = Pathway]
}))
pw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
setorder(pw, Wilcoxon_FDR, Wilcoxon_p)
fwrite(pw, file.path(OUT_ROOT, "GSVA_Hallmark_clusters_pairwise_Wilcoxon.csv"))

# 5c. Grade 2 vs Grade 3 (Wilcoxon; difference = first minus second)
sl_g <- score_long[Grade %in% GRADE_PAIR]
gw <- sl_g[, {
  a <- Score[Grade == GRADE_PAIR[1]]; b <- Score[Grade == GRADE_PAIR[2]]
  .(Median_diff = median(a) - median(b), Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
}, by = Pathway]
gw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
gw[, Higher_in := fifelse(Median_diff > 0, GRADE_PAIR[1], fifelse(Median_diff < 0, GRADE_PAIR[2], "tie"))]
setorder(gw, Wilcoxon_FDR, Wilcoxon_p)
fwrite(gw, file.path(OUT_ROOT, "GSVA_Hallmark_grade_Wilcoxon.csv"))

# 5d. limma on the GSVA scores (second, model-based test)
if (requireNamespace("limma", quietly = TRUE)) {
  lc <- limma_pairs(gv, metadata$Cluster, cluster_levels)
  fwrite(lc$overall, file.path(OUT_ROOT, "GSVA_Hallmark_clusters_limma_overall.csv"))
  fwrite(lc$pairs,   file.path(OUT_ROOT, "GSVA_Hallmark_clusters_limma_pairwise.csv"))
  gsel <- metadata$Grade %in% GRADE_PAIR
  lg <- limma_pairs(gv[, gsel, drop = FALSE], metadata$Grade[gsel], GRADE_PAIR)
  fwrite(lg$pairs, file.path(OUT_ROOT, "GSVA_Hallmark_grade_limma.csv"))
  message("limma significant pathways (FDR < ", FDR_CUT, "): clusters (moderated F) = ", sum(lc$overall$limma_F_FDR < FDR_CUT),
          "; Grade = ", sum(lg$pairs$limma_FDR < FDR_CUT))
} else message("Package limma not installed: the limma tests are skipped.")

message("GSVA significant pathways (FDR < ", FDR_CUT, "): clusters (Kruskal-Wallis) = ", sum(kw$Kruskal_FDR < FDR_CUT),
        "; ", GRADE_PAIR[1], " vs ", GRADE_PAIR[2], " (Wilcoxon) = ", sum(gw$Wilcoxon_FDR < FDR_CUT))

# ------------------------------- 6. GSVA vs ssGSEA ---------------------------
if (file.exists(SS_SCORES) && file.exists(SS_KW)) {
  ss <- fread(SS_SCORES)
  ss_m <- as.matrix(ss[, names(ss)[-1], with = FALSE]); rownames(ss_m) <- ss[[1]]
  common <- intersect(rownames(gv), rownames(ss_m))
  if (all(colnames(gv) %in% colnames(ss_m)) && length(common) > 0) {
    ss_kw <- fread(SS_KW)[, .(Pathway, ssGSEA_Kruskal_FDR = Kruskal_FDR, ssGSEA_Highest_cluster = Highest_cluster)]
    cmp <- data.table(Pathway = common,
                      Spearman_rho_scores = vapply(common, function(p) cor(gv[p, ], ss_m[p, colnames(gv)], method = "spearman"), numeric(1)))
    cmp <- merge(cmp, kw[, .(Pathway, GSVA_Kruskal_FDR = Kruskal_FDR, GSVA_Highest_cluster = Highest_cluster)], by = "Pathway")
    cmp <- merge(cmp, ss_kw, by = "Pathway")
    cmp[, Agreement := fifelse(GSVA_Kruskal_FDR < FDR_CUT & ssGSEA_Kruskal_FDR < FDR_CUT, "significant in both",
                       fifelse(GSVA_Kruskal_FDR < FDR_CUT, "GSVA only",
                       fifelse(ssGSEA_Kruskal_FDR < FDR_CUT, "ssGSEA only", "neither")))]
    cmp[, Same_highest_cluster := GSVA_Highest_cluster == ssGSEA_Highest_cluster]
    setorder(cmp, -Spearman_rho_scores)
    fwrite(cmp, file.path(OUT_ROOT, "GSVA_vs_ssGSEA_Hallmark_comparison.csv"))
    print(table(cmp$Agreement))
    cmp[, Pathway := factor(Pathway, levels = rev(cmp$Pathway))]
    p_cmp <- ggplot(cmp, aes(Spearman_rho_scores, Pathway, fill = Agreement)) +
      geom_col() +
      scale_fill_manual(values = c("significant in both" = "#B2182B", "GSVA only" = "#2166AC",
                                   "ssGSEA only" = "#E08214", "neither" = "grey75")) +
      coord_cartesian(xlim = c(min(0, min(cmp$Spearman_rho_scores)), 1)) +
      labs(x = "Spearman correlation of GSVA and ssGSEA scores across samples", y = NULL,
           title = "GSVA vs ssGSEA (Hallmark)", subtitle = paste0("Colour: Kruskal-Wallis FDR < ", FDR_CUT, " across clusters")) +
      theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
    save_plot(p_cmp, file.path(OUT_ROOT, "GSVA_vs_ssGSEA_Hallmark_comparison"), 8, 9)
  } else message("ssGSEA scores do not match the GSVA samples or pathways: comparison skipped.")
} else message("ssGSEA files from pILC_02 not found: comparison skipped.")

# ------------------------------- 7. FIGURES ----------------------------------
# 7a. Heatmap of all Hallmark pathways (row z-score, samples ordered by cluster then grade)
ord <- order(metadata$Cluster, metadata$Grade)
z <- t(scale(t(gv[, ord, drop = FALSE])))
z[z >  2.5] <-  2.5
z[z < -2.5] <- -2.5
ann <- data.frame(Cluster = metadata$Cluster, Grade = metadata$Grade, row.names = rownames(metadata))
draw_heatmap <- function() {
  pheatmap(z, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           annotation_col = ann, annotation_colors = list(Cluster = cluster_colors, Grade = grade_colors),
           cluster_rows = TRUE, cluster_cols = FALSE, gaps_col = head(cumsum(table(metadata$Cluster)), -1),
           fontsize_row = 8, fontsize_col = 7, border_color = NA,
           main = "GSVA Hallmark scores (row z-score; samples ordered by cluster, then grade)")
}
save_pheatmap(draw_heatmap, file.path(OUT_ROOT, "GSVA_Hallmark_heatmap"), 10, 9)

# 7b. Boxplots of pathways with Kruskal-Wallis FDR < FDR_CUT (cluster) or Wilcoxon FDR < FDR_CUT (grade)
box_layers <- function(cols) list(
  geom_boxplot(outlier.shape = NA, alpha = 0.4, colour = "black", linewidth = 0.8),
  geom_jitter(width = 0.15, size = 2, colour = "black", alpha = 0.8),
  scale_fill_manual(values = cols))
boxplots_for <- function(sl, stats_dt, fdr_col, group_var, cols, levs, stem, ind_dir = NULL) {
  sig <- stats_dt[get(fdr_col) < FDR_CUT][order(get(fdr_col))]
  if (nrow(sig) == 0) { message("No significant pathway for ", stem, " (", fdr_col, " < ", FDR_CUT, ")."); return(invisible(NULL)) }
  d <- sl[Pathway %in% sig$Pathway]
  d[, (group_var) := factor(get(group_var), levels = levs)]
  if (!is.null(ind_dir)) for (i in seq_len(nrow(sig))) {
    p <- ggplot(d[Pathway == sig$Pathway[i]], aes(.data[[group_var]], Score, fill = .data[[group_var]])) + box_layers(cols) +
      labs(title = paste0("HALLMARK_", sig$Pathway[i], "\nFDR-adjusted p = ", fmt_p(sig[[fdr_col]][i])), x = NULL, y = "GSVA Enrichment Score") +
      theme_classic(base_size = 16) +
      theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5, size = 16), axis.title = element_text(face = "bold"))
    save_plot(p, file.path(ind_dir, paste0("GSVA_", sig$Pathway[i])), 6.5, 5.5)
  }
  lab <- setNames(paste0(sig$Pathway, "\nFDR = ", vapply(sig[[fdr_col]], fmt_p, character(1))), sig$Pathway)
  d[, Panel := factor(lab[Pathway], levels = unname(lab))]
  pf <- ggplot(d, aes(.data[[group_var]], Score, fill = .data[[group_var]])) + box_layers(cols) +
    facet_wrap(~ Panel, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = "GSVA Enrichment Score", title = paste0("Hallmark pathways with FDR < ", FDR_CUT)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9), panel.grid.minor = element_blank())
  save_plot(pf, file.path(OUT_ROOT, stem), 12, 3.3 * ceiling(nrow(sig) / 3) + 0.8)
}
boxplots_for(score_long, kw, "Kruskal_FDR", "Cluster", cluster_colors, cluster_levels, "GSVA_boxplots_clusters_facet", OUT_IND)
boxplots_for(sl_g, gw, "Wilcoxon_FDR", "Grade", grade_colors, GRADE_PAIR, "GSVA_boxplots_grade_facet")

message("Done. Outputs in: ", OUT_ROOT)