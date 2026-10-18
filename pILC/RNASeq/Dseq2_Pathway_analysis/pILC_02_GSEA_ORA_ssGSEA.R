###############################################################################
# pILC RNA-seq | script 2: GSEA, ORA and ssGSEA (MSigDB Hallmark)
# Reads the outputs of pILC_01a_DESeq2_clusters.R and pILC_01b_DESeq2_grade.R
#   GSEA  : all tested genes ranked by the unshrunken Wald statistic, every contrast
#   ORA   : significant up / down genes (FDR < 0.05 and |shrunken log2FC| >= 1)
#   ssGSEA: per-sample Hallmark scores; clusters (Kruskal-Wallis) and grades (Wilcoxon)
# Positive NES / log2FC / difference = higher in the FIRST-named group of a contrast.
# Note  : clusters were derived from the same expression data, so p-values describe
#         how the groups differ and are not independent validation.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(clusterProfiler)
  library(GSVA)         # needs >= 1.50 (ssgseaParam)
  library(ggplot2)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
ROOT    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
DE_ROOT <- file.path(ROOT, "DESeq2_stranded")
SHARED  <- file.path(DE_ROOT, "Shared_inputs_for_pathway_analysis")

CLUSTER_RES  <- file.path(DE_ROOT, "DESeq2_Clusters_R_NMF_k3", "DESeq2_all_pairwise_results.tsv")
GRADE_RES    <- file.path(DE_ROOT, "DESeq2_Grade2_vs_Grade3", "DESeq2_Grade2_vs_Grade3_results.tsv")
META_FILE    <- file.path(SHARED, "Sample_metadata.tsv")
VST_FILE     <- file.path(SHARED, "VST_blind_ssGSEA_input.tsv")
HALLMARK_GMT <- file.path(ROOT, "Pathway_analysis", "Dataset", "h.all.v2025.1.Hs.symbols.gmt")
for (f in c(CLUSTER_RES, META_FILE, VST_FILE))
  if (!file.exists(f)) stop("Missing input (run pILC_01a first): ", f)
if (!file.exists(GRADE_RES)) stop("Missing input (run pILC_01b first): ", GRADE_RES)
if (!file.exists(HALLMARK_GMT)) stop("Hallmark gene-set file not found: ", HALLMARK_GMT)

OUT_ROOT <- file.path(ROOT, "Pathway_analysis")                        # no dates in names
OUT_GS   <- file.path(OUT_ROOT, "GSEA")
OUT_OR   <- file.path(OUT_ROOT, "ORA")
OUT_SS   <- file.path(OUT_ROOT, "ssGSEA")
dir.create(file.path(OUT_GS, "enrichment_individual"), recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_OR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SS, recursive = TRUE, showWarnings = FALSE)

FDR_CUT          <- 0.05   # pathway-level FDR (GSEA, ORA, ssGSEA)
MIN_GENES_ORA    <- 10     # ORA: skip gene lists shorter than this
MAX_FACET_PANELS <- 12     # GSEA enrichment facet: panels per contrast (most significant first)

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
grade_colors   <- c("Grade 2" = "#8BC34A", "Grade 3" = "#6A1B9A")
cluster_levels <- names(cluster_colors)

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
fmt_p <- function(x) {
  if (is.na(x)) "NA" else if (x < 1e-4) format(x, digits = 3, scientific = TRUE)
  else format(signif(x, 3), scientific = FALSE)
}
# "Cluster1_vs_Cluster2" -> c("Cluster 1", "Cluster 2");  "Grade2_vs_Grade3" -> c("Grade 2", "Grade 3")
cmp_groups <- function(cmp) gsub("(Cluster|Grade)([0-9]+)", "\\1 \\2", strsplit(cmp, "_vs_")[[1]])

# ------------------------------- 3. LOAD DATA --------------------------------
cluster_dt <- fread(CLUSTER_RES)
grade_dt   <- fread(GRADE_RES)
de <- rbindlist(list(cluster_dt, grade_dt), use.names = TRUE, fill = TRUE)
need <- c("Comparison", "gene_symbol", "Wald_stat", "padj", "log2FoldChange", "Significant_FDR05_LFC1")
if (!all(need %in% names(de))) stop("Result tables lack columns: ", paste(setdiff(need, names(de)), collapse = ", "))
contrast_names  <- unique(de$Comparison)
GRADE_PAIR <- cmp_groups(unique(grade_dt$Comparison))
message("Contrasts: ", paste(contrast_names, collapse = ", "))

gmt      <- strsplit(readLines(HALLMARK_GMT), "\t")
hallmark <- setNames(lapply(gmt, function(x) unique(x[-(1:2)])),
                     sub("^HALLMARK_", "", vapply(gmt, `[`, character(1), 1)))
term2gene <- as.data.frame(data.table(term = rep(names(hallmark), lengths(hallmark)),
                                      gene = unlist(hallmark, use.names = FALSE)))
message(sprintf("%.1f%% of Hallmark genes are among the tested genes",
                100 * mean(unique(unlist(hallmark)) %in% unique(de$gene_symbol))))

# =============================================================================
# 4. GSEA (ranked by the unshrunken Wald statistic)
# =============================================================================
ranks_list <- list()
gsea_list  <- list()
for (cmp in contrast_names) {
  d <- de[Comparison == cmp & !is.na(Wald_stat) & !is.na(gene_symbol) & gene_symbol != ""]
  d <- d[order(-abs(Wald_stat))][!duplicated(gene_symbol)]      # duplicated symbols: keep the larger |stat|
  ranks <- sort(setNames(d$Wald_stat, d$gene_symbol), decreasing = TRUE)
  ranks_list[[cmp]] <- ranks
  set.seed(SEED)
  res <- fgseaMultilevel(pathways = hallmark, stats = ranks, minSize = 15, maxSize = 500, eps = 0)
  res[, Comparison := cmp]
  gsea_list[[cmp]] <- res
}
gsea_all <- rbindlist(gsea_list)
gsea_out <- copy(gsea_all)
gsea_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
setorder(gsea_out, Comparison, padj)
fwrite(gsea_out, file.path(OUT_GS, "GSEA_Hallmark_results.tsv"), sep = "\t", quote = FALSE)

gsea_summary <- gsea_all[, .(Pathways_tested = .N,
                             FDR05 = sum(padj < FDR_CUT, na.rm = TRUE),
                             Positive_NES = sum(padj < FDR_CUT & NES > 0, na.rm = TRUE),
                             Negative_NES = sum(padj < FDR_CUT & NES < 0, na.rm = TRUE)),
                         by = Comparison]
fwrite(gsea_summary, file.path(OUT_GS, "GSEA_Hallmark_summary.tsv"), sep = "\t", quote = FALSE)
message("GSEA significant pathways (FDR < ", FDR_CUT, ") per contrast:"); print(gsea_summary)

# 4a. NES dot plot across contrast_names (pathways significant in at least one contrast)
sig_any <- unique(gsea_all[padj < FDR_CUT, pathway])
if (length(sig_any) > 0) {
  pd <- gsea_all[pathway %in% sig_any]
  pd[, sig := !is.na(padj) & padj < FDR_CUT]
  pd[, FDR_plot := fifelse(is.na(padj), 1, pmax(padj, 1e-10))]
  ord <- pd[, .(m = mean(NES)), by = pathway][order(m), pathway]
  pd[, pathway := factor(pathway, levels = ord)]
  p_dot <- ggplot(pd, aes(Comparison, pathway)) +
    geom_point(aes(size = -log10(FDR_plot), fill = NES, stroke = ifelse(sig, 1.3, 0.3)),
               shape = 21, colour = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(x = NULL, y = NULL, size = "-log10(FDR)", title = "GSEA Hallmark (ranked by DESeq2 Wald statistic)",
         subtitle = "Thick outline: FDR < 0.05\nPositive NES = higher in the first-named group") +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_dot, file.path(OUT_GS, "GSEA_Hallmark_NES_dotplot"),
            max(8, 1.4 * length(contrast_names) + 4), max(4, 0.3 * length(sig_any) + 2))
} else message("No Hallmark pathway reached GSEA FDR < ", FDR_CUT, " in any contrast.")

# 4b. Enrichment plots for every significant pathway: individual files + one facet per contrast
facet_ok <- requireNamespace("patchwork", quietly = TRUE)
if (!facet_ok) message("Package 'patchwork' not installed: GSEA facet figures are skipped (individual plots are still made).")
for (cmp in contrast_names) {
  g   <- cmp_groups(cmp)
  sig <- gsea_all[Comparison == cmp & !is.na(padj) & padj < FDR_CUT][order(padj)]
  if (nrow(sig) == 0) next
  plots <- list()
  for (i in seq_len(nrow(sig))) {
    pth <- sig$pathway[i]
    p <- plotEnrichment(hallmark[[pth]], ranks_list[[cmp]]) +
      labs(title = paste0("HALLMARK_", pth),
           subtitle = sprintf("%s vs %s | NES = %.2f | FDR = %s", g[1], g[2], sig$NES[i], fmt_p(sig$padj[i])),
           x = "Gene rank (Wald statistic)", y = "Enrichment score") +
      theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(size = 9))
    save_plot(p, file.path(OUT_GS, "enrichment_individual", paste0(cmp, "_", pth)), 6, 4.5)
    plots[[pth]] <- p
  }
  if (facet_ok) {
    n_pan <- min(length(plots), MAX_FACET_PANELS)
    pf <- patchwork::wrap_plots(plots[seq_len(n_pan)], ncol = 3) +
      patchwork::plot_annotation(title = paste0("GSEA Hallmark: ", g[1], " vs ", g[2],
                                                " (FDR < ", FDR_CUT, "; top ", n_pan, " of ", length(plots), ")"))
    save_plot(pf, file.path(OUT_GS, paste0("GSEA_enrichment_facet_", cmp)), 15, 4.3 * ceiling(n_pan / 3) + 0.5)
  }
}

# =============================================================================
# 5. ORA (significant genes vs all tested genes)
# =============================================================================
ora_list <- list()
for (cmp in contrast_names) {
  g <- cmp_groups(cmp)
  d <- de[Comparison == cmp & !is.na(padj) & !is.na(gene_symbol) & gene_symbol != ""]
  universe <- unique(d$gene_symbol)                         # background = all genes tested in this contrast
  gene_sets <- list(
    up   = unique(d$gene_symbol[d$Significant_FDR05_LFC1 %in% TRUE & d$log2FoldChange > 0]),
    down = unique(d$gene_symbol[d$Significant_FDR05_LFC1 %in% TRUE & d$log2FoldChange < 0]))
  dir_label <- c(up = paste0("Higher_in_", g[1]), down = paste0("Higher_in_", g[2]))
  for (dn in names(gene_sets)) {
    gl <- gene_sets[[dn]]
    message(sprintf("ORA %s | %s: %d genes (background %d)", cmp, dir_label[[dn]], length(gl), length(universe)))
    ora_list[[length(ora_list) + 1]] <- data.table(Comparison = cmp, Direction = dir_label[[dn]],
                                                   n_input = length(gl), n_background = length(universe),
                                                   tested = length(gl) >= MIN_GENES_ORA)
    if (length(gl) < MIN_GENES_ORA) next
    e <- enricher(gene = gl, universe = universe, TERM2GENE = term2gene, pAdjustMethod = "BH",
                  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
    if (is.null(e) || nrow(e@result) == 0) next
    r <- as.data.table(e@result)
    r[, `:=`(Comparison = cmp, Direction = dir_label[[dn]], n_input = length(gl))]
    ora_list[[length(ora_list) + 1]] <- r
  }
}
ora_inputs  <- rbindlist(Filter(function(x) "tested" %in% names(x), ora_list))
ora_results <- rbindlist(Filter(function(x) "p.adjust" %in% names(x), ora_list), use.names = TRUE, fill = TRUE)
fwrite(ora_inputs, file.path(OUT_OR, "ORA_Hallmark_gene_list_sizes.tsv"), sep = "\t", quote = FALSE)

if (nrow(ora_results) > 0) {
  ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
  ora_results[, FoldEnrichment := ratio(GeneRatio) / ratio(BgRatio)]
  setorder(ora_results, Comparison, Direction, p.adjust)
  fwrite(ora_results, file.path(OUT_OR, "ORA_Hallmark_results.tsv"), sep = "\t", quote = FALSE)

  ora_sig <- ora_results[p.adjust < FDR_CUT]
  if (nrow(ora_sig) > 0) {
    # Short x-axis labels (Cluster -> C) so that neighbouring labels do not overlap
    ora_sig[, Group := paste0(gsub("Cluster", "C", Comparison), "\nhigher in ",
                              gsub("Cluster ", "C", sub("^Higher_in_", "", Direction)))]
    p_ora <- ggplot(ora_sig, aes(Group, reorder(ID, FoldEnrichment))) +
      geom_point(aes(size = Count, fill = -log10(p.adjust)), shape = 21, colour = "black") +
      scale_fill_gradient(low = "white", high = "#B2182B") +
      labs(x = NULL, y = NULL, fill = "-log10(FDR)", size = "Overlap (genes)", title = "ORA Hallmark",
           subtitle = "Significant genes vs tested background; shown: FDR < 0.05") +
      theme_bw(base_size = 12)
    save_plot(p_ora, file.path(OUT_OR, "ORA_Hallmark_dotplot"),
              max(8, 2.0 * length(unique(ora_sig$Group)) + 3),
              max(4, 0.3 * length(unique(ora_sig$ID)) + 2))
  } else message("No Hallmark pathway reached ORA FDR < ", FDR_CUT, ". See ORA_Hallmark_results.tsv.")
} else message("No gene list had >= ", MIN_GENES_ORA, " genes: ORA was skipped (see ORA_Hallmark_gene_list_sizes.tsv).")

# =============================================================================
# 6. ssGSEA (per-sample Hallmark scores)
# =============================================================================
if (packageVersion("GSVA") < "1.50.0") stop("GSVA >= 1.50 is needed (ssgseaParam).")

meta <- fread(META_FILE)
metadata <- as.data.frame(meta)
rownames(metadata) <- metadata$Sample
metadata$Cluster <- factor(metadata$Cluster, levels = cluster_levels)
metadata$Grade   <- factor(metadata$Grade, levels = GRADE_PAIR)

vst_dt <- fread(VST_FILE)
expr <- as.matrix(vst_dt[, names(vst_dt)[-1], with = FALSE])
rownames(expr) <- vst_dt[[1]]
if (!all(rownames(metadata) %in% colnames(expr))) stop("Samples of Sample_metadata.tsv are missing in the VST matrix.")
expr <- expr[, rownames(metadata)]
coverage <- mean(unique(unlist(hallmark)) %in% rownames(expr))
message(sprintf("ssGSEA input: %d genes x %d samples; %.1f%% of Hallmark genes present",
                nrow(expr), ncol(expr), 100 * coverage))
if (coverage < 0.8) warning("Less than 80% of Hallmark genes are in the ssGSEA input: check the gene filter in pILC_01a_DESeq2_clusters.R.")

ss_scores <- gsva(ssgseaParam(exprData = expr, geneSets = hallmark, normalize = TRUE))
fwrite(as.data.table(ss_scores, keep.rownames = "Pathway"), file.path(OUT_SS, "ssGSEA_Hallmark_scores.csv"))

score_long <- melt(as.data.table(t(ss_scores), keep.rownames = "Sample"),
                   id.vars = "Sample", variable.name = "Pathway", value.name = "Score")
score_long[, Pathway := as.character(Pathway)]
score_long[, Cluster := metadata[Sample, "Cluster"]]
score_long[, Grade   := metadata[Sample, "Grade"]]

# 6a. Across the three clusters: Kruskal-Wallis (primary, BH over all pathways) + pairwise Wilcoxon
kw <- score_long[, .(Kruskal_p = kruskal.test(Score ~ Cluster)$p.value), by = Pathway]
kw[, Kruskal_FDR := p.adjust(Kruskal_p, method = "BH")]
med <- dcast(score_long[, .(Median = median(Score)), by = .(Pathway, Cluster)],
             Pathway ~ Cluster, value.var = "Median")
med_cols <- paste0("Median_", sub(" ", "_", cluster_levels))
setnames(med, cluster_levels, med_cols)
kw <- merge(kw, med, by = "Pathway")
kw[, Highest_cluster := cluster_levels[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = med_cols]
setorder(kw, Kruskal_FDR, Kruskal_p)
fwrite(kw, file.path(OUT_SS, "ssGSEA_Hallmark_clusters_KruskalWallis.csv"))

pw <- rbindlist(lapply(combn(cluster_levels, 2, simplify = FALSE), function(pr) {
  score_long[Cluster %in% pr, {
    a <- Score[Cluster == pr[1]]; b <- Score[Cluster == pr[2]]
    .(Comparison = paste0(sub(" ", "", pr[1]), "_vs_", sub(" ", "", pr[2])),
      Median_diff = median(a) - median(b),
      Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
  }, by = Pathway]
}))
pw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]      # over all pathway x pair tests
setorder(pw, Wilcoxon_FDR, Wilcoxon_p)
fwrite(pw, file.path(OUT_SS, "ssGSEA_Hallmark_clusters_pairwise_Wilcoxon.csv"))

# 6b. Grade 2 vs Grade 3: Wilcoxon per pathway (BH over all pathways); difference = first minus second
sl_g <- score_long[!is.na(Grade)]
gw <- sl_g[, {
  a <- Score[Grade == GRADE_PAIR[1]]; b <- Score[Grade == GRADE_PAIR[2]]
  .(Median_first = median(a), Median_second = median(b), Median_diff = median(a) - median(b),
    Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
}, by = Pathway]
setnames(gw, c("Median_first", "Median_second"), paste0("Median_", sub(" ", "_", GRADE_PAIR)))
gw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
gw[, Higher_in := fifelse(Median_diff > 0, GRADE_PAIR[1], fifelse(Median_diff < 0, GRADE_PAIR[2], "tie"))]
setorder(gw, Wilcoxon_FDR, Wilcoxon_p)
fwrite(gw, file.path(OUT_SS, "ssGSEA_Hallmark_grade_Wilcoxon.csv"))
message("ssGSEA significant pathways (FDR < ", FDR_CUT, "): clusters (Kruskal-Wallis) = ",
        sum(kw$Kruskal_FDR < FDR_CUT), "; ", GRADE_PAIR[1], " vs ", GRADE_PAIR[2], " (Wilcoxon) = ",
        sum(gw$Wilcoxon_FDR < FDR_CUT))

# 6c. Heatmap of all Hallmark pathways (row z-score, samples ordered by cluster then grade)
ord <- order(metadata$Cluster, metadata$Grade)
z <- t(scale(t(ss_scores[, ord, drop = FALSE])))
z[z >  2.5] <-  2.5
z[z < -2.5] <- -2.5
ann <- data.frame(Cluster = metadata$Cluster, Grade = metadata$Grade, row.names = rownames(metadata))
draw_ss_heatmap <- function() {
  pheatmap(z, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           annotation_col = ann, annotation_colors = list(Cluster = cluster_colors, Grade = grade_colors),
           cluster_rows = TRUE, cluster_cols = FALSE,
           gaps_col = head(cumsum(table(metadata$Cluster)), -1),
           fontsize_row = 8, fontsize_col = 7, border_color = NA,
           main = "ssGSEA Hallmark scores (row z-score; samples ordered by cluster, then grade)")
}
save_pheatmap(draw_ss_heatmap, file.path(OUT_SS, "ssGSEA_Hallmark_heatmap"), 10, 9)

# 6d. Boxplots of significant pathways: one facet figure + one individual plot per pathway
boxplots_for <- function(sl, stats_dt, fdr_col, group_var, group_colors, group_levels, prefix) {
  sig <- stats_dt[get(fdr_col) < FDR_CUT][order(get(fdr_col))]
  if (nrow(sig) == 0) { message("No significant pathway for ", prefix, " (", fdr_col, " < ", FDR_CUT, ")."); return(invisible(NULL)) }
  ind_dir <- file.path(OUT_SS, paste0(prefix, "_individual"))
  dir.create(ind_dir, recursive = TRUE, showWarnings = FALSE)
  d <- sl[Pathway %in% sig$Pathway]
  d[, (group_var) := factor(get(group_var), levels = group_levels)]
  box_layers <- list(
    geom_boxplot(outlier.shape = NA, alpha = 0.4, colour = "black", linewidth = 0.8),
    geom_jitter(width = 0.15, size = 2, colour = "black", alpha = 0.8),
    scale_fill_manual(values = group_colors))

  for (i in seq_len(nrow(sig))) {                       # individual plots
    pth <- sig$Pathway[i]; fdr <- sig[[fdr_col]][i]
    p <- ggplot(d[Pathway == pth], aes(.data[[group_var]], Score, fill = .data[[group_var]])) +
      box_layers +
      labs(title = paste0("HALLMARK_", pth, "\nFDR-adjusted p = ", fmt_p(fdr)),
           x = NULL, y = "ssGSEA Enrichment Score") +
      theme_classic(base_size = 16) +
      theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            axis.title = element_text(face = "bold"))
    save_plot(p, file.path(ind_dir, paste0("ssGSEA_", pth)), 6.5, 5.5)
  }

  lab <- setNames(paste0(sig$Pathway, "\nFDR = ", vapply(sig[[fdr_col]], fmt_p, character(1))), sig$Pathway)
  d[, Panel := factor(lab[Pathway], levels = unname(lab))]
  pf <- ggplot(d, aes(.data[[group_var]], Score, fill = .data[[group_var]])) +
    box_layers +
    facet_wrap(~ Panel, scales = "free_y", ncol = 3) +
    labs(x = NULL, y = "ssGSEA Enrichment Score",
         title = paste0("Hallmark pathways with FDR < ", FDR_CUT, " (", prefix, ")")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 9), panel.grid.minor = element_blank())
  save_plot(pf, file.path(OUT_SS, paste0(prefix, "_facet")), 12, 3.3 * ceiling(nrow(sig) / 3) + 0.8)
  invisible(sig$Pathway)
}

boxplots_for(score_long, kw, "Kruskal_FDR", "Cluster", cluster_colors, cluster_levels,
             "ssGSEA_boxplots_clusters")
boxplots_for(sl_g, gw, "Wilcoxon_FDR", "Grade", grade_colors, GRADE_PAIR,
             "ssGSEA_boxplots_grade")

message("Done. Outputs in: ", OUT_ROOT)
