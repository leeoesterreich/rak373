###############################################################################
# pILC RNA-seq | script: ssGSEA with your own GMT file (e.g. tams.gmt) across the NMF clusters (Ref: https://doi.org/10.1038/s41698-024-00746-z)
# Reads the outputs of pILC_01a_DESeq2_clusters.R (run it first).
#   Scores : ssGSEA (GSVA ssgseaParam, normalize = TRUE) on the cluster-blind VST matrix,
#            every gene set in the GMT file(s) is scored (no name filter)
#   Tests  : Kruskal-Wallis across clusters + pairwise Wilcoxon (BH-FDR), Grade 2 vs 3 (Wilcoxon),
#            Spearman correlation with tumor cellularity / TILs / pleomorphic component
#   Plots  : score heatmap, boxplots (facet + individual), correlation tiles
# Note   : clusters were derived from the same expression data, so p-values describe how the
#          groups differ and are not independent validation. Scores of immune or stromal signatures
#          also reflect how much of that cell type is in the sample, not only tumor-cell biology.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
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
META_FILE <- file.path(SHARED, "Sample_metadata.tsv")
VST_FILE  <- file.path(SHARED, "VST_blind_ssGSEA_input.tsv")      # gene symbols x samples (cluster-blind VST)

GMT_FILES <- file.path(ROOT, "Pathway_analysis", "Dataset", "tams.gmt")   # one or more GMT files
GMT_TAG   <- tools::file_path_sans_ext(basename(GMT_FILES[1]))             # used in the folder and file names

MIN_SET_GENES <- 5      # a set needs at least this many genes present in the expression matrix
MAX_SET_GENES <- 500
FDR_CUT       <- 0.05

OUT_ROOT <- file.path(ROOT, "Pathway_analysis", paste0("ssGSEA_", GMT_TAG))   # no dates in names
OUT_IND  <- file.path(OUT_ROOT, "boxplots_individual")
dir.create(OUT_IND, recursive = TRUE, showWarnings = FALSE)

for (f in c(META_FILE, VST_FILE)) if (!file.exists(f)) stop("Missing input (run pILC_01a first): ", f)
if (any(!file.exists(GMT_FILES))) stop("GMT file not found: ", paste(GMT_FILES[!file.exists(GMT_FILES)], collapse = ", "))
if (packageVersion("GSVA") < "1.50.0") stop("GSVA >= 1.50 is needed (ssgseaParam).")

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
clin_vars      <- c("Tumor_cellularity", "TILs", "Pleomorphic_component")
GRADE_PAIR     <- c("Grade 2", "Grade 3")

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
pre <- function(x) file.path(OUT_ROOT, paste0(GMT_TAG, "_", x))      # output file prefix

# ------------------------------- 3. LOAD DATA --------------------------------
meta_dt <- fread(META_FILE)
need <- c("Sample", "Cluster", "Grade")
if (length(setdiff(need, names(meta_dt)))) stop("Sample_metadata.tsv lacks: ", paste(setdiff(need, names(meta_dt)), collapse = ", "))
metadata <- as.data.frame(meta_dt)
rownames(metadata) <- metadata$Sample
metadata$Cluster <- factor(metadata$Cluster, levels = cluster_levels)
clin_vars <- intersect(clin_vars, names(metadata))
message("Cluster sizes:"); print(table(metadata$Cluster))

vst_dt <- fread(VST_FILE)
expr <- as.matrix(vst_dt[, names(vst_dt)[-1], with = FALSE])
rownames(expr) <- vst_dt[[1]]
if (!all(rownames(metadata) %in% colnames(expr))) stop("Samples of Sample_metadata.tsv are missing in the VST matrix.")
expr <- expr[, rownames(metadata)]
message(sprintf("Expression matrix: %d genes x %d samples", nrow(expr), ncol(expr)))

# ------------------------------- 4. ssGSEA -----------------------------------
sets <- list()
for (f in GMT_FILES) {
  gmt <- strsplit(readLines(f), "\t")
  sets <- c(sets, setNames(lapply(gmt, function(x) unique(x[-(1:2)])), vapply(gmt, `[`, character(1), 1)))
}
set_size <- data.table(Gene_set = names(sets), Genes_in_set = lengths(sets),
                       Genes_in_matrix = vapply(sets, function(g) sum(g %in% rownames(expr)), integer(1)))
set_size[, Used := Genes_in_matrix >= MIN_SET_GENES & Genes_in_matrix <= MAX_SET_GENES]
fwrite(set_size, pre("gene_set_sizes.tsv"), sep = "\t", quote = FALSE)
print(set_size)
sets <- lapply(sets[set_size$Used], function(g) intersect(g, rownames(expr)))
if (length(sets) == 0) stop("No gene set has between ", MIN_SET_GENES, " and ", MAX_SET_GENES, " genes in the expression matrix.")

all_scores <- gsva(ssgseaParam(exprData = expr, geneSets = sets, normalize = TRUE))   # sets x samples
score_names <- rownames(all_scores)

sample_tab <- data.table(Sample = colnames(all_scores), t(all_scores))
sample_tab <- merge(meta_dt[, c("Sample", "Cluster", "Grade", clin_vars), with = FALSE], sample_tab, by = "Sample")
fwrite(sample_tab, pre("ssGSEA_scores_per_sample.tsv"), sep = "\t", quote = FALSE, na = "")

score_long <- melt(as.data.table(t(all_scores), keep.rownames = "Sample"),
                   id.vars = "Sample", variable.name = "Score", value.name = "Value")
score_long[, Score   := as.character(Score)]
score_long[, Cluster := metadata[Sample, "Cluster"]]
score_long[, Grade   := metadata[Sample, "Grade"]]

# ------------------------------- 5. STATISTICS -------------------------------
# 5a. Kruskal-Wallis across the clusters (BH over all sets) + median per cluster
kw <- score_long[, .(Kruskal_p = kruskal.test(Value ~ Cluster)$p.value), by = Score]
kw[, Kruskal_FDR := p.adjust(Kruskal_p, method = "BH")]
med <- dcast(score_long[, .(Median = median(Value)), by = .(Score, Cluster)], Score ~ Cluster, value.var = "Median")
med_cols <- paste0("Median_", sub(" ", "_", cluster_levels))
setnames(med, cluster_levels, med_cols)
kw <- merge(kw, med, by = "Score")
kw[, Highest_cluster := cluster_levels[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = med_cols]
setorder(kw, Kruskal_FDR, Kruskal_p)
fwrite(kw, pre("clusters_KruskalWallis.tsv"), sep = "\t", quote = FALSE)
print(kw)

# 5b. Pairwise Wilcoxon (difference = first minus second; BH over all set x pair tests)
pw <- rbindlist(lapply(combn(cluster_levels, 2, simplify = FALSE), function(pr) {
  score_long[Cluster %in% pr, {
    a <- Value[Cluster == pr[1]]; b <- Value[Cluster == pr[2]]
    .(Comparison = paste0(sub(" ", "", pr[1]), "_vs_", sub(" ", "", pr[2])),
      Median_diff = median(a) - median(b),
      Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
  }, by = Score]
}))
pw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
setorder(pw, Wilcoxon_FDR, Wilcoxon_p)
fwrite(pw, pre("clusters_pairwise_Wilcoxon.tsv"), sep = "\t", quote = FALSE)

# 5c. Grade 2 vs Grade 3 (Wilcoxon; difference = first minus second)
sl_g <- score_long[Grade %in% GRADE_PAIR]
if (all(table(factor(sl_g[Score == score_names[1], Grade], levels = GRADE_PAIR)) >= 3)) {
  gw <- sl_g[, {
    a <- Value[Grade == GRADE_PAIR[1]]; b <- Value[Grade == GRADE_PAIR[2]]
    .(Median_diff = median(a) - median(b), Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
  }, by = Score]
  gw[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
  gw[, Higher_in := fifelse(Median_diff > 0, GRADE_PAIR[1], fifelse(Median_diff < 0, GRADE_PAIR[2], "tie"))]
  setorder(gw, Wilcoxon_FDR, Wilcoxon_p)
  fwrite(gw, pre("grade_Wilcoxon.tsv"), sep = "\t", quote = FALSE)
} else message("Fewer than 3 samples in a grade group: grade test skipped.")

# 5d. Spearman correlation with tumor cellularity, TILs and pleomorphic component (BH over all tests)
cor_tab <- rbindlist(lapply(score_names, function(s) rbindlist(lapply(clin_vars, function(v) {
  x <- all_scores[s, rownames(metadata)]; y <- as.numeric(metadata[[v]])
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.table(Score = s, Variable = v, n = sum(ok), Spearman_rho = unname(ct$estimate), p = ct$p.value)
}))))
if (nrow(cor_tab) > 0) {
  cor_tab[, FDR := p.adjust(p, method = "BH")]
  setorder(cor_tab, FDR, p)
  fwrite(cor_tab, pre("clinical_correlation_Spearman.tsv"), sep = "\t", quote = FALSE)
}

# ------------------------------- 6. FIGURES ----------------------------------
ann_col <- cbind(data.frame(Cluster = metadata$Cluster, row.names = rownames(metadata)),
                 metadata[, clin_vars, drop = FALSE])
ann_col <- ann_col[, colSums(!is.na(ann_col)) > 0, drop = FALSE]
# samples grouped by cluster, ordered by similarity within each cluster
ord <- unlist(lapply(levels(metadata$Cluster), function(g) {
  cols <- rownames(metadata)[metadata$Cluster == g]
  if (length(cols) > 2) cols[hclust(dist(t(all_scores[, cols, drop = FALSE])), "ward.D2")$order] else cols
}))
gaps_col <- as.integer(head(cumsum(table(metadata[ord, "Cluster"])), -1))

# 6a. Score heatmap (each set as a z-score across samples; needs at least 2 sets)
if (nrow(all_scores) >= 2) {
  sm <- t(scale(t(all_scores[, ord, drop = FALSE])))
  sm[sm >  2.5] <-  2.5
  sm[sm < -2.5] <- -2.5
  draw_scores <- function() {
    pheatmap(sm, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
             annotation_col = ann_col[ord, , drop = FALSE], annotation_colors = list(Cluster = cluster_colors),
             cluster_rows = nrow(sm) > 2, cluster_cols = FALSE, gaps_col = gaps_col,
             clustering_method = "ward.D2", fontsize_row = 8, fontsize_col = 7, border_color = NA,
             main = paste0(GMT_TAG, " ssGSEA scores; samples grouped by NMF cluster"))
  }
  save_pheatmap(draw_scores, pre("ssGSEA_score_heatmap"), 10, max(4.5, 0.35 * nrow(sm) + 3))
}

# 6b. Boxplots by cluster: one facet figure + one individual plot per set with Kruskal FDR < FDR_CUT
lab <- kw[, .(Score, Label = paste0(Score, "\nKruskal-Wallis FDR = ", fmt_p(Kruskal_FDR)))]
sl  <- merge(score_long, lab, by = "Score")
sl[, Label := factor(Label, levels = lab$Label)]
box_layer <- list(geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6),
                  geom_jitter(width = 0.12, size = 1.8, shape = 21, fill = "white"),
                  scale_fill_manual(values = cluster_colors),
                  theme_bw(base_size = 11),
                  theme(legend.position = "none", panel.grid.minor = element_blank()))
p_facet <- ggplot(sl, aes(Cluster, Value, fill = Cluster)) + box_layer +
  facet_wrap(~ Label, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "ssGSEA score") +
  theme(strip.text = element_text(size = 8), axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_facet, pre("ssGSEA_by_cluster_facet"), 11, 3.2 * ceiling(length(score_names) / 3) + 0.6)

for (s in kw[Kruskal_FDR < FDR_CUT, Score]) {
  p <- ggplot(score_long[Score == s], aes(Cluster, Value, fill = Cluster)) + box_layer +
    labs(title = s, subtitle = paste0("Kruskal-Wallis FDR = ", fmt_p(kw[Score == s, Kruskal_FDR])),
         x = NULL, y = "ssGSEA score") +
    theme(plot.title = element_text(face = "bold"))
  save_plot(p, file.path(OUT_IND, paste0(GMT_TAG, "_", gsub("[^A-Za-z0-9_]+", "_", s))), 5, 4.5)
}

# 6c. Correlation with clinical variables (* = FDR < FDR_CUT)
if (nrow(cor_tab) > 0) {
  cor_tab[, Score := factor(Score, levels = rev(kw$Score))]
  p_cor <- ggplot(cor_tab, aes(Variable, Score, fill = Spearman_rho)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(FDR < FDR_CUT, "*", "")), size = 6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "Spearman rho", title = "ssGSEA scores vs clinical variables",
         subtitle = paste0("* FDR < ", FDR_CUT, " (BH over all set x variable tests)")) +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_cor, pre("clinical_correlation"), 6.5, max(3.5, 0.3 * length(score_names) + 2))
}

message("Done. Outputs in: ", OUT_ROOT)