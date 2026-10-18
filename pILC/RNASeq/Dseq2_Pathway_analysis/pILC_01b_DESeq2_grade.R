###############################################################################
# pILC RNA-seq | script 1b: DESeq2 Grade 2 vs Grade 3 (Wald)
# Run pILC_01a_DESeq2_clusters.R first: it writes Sample_metadata.tsv (cluster, grade,
# tumor cellularity, TILs, pleomorphic component) that this script reads.
# Genes : ALL annotated genes that pass the count prefilter (no biotype or gene-name filter)
# Model : ~ Grade only (no covariates); the grade x cluster table checks for confounding.
# Positive log2FC = higher in the first-named grade (Grade 2).
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
PILC_ROOT  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"
ROOT       <- file.path(PILC_ROOT, "RNASeq")
VST_DIR    <- file.path(ROOT, "4_Count_file_HTSeq_stranded", "VST_data_GENCODEv47")
COUNT_FILE <- file.path(VST_DIR, "htseq_raw_counts_with_gene_names.csv")  # column "Gene" + one column per sample
ANNOT_FILE <- file.path(VST_DIR, "gene_annotation_from_GTF.tsv")

OUT_ROOT  <- file.path(ROOT, "DESeq2_stranded")                       # no dates in names
META_FILE <- file.path(OUT_ROOT, "Shared_inputs_for_pathway_analysis", "Sample_metadata.tsv")  # from script 1a
OUT_GR    <- file.path(OUT_ROOT, "DESeq2_Grade2_vs_Grade3")
dir.create(OUT_GR, recursive = TRUE, showWarnings = FALSE)

MIN_COUNT         <- 10    # prefilter: >= MIN_COUNT counts ...
MIN_SAMPLES       <- 3     # ... in >= MIN_SAMPLES samples (within the grade samples)
PADJ_CUT          <- 0.05  # FDR cutoff
LFC_CUT           <- 1     # |shrunken log2 fold change| cutoff
N_LABEL           <- 10    # genes labelled per side in the volcano plot
MAX_HEATMAP_GENES <- 100   # heatmap shows at most this many genes that pass BOTH cutoffs (use Inf for all)
PCA_NTOP          <- 2000  # PCA uses the PCA_NTOP most variable genes 

# The FIRST name is the numerator (positive log2FC = higher in the first)
GRADE_PAIR <- c("Grade 2", "Grade 3")

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
# Grade colours: green / purple like the earlier grade figure (shades differ from the cluster colours)
grade_colors   <- c("Grade 2" = "#8BC34A", "Grade 3" = "#6A1B9A")
cluster_levels <- names(cluster_colors)

# ------------------------------- 2. HELPERS ----------------------------------
# "TP13.1150_S40" / "TP13-1150" / "TP09-09" -> "TP13.1150" / "TP09.9" (same style as count columns)
norm_id <- function(x) {
  x <- gsub("-", ".", sub("_S[0-9]+$", "", trimws(as.character(x))))
  sub("\\.0+([0-9])", ".\\1", x)
}
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

# ------------------------------- 3. SAMPLE TABLE (from script 1a) -------------
if (!file.exists(META_FILE)) stop("Run pILC_01a_DESeq2_clusters.R first: not found: ", META_FILE)
meta_dt <- fread(META_FILE)
need <- c("Sample", "Cluster", "Grade", "Grade_description", "Tumor_cellularity", "TILs", "Pleomorphic_component")
if (length(setdiff(need, names(meta_dt))))
  stop("Sample_metadata.tsv lacks: ", paste(setdiff(need, names(meta_dt)), collapse = ", "),
       ". Re-run pILC_01a_DESeq2_clusters.R.")
metadata <- as.data.frame(meta_dt[, ..need])
rownames(metadata) <- metadata$Sample
metadata$Sample  <- NULL
metadata$Cluster <- factor(metadata$Cluster, levels = cluster_levels)

# Heatmap annotation: Cluster (and Grade when given) plus the three continuous clinical variables.
# A variable with no values at all is left out, because pheatmap cannot colour an empty column.
clin_vars <- c("Tumor_cellularity", "TILs", "Pleomorphic_component")
make_ann <- function(md, first_cols) {
  a <- md[, clin_vars, drop = FALSE]
  a <- a[, colSums(!is.na(a)) > 0, drop = FALSE]
  cbind(md[, first_cols, drop = FALSE], a)
}

# ------------------------------- 4. COUNTS + GENES ---------------------------
g_samples <- rownames(metadata)[metadata$Grade %in% GRADE_PAIR]

count_dt <- fread(COUNT_FILE)
if (!("Gene" %in% names(count_dt))) stop("Count file must contain a column named 'Gene'.")
if (anyDuplicated(count_dt$Gene)) stop("Duplicate gene labels in the count file.")
sample_cols <- setdiff(names(count_dt), "Gene")
idx <- match(g_samples, norm_id(sample_cols))
if (anyNA(idx))
  stop("Grade samples not found in counts: ", paste(g_samples[is.na(idx)], collapse = ", "),
       "\nCount columns look like: ", paste(head(sample_cols, 5), collapse = ", "))

count_matrix <- as.matrix(count_dt[, sample_cols[idx], with = FALSE])
rownames(count_matrix) <- count_dt$Gene
colnames(count_matrix) <- g_samples
storage.mode(count_matrix) <- "numeric"
if (anyNA(count_matrix)) stop("Missing count values detected.")
if (any(count_matrix < 0)) stop("Negative count values detected.")
if (any(count_matrix != round(count_matrix))) stop("Non-integer counts detected.")
storage.mode(count_matrix) <- "integer"
count_matrix <- count_matrix[!grepl("^__", rownames(count_matrix)), , drop = FALSE]  # HTSeq summary rows

gene_annot <- fread(ANNOT_FILE)[, .(Gene, ensembl_id, gene_symbol, gene_biotype)]
if (anyDuplicated(gene_annot$Gene)) stop("Duplicate Gene labels in the annotation file.")
add_annotation <- function(dt) merge(dt, gene_annot, by = "Gene", all.x = TRUE, sort = FALSE)
# A gene "has a name" when the symbol is not empty and not just an Ensembl ID; gene_label is for plots only.
gene_annot[, has_symbol := !is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSG[0-9]+", gene_symbol)]
gene_annot[, gene_label := fifelse(has_symbol, gene_symbol, Gene)]

# ------------------------------- 7. FUNCTIONS FOR DESeq2 RESULTS -------------
# One Wald contrast: shrunken log2FC, unshrunken Wald statistic (for GSEA ranking), flags
run_wald_contrast <- function(dds_wald, cv, cmp_name) {
  raw <- results(dds_wald, contrast = cv, alpha = PADJ_CUT, independentFiltering = TRUE)
  wald_stat <- setNames(raw$stat, rownames(raw))
  lfc_type  <- "shrunken_normal"
  res <- tryCatch(
    lfcShrink(dds_wald, contrast = cv, res = raw, type = "normal"),
    error = function(e) {
      warning("LFC shrinkage failed for ", cmp_name, ": ", conditionMessage(e))
      lfc_type <<- "unshrunken"
      raw
    })
  d <- as.data.table(as.data.frame(res), keep.rownames = "Gene")
  if ("stat" %in% names(d)) d[, stat := NULL]
  d[, Wald_stat := wald_stat[Gene]]
  d[, Comparison := cmp_name]
  d[, LFC_type := lfc_type]
  d[, Significant_FDR05 := !is.na(padj) & padj < PADJ_CUT]
  d[, Significant_FDR05_LFC1 := Significant_FDR05 & abs(log2FoldChange) >= LFC_CUT]
  d[, Direction := fifelse(Significant_FDR05_LFC1 & log2FoldChange > 0, paste0("Higher_in_", cv[2]),
                   fifelse(Significant_FDR05_LFC1 & log2FoldChange < 0, paste0("Higher_in_", cv[3]),
                           "Not_significant"))]
  d <- add_annotation(d)
  d[, abs_lfc := abs(log2FoldChange)]
  setorder(d, padj, pvalue, -abs_lfc, na.last = TRUE)
  d[, abs_lfc := NULL]
  d
}

summarise_result <- function(x, cv, cmp_name) {
  data.table(Analysis = cmp_name, Group_A = cv[2], Group_B = cv[3],
             Tested_genes = sum(!is.na(x$padj)),
             FDR05_genes = sum(x$Significant_FDR05),
             FDR05_LFC1_genes = sum(x$Significant_FDR05_LFC1),
             Higher_in_Group_A = sum(x$Significant_FDR05_LFC1 & x$log2FoldChange > 0),
             Higher_in_Group_B = sum(x$Significant_FDR05_LFC1 & x$log2FoldChange < 0))
}

plot_volcano <- function(x, cv, group_colors, stem) {
  x <- copy(x[!is.na(padj)])
  x[, minus_log10_padj := -log10(pmax(padj, .Machine$double.xmin))]
  g1 <- cv[2]; g2 <- cv[3]
  up_lab <- paste0("Higher in ", g1); down_lab <- paste0("Higher in ", g2)
  x[, Group := "Not significant"]
  x[Significant_FDR05_LFC1 & log2FoldChange > 0, Group := up_lab]
  x[Significant_FDR05_LFC1 & log2FoldChange < 0, Group := down_lab]
  cand <- x[Group != "Not significant"]
  lab  <- rbind(cand[log2FoldChange > 0][order(padj)][seq_len(min(N_LABEL, .N))],
                cand[log2FoldChange < 0][order(padj)][seq_len(min(N_LABEL, .N))])
  cols <- c("Not significant" = "#BDBDBD",
            setNames(group_colors[g1], up_lab), setNames(group_colors[g2], down_lab))
  p <- ggplot(x, aes(log2FoldChange, minus_log10_padj, color = Group)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", color = "grey35") +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", color = "grey35") +
    geom_text_repel(data = lab, aes(label = gene_label), color = "black", size = 3,
                    max.overlaps = Inf, box.padding = 0.35, min.segment.length = 0,
                    segment.color = "grey50", segment.size = 0.3, show.legend = FALSE) +
    scale_color_manual(values = cols, breaks = c(up_lab, down_lab, "Not significant")) +
    labs(title = paste0(g1, " vs ", g2),
         subtitle = paste0(nrow(x), " genes tested; ", nrow(cand),
                           " with |shrunken log2FC| >= ", LFC_CUT, " and adjusted p (FDR) < ", PADJ_CUT,
                           "; top ", N_LABEL, " per side labelled"),
         x = paste0("Log2 fold change (", g1, " / ", g2, ")"),
         y = "-log10 adjusted p-value", color = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title.position = "plot", plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  save_plot(p, stem, 9, 7.5)
}

# Heatmap of selected genes (VST, row z-score clipped at +/-2.5); samples = rows of ann_col
plot_top_heatmap <- function(vst_mat, genes, ann_col, ann_colors, main, stem, w, h) {
  hm <- t(scale(t(vst_mat[genes, rownames(ann_col), drop = FALSE])))
  hm[hm >  2.5] <-  2.5
  hm[hm < -2.5] <- -2.5
  rownames(hm) <- make.unique(gene_annot$gene_label[match(genes, gene_annot$Gene)])
  hm <- hm[complete.cases(hm), , drop = FALSE]
  draw <- function() {
    pheatmap(hm, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
             annotation_col = ann_col, annotation_colors = ann_colors,
             cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "ward.D2",
             fontsize_row = 5, fontsize_col = 7, border_color = NA, main = main)
  }
  save_pheatmap(draw, stem, w, h)
}

plot_pca <- function(dds, color_var, shape_var = NULL, colors, subtitle, stem) {
  pca <- plotPCA(vst(dds, blind = TRUE), intgroup = c(color_var, shape_var), ntop = PCA_NTOP, returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"), 1)
  p <- if (is.null(shape_var)) {
    ggplot(pca, aes(PC1, PC2, color = .data[[color_var]], label = name))
  } else {
    ggplot(pca, aes(PC1, PC2, color = .data[[color_var]], shape = .data[[shape_var]], label = name))
  }
  p <- p + geom_point(size = 4) +
    geom_text_repel(size = 2.7, show.legend = FALSE, max.overlaps = Inf) +
    scale_color_manual(values = colors, drop = FALSE) +
    labs(title = "pILC RNA-seq PCA", subtitle = subtitle,
         x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  save_plot(p, stem, 9, 7)
}

# =============================================================================
# DESeq2 Grade 2 vs Grade 3
# =============================================================================
g_samples <- rownames(metadata)[metadata$Grade %in% GRADE_PAIR]
gm <- metadata[g_samples, ]
gm$Grade <- factor(gm$Grade, levels = GRADE_PAIR)
n_g <- table(gm$Grade)
message("[Grade] samples used: ", paste(names(n_g), n_g, sep = " = ", collapse = ", "))
excluded <- metadata[!rownames(metadata) %in% g_samples, c("Cluster", "Grade_description")]
if (nrow(excluded)) { message("[Grade] excluded samples (other or missing grade):"); print(excluded) }
if (any(n_g < 3)) stop("Fewer than 3 samples in a grade group: the grade comparison is not possible.")

# Confounding check: does grade line up with cluster or tumor cellularity?
xt <- table(Grade = gm$Grade, Cluster = gm$Cluster)
message("[Grade] grade x cluster table:"); print(xt)
fisher_p <- tryCatch(fisher.test(xt)$p.value, error = function(e) NA_real_)
message("[Grade] grade x cluster association, Fisher exact p = ", signif(fisher_p, 3))
cell_med <- tapply(gm$Tumor_cellularity, gm$Grade, median, na.rm = TRUE)
cell_p   <- tryCatch(suppressWarnings(wilcox.test(Tumor_cellularity ~ Grade, data = gm)$p.value),
                     error = function(e) NA_real_)
message("[Grade] median tumor cellularity (%): ",
        paste(names(cell_med), cell_med, sep = " = ", collapse = ", "),
        "; Wilcoxon p = ", signif(cell_p, 3))
fwrite(as.data.table(xt), file.path(OUT_GR, "Grade_vs_Cluster_crosstab.tsv"), sep = "\t")
fwrite(data.table(Sample = rownames(gm), gm), file.path(OUT_GR, "Grade_analysis_samples.tsv"), sep = "\t")

# DESeq2 on the grade subset only (size factors, dispersions and prefilter from these samples)
dds_g0 <- DESeqDataSetFromMatrix(countData = count_matrix[, g_samples, drop = FALSE],
                                 colData = gm[, c("Grade", "Cluster")], design = ~ Grade)
dds_g  <- dds_g0[rowSums(counts(dds_g0) >= MIN_COUNT) >= MIN_SAMPLES, ]
message("[Grade] prefilter (>= ", MIN_COUNT, " counts in >= ", MIN_SAMPLES, " samples): ",
        nrow(dds_g), " of ", nrow(dds_g0), " genes retained")
dds_g_wald <- DESeq(dds_g, quiet = TRUE)

gr_name <- paste0(gsub(" ", "", GRADE_PAIR[1]), "_vs_", gsub(" ", "", GRADE_PAIR[2]))
gr_cv   <- c("Grade", GRADE_PAIR[1], GRADE_PAIR[2])
gr_res  <- run_wald_contrast(dds_g_wald, gr_cv, gr_name)
fwrite(gr_res, file.path(OUT_GR, paste0("DESeq2_", gr_name, "_results.tsv")),
       sep = "\t", quote = FALSE, na = "")
gr_summary <- summarise_result(gr_res, gr_cv, gr_name)
fwrite(gr_summary, file.path(OUT_GR, "DESeq2_analysis_summary.tsv"), sep = "\t", quote = FALSE, na = "")
print(gr_summary)

plot_volcano(gr_res, gr_cv, grade_colors, file.path(OUT_GR, paste0("DESeq2_volcano_", gr_name)))

# Heatmap genes: FDR < PADJ_CUT and |shrunken log2FC| >= LFC_CUT, ranked by FDR (at most MAX_HEATMAP_GENES)
hm_gr <- gr_res[Significant_FDR05_LFC1 == TRUE]
message("[Grade heatmap] genes passing both cutoffs: ", nrow(hm_gr))
if (nrow(hm_gr) < 2) {
  warning("Fewer than 2 genes pass both cutoffs: grade heatmap skipped.")
} else {
  hm_genes <- head(hm_gr$Gene, MAX_HEATMAP_GENES)
  ann_g <- make_ann(gm, c("Grade", "Cluster"))
  plot_top_heatmap(assay(vst(dds_g_wald, blind = FALSE)), hm_genes, ann_g,
                   list(Grade = grade_colors, Cluster = cluster_colors),
                   paste0(length(hm_genes), " genes: FDR < ", PADJ_CUT, " and |log2FC| >= ", LFC_CUT,
                          " (", GRADE_PAIR[1], " vs ", GRADE_PAIR[2], ")"),
                   file.path(OUT_GR, paste0("DESeq2_", gr_name, "_DE_genes_heatmap")), 10, 10)
}

plot_pca(dds_g_wald, "Grade", "Cluster", grade_colors,
         paste0("Grade 2 and Grade 3 samples only; top ", PCA_NTOP, " most variable genes; colour = grade, shape = cluster"),
         file.path(OUT_GR, "DESeq2_PCA_grade"))

message("Done. Outputs in: ", OUT_GR)