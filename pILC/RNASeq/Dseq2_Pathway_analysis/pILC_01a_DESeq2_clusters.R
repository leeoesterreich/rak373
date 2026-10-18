###############################################################################
# pILC RNA-seq | script 1a: DESeq2 across the three R consensus-NMF clusters
#   (global LRT + pairwise Wald). Run this first; the Grade 2 vs Grade 3 analysis is in
#   pILC_01b_DESeq2_grade.R, and GSEA / ORA / ssGSEA is in pILC_02_GSEA_ORA_ssGSEA.R.
# Genes : ALL annotated genes that pass the count prefilter (no biotype or gene-name filter)
# Also writes the sample table (cluster, grade, clinical variables) and a cluster-blind VST
# matrix that the other two scripts read.
# Note  : clusters were derived from the same expression data, so p-values
#         describe how the groups differ and are not independent validation.
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
PILC_ROOT     <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"
ROOT          <- file.path(PILC_ROOT, "RNASeq")
VST_DIR       <- file.path(ROOT, "4_Count_file_HTSeq_stranded", "VST_data_GENCODEv47")
COUNT_FILE    <- file.path(VST_DIR, "htseq_raw_counts_with_gene_names.csv")  # column "Gene" + one column per sample
ANNOT_FILE    <- file.path(VST_DIR, "gene_annotation_from_GTF.tsv")
CLUSTER_FILE  <- file.path(ROOT, "NMF", "R_NMFConsensus_50runs", "Cluster_membership_k3_GPlabels.csv")
CLINICAL_FILE <- file.path(PILC_ROOT, "Clinical_correlation", "pILC_ClinicoPathological_details.csv")

OUT_ROOT <- file.path(ROOT, "DESeq2_stranded")                       # no dates in names
OUT_CL   <- file.path(OUT_ROOT, "DESeq2_Clusters_R_NMF_k3")
OUT_SH   <- file.path(OUT_ROOT, "Shared_inputs_for_pathway_analysis")  # read by script 2
for (d in c(OUT_CL, OUT_SH)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

MIN_COUNT        <- 10    # prefilter: >= MIN_COUNT counts ...
MIN_SAMPLES      <- 3     # ... in >= MIN_SAMPLES samples (computed within each analysis)
PADJ_CUT         <- 0.05  # FDR cutoff
LFC_CUT          <- 1     # |shrunken log2 fold change| cutoff
N_LABEL          <- 10    # genes labelled per side in each volcano plot
MAX_HEATMAP_GENES <- 100  # heatmaps show at most this many genes that pass BOTH cutoffs (use Inf for all)
SSGSEA_MIN_TOTAL <- 10    # light filter (total counts) for the ssGSEA input matrix
PCA_NTOP         <- 2000  # PCA uses the PCA_NTOP most variable genes

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
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

# ------------------------------- 3. CLUSTERS ---------------------------------
cl <- fread(CLUSTER_FILE)
message("Cluster file columns: ", paste(names(cl), collapse = ", "))
print(head(cl, 3))
sample_col  <- names(cl)[1]
cluster_col <- setdiff(grep("cluster|member|group", names(cl), ignore.case = TRUE, value = TRUE),
                       sample_col)[1]
if (is.na(cluster_col)) cluster_col <- names(cl)[2]
message("Using sample column '", sample_col, "' and cluster column '", cluster_col, "'")

cluster_num <- gsub("[^0-9]", "", as.character(cl[[cluster_col]]))   # 1 / C1 / Cluster 1 -> 1
if (any(cluster_num == "")) stop("Some cluster labels contain no number: check CLUSTER_FILE.")
cluster_map <- setNames(paste0("Cluster ", cluster_num), norm_id(cl[[sample_col]]))
if (anyDuplicated(names(cluster_map))) stop("Duplicated sample IDs in CLUSTER_FILE.")
if (!all(cluster_map %in% cluster_levels))
  stop("Unexpected cluster labels: ", paste(unique(cluster_map), collapse = ", "))

message("Cluster sizes (printed for information; no size check):")
print(table(factor(cluster_map, levels = cluster_levels)))

# ------------------------------- 4. CLINICAL (grade and other variables for the sample table)
clin <- fread(CLINICAL_FILE)
if (!("Sample" %in% names(clin))) stop("The clinical file needs a column named 'Sample'.")
clin[, Sample_key := norm_id(Sample)]
if (anyDuplicated(clin$Sample_key)) stop("Duplicated samples in the clinical file.")

grade_col <- grep("^Grade", names(clin), value = TRUE)[1]
if (is.na(grade_col)) stop("No column starting with 'Grade' in the clinical file.")
grade_raw <- trimws(clin[[grade_col]])
clin[, Grade := fifelse(grepl("^Grade[[:space:]]*(II|2)$", grade_raw, ignore.case = TRUE), "Grade 2",
                fifelse(grepl("^Grade[[:space:]]*(III|3)$", grade_raw, ignore.case = TRUE), "Grade 3",
                        NA_character_))]
message("Grade column '", grade_col, "': original entries (rows) vs. grade used here")
print(table(grade_raw, clin$Grade, useNA = "ifany"))

# First column whose name matches the pattern (case-insensitive); a "%" sign is removed before
# converting to a number. The chosen column is printed so you can confirm it is the right one.
num_col <- function(pattern) {
  cc <- grep(pattern, names(clin), value = TRUE, ignore.case = TRUE)[1]
  message("Clinical column for '", pattern, "': ", ifelse(is.na(cc), "NOT FOUND (annotation skipped)", cc))
  if (is.na(cc)) rep(NA_real_, nrow(clin))
  else suppressWarnings(as.numeric(gsub("%", "", trimws(as.character(clin[[cc]])))))
}
clin[, Tumor_cellularity := num_col("^tumor_cellularity")]
clin[, TILs := num_col("^TILs")]
clin[, Pleomorphic_component := num_col("pleomorph")]

m <- match(names(cluster_map), clin$Sample_key)
if (anyNA(m)) message("No clinical row for: ", paste(names(cluster_map)[is.na(m)], collapse = ", "))
metadata <- data.frame(Cluster = factor(unname(cluster_map), levels = cluster_levels),
                       Grade = clin$Grade[m],
                       Grade_description = grade_raw[m],
                       Tumor_cellularity = clin$Tumor_cellularity[m],
                       TILs = clin$TILs[m],
                       Pleomorphic_component = clin$Pleomorphic_component[m],
                       row.names = names(cluster_map))
fwrite(data.table(Sample = rownames(metadata), metadata), file.path(OUT_SH, "Sample_metadata.tsv"), sep = "\t")

# Heatmap annotation: Cluster (and Grade when given) plus the three continuous clinical variables.
# A variable with no values at all is left out, because pheatmap cannot colour an empty column.
clin_vars <- c("Tumor_cellularity", "TILs", "Pleomorphic_component")
make_ann <- function(md, first_cols) {
  a <- md[, clin_vars, drop = FALSE]
  a <- a[, colSums(!is.na(a)) > 0, drop = FALSE]
  cbind(md[, first_cols, drop = FALSE], a)
}

# ------------------------------- 5. COUNTS + GENES ---------------------------
count_dt <- fread(COUNT_FILE)
if (!("Gene" %in% names(count_dt))) stop("Count file must contain a column named 'Gene'.")
if (anyDuplicated(count_dt$Gene)) stop("Duplicate gene labels in the count file.")

sample_cols <- setdiff(names(count_dt), "Gene")
idx <- match(names(cluster_map), norm_id(sample_cols))
if (anyNA(idx))
  stop("Cluster samples not found in counts: ", paste(names(cluster_map)[is.na(idx)], collapse = ", "),
       "\nCount columns look like: ", paste(head(sample_cols, 5), collapse = ", "))
extra <- setdiff(norm_id(sample_cols), names(cluster_map))
if (length(extra)) stop("Count samples without a cluster: ", paste(extra, collapse = ", "))

count_matrix <- as.matrix(count_dt[, sample_cols[idx], with = FALSE])
rownames(count_matrix) <- count_dt$Gene
colnames(count_matrix) <- names(cluster_map)          # same order as metadata
storage.mode(count_matrix) <- "numeric"
if (anyNA(count_matrix)) stop("Missing count values detected.")
if (any(count_matrix < 0)) stop("Negative count values detected.")
if (any(count_matrix != round(count_matrix))) stop("Non-integer counts detected.")
storage.mode(count_matrix) <- "integer"
count_matrix <- count_matrix[!grepl("^__", rownames(count_matrix)), , drop = FALSE]  # HTSeq summary rows

gene_annot <- fread(ANNOT_FILE)[, .(Gene, ensembl_id, gene_symbol, gene_biotype)]
if (anyDuplicated(gene_annot$Gene)) stop("Duplicate Gene labels in the annotation file.")
add_annotation <- function(dt) merge(dt, gene_annot, by = "Gene", all.x = TRUE, sort = FALSE)

# A gene "has a name" when the symbol is not empty and not just an Ensembl ID.
# No gene is removed here: DESeq2 runs on all genes. gene_label is only used for plots.
gene_annot[, has_symbol := !is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSG[0-9]+", gene_symbol)]
gene_annot[, gene_label := fifelse(has_symbol, gene_symbol, Gene)]
message("Genes in count file: ", nrow(count_dt), " | genes used for DESeq2 (all): ", nrow(count_matrix),
        " | with a gene name: ", sum(gene_annot$has_symbol[match(rownames(count_matrix), gene_annot$Gene)], na.rm = TRUE))

# ------------------------------- 6. SHARED INPUT FOR ssGSEA ------------------
# All genes with a light filter; VST blind to cluster/grade. Only genes with a gene symbol are kept
# for ssGSEA, because the Hallmark sets are written as gene symbols.
dds_ss <- DESeqDataSetFromMatrix(countData = count_matrix[rowSums(count_matrix) >= SSGSEA_MIN_TOTAL, ],
                                 colData = metadata[, "Cluster", drop = FALSE], design = ~ 1)
expr <- assay(vst(dds_ss, blind = TRUE))
sym  <- gene_annot$gene_symbol[match(rownames(expr), gene_annot$Gene)]
keep <- gene_annot$has_symbol[match(rownames(expr), gene_annot$Gene)] %in% TRUE
expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
o    <- order(rowMeans(expr), decreasing = TRUE)      # duplicated symbols: keep the highest mean
expr <- expr[o, , drop = FALSE]; sym <- sym[o]
expr <- expr[!duplicated(sym), , drop = FALSE]
rownames(expr) <- sym[!duplicated(sym)]
fwrite(as.data.table(expr, keep.rownames = "Gene_symbol"),
       file.path(OUT_SH, "VST_blind_ssGSEA_input.tsv"), sep = "\t", quote = FALSE)
message("ssGSEA input matrix written: ", nrow(expr), " genes x ", ncol(expr), " samples")

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
# order_by = NULL : samples are clustered freely (dendrogram on top).
# order_by = "Cluster" (a column of ann_col): samples are grouped by that column, ordered by similarity
#   (Ward clustering) within each group, with a gap between groups and no column dendrogram.
plot_top_heatmap <- function(vst_mat, genes, ann_col, ann_colors, main, stem, w, h, order_by = NULL) {
  hm <- t(scale(t(vst_mat[genes, rownames(ann_col), drop = FALSE])))
  hm[hm >  2.5] <-  2.5
  hm[hm < -2.5] <- -2.5
  rownames(hm) <- make.unique(gene_annot$gene_label[match(genes, gene_annot$Gene)])
  hm <- hm[complete.cases(hm), , drop = FALSE]
  gaps <- NULL
  if (!is.null(order_by)) {
    grp <- as.character(ann_col[colnames(hm), order_by])
    ord <- unlist(lapply(unique(sort(grp)), function(g) {
      cols <- colnames(hm)[grp == g]
      if (length(cols) > 2) cols[hclust(dist(t(hm[, cols, drop = FALSE])), "ward.D2")$order] else cols
    }))
    hm   <- hm[, ord, drop = FALSE]
    gaps <- as.integer(head(cumsum(table(factor(as.character(ann_col[ord, order_by]),
                                     levels = unique(as.character(ann_col[ord, order_by]))))), -1))
  }
  draw <- function() {
    pheatmap(hm, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
             annotation_col = ann_col, annotation_colors = ann_colors,
             cluster_rows = TRUE, cluster_cols = is.null(order_by), gaps_col = gaps,
             clustering_method = "ward.D2",
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
# A) DESeq2 across the three clusters
# =============================================================================
dds_cl0 <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = metadata[, "Cluster", drop = FALSE], design = ~ Cluster)
dds_cl  <- dds_cl0[rowSums(counts(dds_cl0) >= MIN_COUNT) >= MIN_SAMPLES, ]
message("[Clusters] prefilter (>= ", MIN_COUNT, " counts in >= ", MIN_SAMPLES, " samples): ",
        nrow(dds_cl), " of ", nrow(dds_cl0), " genes retained")

# Global likelihood-ratio test: does a gene differ in any cluster?
lrt <- as.data.table(as.data.frame(results(DESeq(dds_cl, test = "LRT", reduced = ~ 1, quiet = TRUE))),
                     keep.rownames = "Gene")
lrt[, c("log2FoldChange", "lfcSE") := NULL]    # an LRT shows only one coefficient: removed to avoid misreading
lrt <- add_annotation(lrt)
lrt[, Significant_FDR05 := !is.na(padj) & padj < PADJ_CUT]
setorder(lrt, padj, pvalue, -baseMean, na.last = TRUE)
fwrite(lrt, file.path(OUT_CL, "DESeq2_global_LRT_any_cluster.tsv"), sep = "\t", quote = FALSE, na = "")

# Pairwise Wald tests (positive log2FC = higher in the first-named cluster)
dds_cl_wald <- DESeq(dds_cl, quiet = TRUE)
cl_comparisons <- list(
  Cluster1_vs_Cluster2 = c("Cluster", "Cluster 1", "Cluster 2"),
  Cluster1_vs_Cluster3 = c("Cluster", "Cluster 1", "Cluster 3"),
  Cluster2_vs_Cluster3 = c("Cluster", "Cluster 2", "Cluster 3"))
cl_results <- lapply(names(cl_comparisons), function(nm) run_wald_contrast(dds_cl_wald, cl_comparisons[[nm]], nm))
names(cl_results) <- names(cl_comparisons)
fwrite(rbindlist(cl_results, use.names = TRUE, fill = TRUE),
       file.path(OUT_CL, "DESeq2_all_pairwise_results.tsv"), sep = "\t", quote = FALSE, na = "")

cl_summary <- rbind(
  data.table(Analysis = "Global_LRT_any_cluster", Group_A = "All", Group_B = "All",
             Tested_genes = sum(!is.na(lrt$padj)), FDR05_genes = sum(lrt$Significant_FDR05)),
  rbindlist(lapply(names(cl_comparisons), function(nm)
    summarise_result(cl_results[[nm]], cl_comparisons[[nm]], nm))),
  fill = TRUE)
fwrite(cl_summary, file.path(OUT_CL, "DESeq2_analysis_summary.tsv"), sep = "\t", quote = FALSE, na = "")
print(cl_summary)

for (nm in names(cl_comparisons))
  plot_volcano(cl_results[[nm]], cl_comparisons[[nm]], cluster_colors,
               file.path(OUT_CL, paste0("DESeq2_volcano_", nm)))

# Heatmap genes: global LRT FDR < PADJ_CUT AND |shrunken log2FC| >= LFC_CUT with FDR < PADJ_CUT in at
# least one pairwise comparison. Ranked by LRT FDR; at most MAX_HEATMAP_GENES are drawn (VST: display only).
pass_pair <- unique(unlist(lapply(cl_results, function(x) x$Gene[x$Significant_FDR05_LFC1])))
hm_cl <- lrt[Significant_FDR05 == TRUE & Gene %in% pass_pair]
message("[Cluster heatmap] genes passing both cutoffs: ", nrow(hm_cl))
if (nrow(hm_cl) < 2) {
  warning("Fewer than 2 genes pass both cutoffs: cluster heatmap skipped.")
} else {
  hm_genes <- head(hm_cl$Gene, MAX_HEATMAP_GENES)
  ann_cl <- make_ann(metadata, "Cluster")
  plot_top_heatmap(assay(vst(dds_cl_wald, blind = FALSE)), hm_genes, ann_cl,
                   list(Cluster = cluster_colors),
                   paste0(length(hm_genes), " genes differing between clusters (FDR < ", PADJ_CUT,
                          ", |log2FC| >= ", LFC_CUT, ")"),
                   file.path(OUT_CL, "DESeq2_cluster_DE_genes_heatmap"), 11, 13)
  fwrite(hm_cl[Gene %in% hm_genes], file.path(OUT_CL, "DESeq2_cluster_DE_genes_heatmap_genes.tsv"),
         sep = "\t", quote = FALSE, na = "")

  # Same genes, samples grouped by NMF cluster (display choice: the dendrogram above shows how the
  # samples group on their own; this version shows where each cluster's samples sit within their group)
  plot_top_heatmap(assay(vst(dds_cl_wald, blind = FALSE)), hm_genes, ann_cl,
                   list(Cluster = cluster_colors),
                   paste0(length(hm_genes), " genes differing between clusters (FDR < ", PADJ_CUT,
                          ", |log2FC| >= ", LFC_CUT, ")\nSamples grouped by NMF cluster"),
                   file.path(OUT_CL, "DESeq2_cluster_DE_genes_heatmap_ordered_by_cluster"), 11, 13,
                   order_by = "Cluster")
}

# PCA (VST blind to the clusters): does each sample sit with its cluster?
plot_pca(dds_cl_wald, "Cluster", NULL, cluster_colors,
         paste0("Top ", PCA_NTOP, " most variable genes, variance-stabilised expression"),
         file.path(OUT_CL, "DESeq2_PCA_clusters"))

message("Done. Outputs in: ", OUT_ROOT)