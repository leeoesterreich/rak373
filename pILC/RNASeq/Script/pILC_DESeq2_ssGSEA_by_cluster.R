###############################################################################
# pILC RNA-seq: (1) DESeq2 across clusters (protein-coding genes with symbol)
#               (2) ssGSEA (Hallmark) across clusters
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(tidyr)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. PATHS ------------------------------------
ROOT         <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
COUNTS_DIR   <- file.path(ROOT, "4_Count_file_pILC")  # one raw-count file per sample (no header: ENSG ID, count)
VST_FILE     <- file.path(ROOT, "NMF", "vst_transformed_data.csv")
CLUSTER_FILE <- file.path(ROOT, "NMF", "GenePattern_NMFConsensus",
                          "vst_transformed_data.consensus.k.3.gct")
HALLMARK_GMT <- file.path(ROOT, "GSEA", "Dataset", "h.all.v2025.1.Hs.symbols.gmt")

OUT_DIR <- file.path(ROOT, "Downstream", "Cluster_DESeq2_ssGSEA")  # no dates in names
dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

# --------------------------- 2. CLUSTER ANNOTATION ---------------------------
# GCT format: line 1 = "#1.2", line 2 = dimensions, line 3 = header
clusters <- read.delim(CLUSTER_FILE, skip = 2, check.names = FALSE,
                       stringsAsFactors = FALSE)
clusters <- data.frame(sample  = clusters$Name,
                       cluster = clusters$membership.ordered,
                       stringsAsFactors = FALSE)
clusters$cluster <- factor(paste0("C", clusters$cluster),
                           levels = sort(unique(paste0("C", clusters$cluster))))
rownames(clusters) <- clusters$sample
print(table(clusters$cluster))  # check cluster sizes (power depends on this)

# ------------------------- 3. RAW COUNTS + GENE FILTER -----------------------
# Read the per-sample files (2 columns, no header: gene ID, raw count)
count_files <- list.files(COUNTS_DIR, pattern = "\\.txt$", full.names = TRUE)

# Harmonize IDs: file "TP03-117.txt" <-> cluster ID "TP03.117_S25"
file_key    <- sub("\\.txt$", "", basename(count_files))
cluster_key <- gsub("\\.", "-", sub("_S[0-9]+$", "", clusters$sample))

missing <- setdiff(cluster_key, file_key)
if (length(missing) > 0)
  stop("No count file found for: ", paste(missing, collapse = ", "))
count_files <- count_files[match(cluster_key, file_key)]  # same order as metadata

read_one <- function(f) {
  x <- read.delim(f, header = FALSE, col.names = c("gene", "count"),
                  stringsAsFactors = FALSE)
  x[grepl("^ENSG", x$gene), ]          # drop any summary lines (e.g. __no_feature)
}
count_list <- lapply(count_files, read_one)

# All files must contain the same genes in the same order
stopifnot(all(sapply(count_list, function(x) identical(x$gene, count_list[[1]]$gene))))

counts <- sapply(count_list, function(x) x$count)
rownames(counts) <- count_list[[1]]$gene
colnames(counts) <- clusters$sample
mode(counts) <- "integer"
message("Count matrix: ", nrow(counts), " genes x ", ncol(counts), " samples")

# Gene annotation: detect whether row names are Ensembl IDs or symbols
if (all(grepl("^ENSG", rownames(counts)))) {
  ids <- sub("\\..*$", "", rownames(counts))             # remove version suffix
  anno <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ids),
                                keytype = "ENSEMBL",
                                columns = c("SYMBOL", "GENETYPE"))
  anno <- anno[!duplicated(anno$ENSEMBL), ]
  symbol <- anno$SYMBOL[match(ids, anno$ENSEMBL)]
  gtype  <- anno$GENETYPE[match(ids, anno$ENSEMBL)]
} else {
  anno <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(rownames(counts)),
                                keytype = "SYMBOL", columns = "GENETYPE")
  anno <- anno[!duplicated(anno$SYMBOL), ]
  symbol <- rownames(counts)
  gtype  <- anno$GENETYPE[match(symbol, anno$SYMBOL)]
}

# Keep protein-coding genes that have a gene symbol
keep <- !is.na(symbol) & symbol != "" & gtype %in% "protein-coding"
counts <- counts[keep, ]
symbol <- symbol[keep]

# Collapse duplicated symbols: keep the row with the highest mean count
ord    <- order(rowMeans(counts), decreasing = TRUE)
counts <- counts[ord, ]
symbol <- symbol[ord]
counts <- counts[!duplicated(symbol), ]
rownames(counts) <- symbol[!duplicated(symbol)]
message("Protein-coding genes with symbol retained: ", nrow(counts))

# ----------------------------- 4. DESeq2 -------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = clusters["cluster"],
                              design    = ~ cluster)

# Low-count filter: >=10 counts in at least as many samples as the smallest cluster
min_n <- min(table(clusters$cluster))
dds   <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
message("Genes after low-count filter: ", nrow(dds))

# Overall test (LRT): does a gene differ anywhere across the clusters?
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- as.data.frame(results(dds_lrt))
res_lrt$gene <- rownames(res_lrt)
write.csv(res_lrt[order(res_lrt$padj), ],
          file.path(OUT_DIR, "DESeq2_LRT_overall_across_clusters.csv"),
          row.names = FALSE)

# Pairwise Wald contrasts (which clusters differ, and in which direction)
dds_wald <- DESeq(dds)
lv <- levels(clusters$cluster)
pairs <- combn(lv, 2, simplify = FALSE)

pairwise <- lapply(pairs, function(p) {
  # contrast = (second cluster) vs (first cluster); positive LFC = higher in second
  r <- as.data.frame(results(dds_wald, contrast = c("cluster", p[2], p[1])))
  r$gene       <- rownames(r)
  r$comparison <- paste0(p[2], "_vs_", p[1])
  r
})
pairwise <- bind_rows(pairwise)
write.csv(pairwise, file.path(OUT_DIR, "DESeq2_pairwise_Wald_contrasts.csv"),
          row.names = FALSE)

# Number of significant genes per contrast (padj < 0.05, |log2FC| > 1)
sig_summary <- pairwise %>%
  group_by(comparison) %>%
  summarise(n_tested = sum(!is.na(padj)),
            n_sig    = sum(padj < 0.05 & abs(log2FoldChange) > 1, na.rm = TRUE),
            n_up     = sum(padj < 0.05 & log2FoldChange >  1, na.rm = TRUE),
            n_down   = sum(padj < 0.05 & log2FoldChange < -1, na.rm = TRUE))
write.csv(sig_summary, file.path(OUT_DIR, "DESeq2_significant_gene_summary.csv"),
          row.names = FALSE)
print(sig_summary)

# Volcano plot per contrast (PNG + PDF)
# Cutoffs: |log2 fold change| > 1 and adjusted p (FDR) < 0.05
# Labels: top 10 up and top 10 down genes per contrast, ranked by adjusted p
LFC_CUT  <- 1
PADJ_CUT <- 0.05

volc <- pairwise %>%
  filter(!is.na(padj)) %>%
  mutate(padj_plot = pmax(padj, 1e-300),               # avoid infinite -log10(0)
         status = case_when(padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "Up",
                            padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "Down",
                            TRUE ~ "NS"))

volc_labels <- volc %>%
  filter(status != "NS") %>%
  group_by(comparison, status) %>%
  arrange(padj, desc(abs(log2FoldChange)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

p_volc <- ggplot(volc, aes(log2FoldChange, -log10(padj_plot), colour = status)) +
  geom_point(size = 0.7, alpha = 0.6) +
  scale_colour_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey75")) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed") +
  geom_text_repel(data = volc_labels, aes(label = gene), colour = "black",
                  size = 3, max.overlaps = Inf, box.padding = 0.4,
                  segment.size = 0.2, show.legend = FALSE) +
  facet_wrap(~ comparison) +
  theme_bw() +
  labs(x = "log2 fold change", y = "-log10(adjusted p)", colour = NULL,
       subtitle = "Cutoffs: |log2 fold change| > 1; adjusted p (FDR) < 0.05",
       caption = "Positive log2 fold change = higher in the first-named cluster of each contrast; top 10 genes per direction labelled")
ggsave(file.path(OUT_DIR, "plots", "DESeq2_volcano_by_contrast.png"),
       p_volc, width = 15, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, "plots", "DESeq2_volcano_by_contrast.pdf"),
       p_volc, width = 15, height = 6)

# ------------- 4b. HEATMAP: top 100 DESeq2 genes across clusters -------------
# Genes are the 100 most significant in the overall LRT (differ across clusters).
# Values are VST from the same filtered protein-coding DESeq2 object, row z-scored.
cl_cols <- setNames(c("#1B9E77", "#D95F02", "#7570B3")[seq_along(lv)], lv)

top100 <- res_lrt %>%
  filter(!is.na(padj)) %>%
  arrange(padj, pvalue) %>%
  slice_head(n = 100)
write.csv(top100, file.path(OUT_DIR, "DESeq2_top100_LRT_genes.csv"), row.names = FALSE)

vst_deg   <- assay(vst(dds, blind = FALSE))      # protein-coding, filtered genes
ord_c     <- order(clusters$cluster)
z_deg     <- t(scale(t(vst_deg[top100$gene, ord_c])))

ht_deg <- Heatmap(z_deg, name = "z-score",
                  top_annotation = HeatmapAnnotation(
                    Cluster = clusters$cluster[ord_c], col = list(Cluster = cl_cols)),
                  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
                  # sample dendrogram (columns): do the samples regroup by cluster?
                  cluster_columns = TRUE, show_column_dend = TRUE,
                  clustering_method_columns = "ward.D2",
                  column_dend_height = unit(25, "mm"),
                  # gene dendrogram (rows)
                  cluster_rows = TRUE, show_row_dend = TRUE,
                  clustering_method_rows = "ward.D2",
                  row_dend_width = unit(25, "mm"),
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 7),
                  column_title = "Top 100 DESeq2 (LRT) genes across clusters")
png(file.path(OUT_DIR, "plots", "DESeq2_top100_genes_heatmap.png"),
    width = 2800, height = 3600, res = 300)
draw(ht_deg); dev.off()
pdf(file.path(OUT_DIR, "plots", "DESeq2_top100_genes_heatmap.pdf"), width = 9.5, height = 12)
draw(ht_deg); dev.off()

# ----------------------------- 5. ssGSEA -------------------------------------
vst <- read.csv(VST_FILE, row.names = 1, check.names = FALSE)

# Column names in the VST file may be written differently (e.g. "X" prefix,
# dots vs hyphens, "_S25" suffix), so match on a harmonized key
make_key <- function(x) gsub("[.]", "-", sub("_S[0-9]+$", "", sub("^X", "", x)))
idx <- match(make_key(clusters$sample), make_key(colnames(vst)))
if (anyNA(idx))
  stop("Samples not found in VST file: ",
       paste(clusters$sample[is.na(idx)], collapse = ", "),
       "\nVST columns look like: ", paste(head(colnames(vst), 5), collapse = ", "))
vst <- as.matrix(vst[, idx])              # same sample order as metadata
colnames(vst) <- clusters$sample

# Read Hallmark GMT (name <tab> description <tab> genes...)
gmt_lines <- strsplit(readLines(HALLMARK_GMT), "\t")
hallmark  <- setNames(lapply(gmt_lines, function(x) x[-(1:2)]),
                      sapply(gmt_lines, `[`, 1))

# ssGSEA: one enrichment score per sample per pathway
ss_param  <- ssgseaParam(exprData = vst, geneSets = hallmark, normalize = TRUE)
ss_scores <- gsva(ss_param)
rownames(ss_scores) <- sub("^HALLMARK_", "", rownames(ss_scores))
write.csv(ss_scores, file.path(OUT_DIR, "ssGSEA_Hallmark_scores.csv"))

# ------------------- 6. STATISTICS: scores across clusters -------------------
score_long <- as.data.frame(t(ss_scores)) %>%
  tibble::rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "pathway", values_to = "score") %>%
  left_join(clusters[, c("sample", "cluster")], by = "sample")

# Kruskal-Wallis across all clusters (primary family, BH-FDR over 50 pathways)
kw <- score_long %>%
  group_by(pathway) %>%
  summarise(kruskal_p = kruskal.test(score ~ cluster)$p.value,
            .groups = "drop") %>%
  mutate(kruskal_FDR = p.adjust(kruskal_p, method = "BH"))

# Pairwise Wilcoxon tests (follow-up), FDR over all pathway x pair tests
pw <- bind_rows(lapply(pairs, function(p) {
  score_long %>%
    filter(cluster %in% p) %>%
    group_by(pathway) %>%
    summarise(comparison = paste0(p[2], "_vs_", p[1]),
              mean_diff  = mean(score[cluster == p[2]]) - mean(score[cluster == p[1]]),
              wilcox_p   = wilcox.test(score[cluster == p[2]],
                                       score[cluster == p[1]])$p.value,
              .groups = "drop")
})) %>%
  mutate(wilcox_FDR = p.adjust(wilcox_p, method = "BH"))

write.csv(kw, file.path(OUT_DIR, "ssGSEA_Hallmark_KruskalWallis.csv"), row.names = FALSE)
write.csv(pw, file.path(OUT_DIR, "ssGSEA_Hallmark_pairwise_Wilcoxon.csv"), row.names = FALSE)

# ------------------------------ 7. FIGURES -----------------------------------
# 7a. Heatmap: row z-scored ssGSEA scores, columns ordered by cluster
ord_s   <- order(clusters$cluster)
z       <- t(scale(t(ss_scores)))[, ord_s]
top_ann <- HeatmapAnnotation(Cluster = clusters$cluster[ord_s],
                             col = list(Cluster = cl_cols))
ht <- Heatmap(z, name = "z-score", top_annotation = top_ann,
              col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
              cluster_columns = FALSE, show_column_names = TRUE,
              column_split = clusters$cluster[ord_s],
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 7))
png(file.path(OUT_DIR, "plots", "ssGSEA_Hallmark_heatmap.png"),
    width = 3000, height = 2600, res = 300)
draw(ht); dev.off()
pdf(file.path(OUT_DIR, "plots", "ssGSEA_Hallmark_heatmap.pdf"), width = 10, height = 8.5)
draw(ht); dev.off()

# 7b. One boxplot per pathway significant across clusters (Kruskal FDR < 0.05),
#     titled with the pathway name and its FDR-adjusted p-value
fmt_p <- function(x) {
  if (x < 1e-4) format(x, digits = 3, scientific = TRUE)
  else format(signif(x, 3), scientific = FALSE)
}

sig_paths <- kw$pathway[kw$kruskal_FDR < 0.05]
if (length(sig_paths) > 0) {
  box_dir <- file.path(OUT_DIR, "plots", "ssGSEA_boxplots")
  dir.create(box_dir, recursive = TRUE, showWarnings = FALSE)

  for (pth in sig_paths) {
    fdr <- kw$kruskal_FDR[kw$pathway == pth]

    p_box <- score_long %>%
      filter(pathway == pth) %>%
      ggplot(aes(cluster, score, fill = cluster)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4, colour = "black", linewidth = 0.8) +
      geom_jitter(width = 0.15, size = 2, colour = "black", alpha = 0.8) +
      scale_fill_manual(values = cl_cols) +
      scale_x_discrete(labels = sub("^C", "", lv)) +
      labs(title = paste0("HALLMARK_", pth, "\nFDR-adjusted p = ", fmt_p(fdr)),
           x = "Cluster", y = "ssGSEA Enrichment Score") +
      theme_classic(base_size = 16) +
      theme(legend.position = "none",
            plot.title  = element_text(face = "bold", hjust = 0.5, size = 16),
            axis.title  = element_text(face = "bold"))

    ggsave(file.path(box_dir, paste0("ssGSEA_", pth, ".png")),
           p_box, width = 7.5, height = 5.5, dpi = 300)
    ggsave(file.path(box_dir, paste0("ssGSEA_", pth, ".pdf")),
           p_box, width = 7.5, height = 5.5)
  }
} else {
  message("No Hallmark pathway reached Kruskal-Wallis FDR < 0.05.")
}

message("Done. Outputs in: ", OUT_DIR)
