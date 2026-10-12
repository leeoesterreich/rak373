###############################################################################
# Which pILC tumor cluster does each PDO resemble?
#   PDOs: BO-046, BO-114, BO-189, BO-268
#   Tumors: pILC samples with NMF clusters C1, C2, C3 (k = 3)
#
# IDEA
#   ssGSEA scores are only comparable between samples that go through the SAME
#   steps, so tumors and PDOs are transformed (VST) and scored (ssGSEA) in ONE
#   run. The clusters and all statistics come from the TUMORS ONLY. Each PDO is
#   then compared with the mean pathway profile of each cluster.
#
# WORKFLOW
#   1. Tumor raw counts (4_Count_file_pILC) + PDO raw counts (PDO_raw_counts.csv)
#   2. Merge on shared genes (Ensembl IDs without version)
#   3. Keep protein-coding genes with a symbol, one row per symbol
#   4. Count filter, VST on all samples together (blind = TRUE)
#   5. Optional variance filter (VAR_FILTER_PCT; 0 = off)
#   6. ssGSEA Hallmark (normalize = TRUE) on all samples together
#   7. Tumors only: Kruskal-Wallis across clusters (BH-FDR) -> significant pathways
#   8. Cluster similarity of each PDO:
#        - pick the pathways that separate the clusters (tumors only)
#        - z-score each pathway using the TUMOR mean and SD
#        - cluster profile = mean z-score of the tumors in that cluster
#        - Spearman correlation (and Euclidean distance) between a PDO and
#          each cluster profile; the highest correlation = closest cluster
#   9. Two versions: all pathways, and without immune/inflammatory pathways
#      (organoids contain no immune cells, so these pathways are low in PDOs
#      and could bias the result)
#  10. Leave-one-out test on the TUMORS: how often is a tumor assigned to its
#      own cluster? This tells you how much to trust the PDO assignment.
#  11. Figures: boxplots (tumors by cluster, PDOs as separate points and lines),
#      heatmap (tumors by cluster + PDO block), PCA of ssGSEA scores
#
# NOTES
#   - PDOs are NOT part of any cluster test and NEVER change the clusters or
#     the cluster profiles. No p-values are given for the PDOs (n = 1 each).
#   - The result is "closest to", not proof of origin.
#   - Tumor VST is recomputed together with the PDOs, so tumor scores differ
#     slightly from the earlier tumor-only run. Step 9b repeats the pipeline on
#     the tumors alone and saves how much the scores and the significant
#     pathways change (Joint_vs_tumor_only_stability.csv).
#   - By default the pathways used to compare PDOs with clusters are the fixed
#     list from the earlier tumor-only analysis (FIXED_PATHWAYS), so the
#     choice does not depend on the PDOs. The leave-one-out accuracy is then
#     slightly optimistic, because that list was chosen using all tumors.
#   - No dates in output file names. Seed 42.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

SEED <- 42
set.seed(SEED)

# ------------------------------ 1. SETTINGS ----------------------------------
ROOT         <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
COUNTS_DIR   <- file.path(ROOT, "4_Count_file_pILC")
CLUSTER_FILE <- file.path(ROOT, "NMF", "GenePattern_NMFConsensus",
                          "vst_transformed_data.consensus.k.3.gct")
HALLMARK_GMT <- file.path(ROOT, "GSEA", "Dataset", "h.all.v2025.1.Hs.symbols.gmt")
PDO_COUNTS   <- file.path(ROOT, "PDO", "PDO_raw_counts.csv")

OUT_DIR <- file.path(ROOT, "PDO", "Joint_with_tumors")
dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

MIN_COUNT      <- 10    # count filter: >= MIN_COUNT counts ...
VAR_FILTER_PCT <- 0.10  # remove the bottom 10% of genes by VST variance (0 = off)
FDR_CUT        <- 0.05  # pathways separating the clusters (Kruskal-Wallis, tumors)
MIN_PATHS      <- 5     # if fewer pathways pass, use the top MIN_PATHS by p-value

# Pathways driven by immune cells / microenvironment (absent in organoids)
IMMUNE_PATHWAYS <- c("ALLOGRAFT_REJECTION", "INTERFERON_ALPHA_RESPONSE",
                     "INTERFERON_GAMMA_RESPONSE", "IL6_JAK_STAT3_SIGNALING",
                     "TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE",
                     "COMPLEMENT", "IL2_STAT5_SIGNALING")

# Pathways that separated the tumor clusters in the EARLIER tumor-only analysis
# (Kruskal-Wallis FDR < 0.05). Using this fixed list keeps the pathway choice
# independent of the PDOs and of the joint VST. Set to NULL to recompute the
# list from the joint run instead.
FIXED_PATHWAYS <- c("ALLOGRAFT_REJECTION", "ANDROGEN_RESPONSE",
                    "IL6_JAK_STAT3_SIGNALING", "INTERFERON_ALPHA_RESPONSE",
                    "INTERFERON_GAMMA_RESPONSE", "KRAS_SIGNALING_UP",
                    "NOTCH_SIGNALING", "REACTIVE_OXYGEN_SPECIES_PATHWAY",
                    "SPERMATOGENESIS", "TNFA_SIGNALING_VIA_NFKB",
                    "WNT_BETA_CATENIN_SIGNALING")

# --------------------------- 2. TUMOR CLUSTERS -------------------------------
clusters <- read.delim(CLUSTER_FILE, skip = 2, check.names = FALSE,
                       stringsAsFactors = FALSE)
clusters <- data.frame(sample  = clusters$Name,
                       cluster = clusters$membership.ordered,
                       stringsAsFactors = FALSE)
clusters$cluster <- factor(paste0("C", clusters$cluster),
                           levels = sort(unique(paste0("C", clusters$cluster))))
rownames(clusters) <- clusters$sample
print(table(clusters$cluster))
lv <- levels(clusters$cluster)

# ------------------------- 3. TUMOR RAW COUNTS -------------------------------
count_files <- list.files(COUNTS_DIR, pattern = "\\.txt$", full.names = TRUE)
file_key    <- sub("\\.txt$", "", basename(count_files))
cluster_key <- gsub("\\.", "-", sub("_S[0-9]+$", "", clusters$sample))
missing <- setdiff(cluster_key, file_key)
if (length(missing) > 0) stop("No count file found for: ", paste(missing, collapse = ", "))
count_files <- count_files[match(cluster_key, file_key)]

read_one <- function(f) {
  x <- read.delim(f, header = FALSE, col.names = c("gene", "count"),
                  stringsAsFactors = FALSE)
  x[grepl("^ENSG", x$gene), ]
}
count_list <- lapply(count_files, read_one)
stopifnot(all(sapply(count_list, function(x) identical(x$gene, count_list[[1]]$gene))))
tum <- sapply(count_list, function(x) x$count)
rownames(tum) <- count_list[[1]]$gene
colnames(tum) <- clusters$sample
mode(tum) <- "integer"

strip_version <- function(x) sub("\\..*$", "", x)
rownames(tum) <- strip_version(rownames(tum))
tum <- tum[!duplicated(rownames(tum)), ]

# --------------------------- 4. PDO RAW COUNTS -------------------------------
pdo_raw <- read.csv(PDO_COUNTS, check.names = FALSE, stringsAsFactors = FALSE)
pdo <- as.matrix(pdo_raw[, setdiff(colnames(pdo_raw), "gene")])
rownames(pdo) <- strip_version(pdo_raw$gene)
pdo <- pdo[!duplicated(rownames(pdo)), ]
mode(pdo) <- "integer"
pdo_names <- colnames(pdo)

# ------------------------------ 5. MERGE -------------------------------------
shared <- intersect(rownames(tum), rownames(pdo))
message("Genes: tumors ", nrow(tum), ", PDOs ", nrow(pdo), ", shared ", length(shared))
counts <- cbind(tum[shared, ], pdo[shared, ])
group  <- c(as.character(clusters$cluster), rep("PDO", length(pdo_names)))
names(group) <- colnames(counts)
tumor_samples <- clusters$sample

# ------------------------ 6. GENE CLEANING (SYMBOLS) -------------------------
anno <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(rownames(counts)),
                              keytype = "ENSEMBL", columns = c("SYMBOL", "GENETYPE"))
anno   <- anno[!duplicated(anno$ENSEMBL), ]
symbol <- anno$SYMBOL[match(rownames(counts), anno$ENSEMBL)]
gtype  <- anno$GENETYPE[match(rownames(counts), anno$ENSEMBL)]
keep   <- !is.na(symbol) & symbol != "" & gtype %in% "protein-coding"
counts <- counts[keep, ]
symbol <- symbol[keep]
ord    <- order(rowMeans(counts), decreasing = TRUE)
counts <- counts[ord, ]
symbol <- symbol[ord]
counts <- counts[!duplicated(symbol), ]
rownames(counts) <- symbol[!duplicated(symbol)]
n_coding <- nrow(counts)

# ----------------------- 7. COUNT FILTER + VST (ALL TOGETHER) ----------------
min_n  <- min(table(clusters$cluster))   # same rule as the DESeq2 script
counts <- counts[rowSums(counts >= MIN_COUNT) >= min_n, ]
n_count <- nrow(counts)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = data.frame(group = group, row.names = colnames(counts)),
                              design    = ~ 1)
vst_mat <- assay(vst(dds, blind = TRUE))

if (VAR_FILTER_PCT > 0) {
  gv      <- apply(vst_mat, 1, var)
  vst_mat <- vst_mat[gv > quantile(gv, probs = VAR_FILTER_PCT), ]
}
write.csv(data.frame(step  = c("Shared genes", "Protein-coding with symbol",
                               "After count filter", "After variance filter"),
                     genes = c(length(shared), n_coding, n_count, nrow(vst_mat))),
          file.path(OUT_DIR, "Joint_gene_filtering_summary.csv"), row.names = FALSE)

# ------------------------------ 8. ssGSEA ------------------------------------
gmt_lines <- strsplit(readLines(HALLMARK_GMT), "\t")
hallmark  <- setNames(lapply(gmt_lines, function(x) x[-(1:2)]),
                      sapply(gmt_lines, `[`, 1))
ss <- gsva(ssgseaParam(exprData = vst_mat, geneSets = hallmark, normalize = TRUE))
rownames(ss) <- sub("^HALLMARK_", "", rownames(ss))
write.csv(ss, file.path(OUT_DIR, "Joint_ssGSEA_Hallmark_scores_tumors_and_PDOs.csv"))

ss_tum <- ss[, tumor_samples]
ss_pdo <- ss[, pdo_names, drop = FALSE]

# ------------- 9. TUMORS ONLY: KRUSKAL-WALLIS ACROSS CLUSTERS ----------------
kw_p <- apply(ss_tum, 1, function(x) kruskal.test(x ~ clusters$cluster)$p.value)
kw <- data.frame(pathway = names(kw_p), kruskal_p = kw_p,
                 kruskal_FDR = p.adjust(kw_p, method = "BH"), row.names = NULL) %>%
  arrange(kruskal_p)
write.csv(kw, file.path(OUT_DIR, "Joint_tumors_KruskalWallis.csv"), row.names = FALSE)
sig_paths <- kw$pathway[kw$kruskal_FDR < FDR_CUT]
message("Significant pathways across tumor clusters (FDR < ", FDR_CUT, "): ",
        length(sig_paths), "\n  ", paste(sig_paths, collapse = ", "),
        "\n  (compare with your earlier tumor-only list)")

# ----------- 9b. STABILITY: does adding the PDOs change the tumor result? ----
# Same pipeline on the TUMORS ALONE, then compare with the joint run.
joint_sig <- sig_paths
cnt_t <- counts[, tumor_samples]
cnt_t <- cnt_t[rowSums(cnt_t >= MIN_COUNT) >= min_n, ]
dds_t <- DESeqDataSetFromMatrix(countData = cnt_t,
                                colData   = data.frame(id = colnames(cnt_t),
                                                       row.names = colnames(cnt_t)),
                                design    = ~ 1)
vst_t <- assay(vst(dds_t, blind = TRUE))
if (VAR_FILTER_PCT > 0) {
  gv_t  <- apply(vst_t, 1, var)
  vst_t <- vst_t[gv_t > quantile(gv_t, probs = VAR_FILTER_PCT), ]
}
ss_t <- gsva(ssgseaParam(exprData = vst_t, geneSets = hallmark, normalize = TRUE))
rownames(ss_t) <- sub("^HALLMARK_", "", rownames(ss_t))
ss_t <- ss_t[rownames(ss_tum), colnames(ss_tum)]

kw_t_p <- apply(ss_t, 1, function(x) kruskal.test(x ~ clusters$cluster)$p.value)
stab <- data.frame(
  pathway = rownames(ss_tum),
  spearman_tumor_scores_joint_vs_tumor_only =
    unname(sapply(rownames(ss_tum), function(pw) cor(ss_t[pw, ], ss_tum[pw,], method = "spearman"))),
  KW_FDR_tumor_only = p.adjust(kw_t_p, method = "BH"),
  KW_FDR_joint      = kw$kruskal_FDR[match(rownames(ss_tum), kw$pathway)],
  row.names = NULL)
stab$sig_tumor_only <- stab$KW_FDR_tumor_only < FDR_CUT
stab$sig_joint      <- stab$KW_FDR_joint < FDR_CUT
write.csv(stab, file.path(OUT_DIR, "Joint_vs_tumor_only_stability.csv"), row.names = FALSE)
message("Stability: median Spearman of tumor scores, joint vs tumor-only = ",
        round(median(stab$spearman_tumor_scores_joint_vs_tumor_only), 3),
        "; significant in tumor-only: ", sum(stab$sig_tumor_only),
        ", in joint: ", sum(stab$sig_joint),
        ", in both: ", sum(stab$sig_tumor_only & stab$sig_joint))

# Pathways shown in the boxplots: the fixed list if used, otherwise the joint result
if (!is.null(FIXED_PATHWAYS)) sig_paths <- intersect(FIXED_PATHWAYS, rownames(ss_tum))

# ------------------- 10. CLUSTER SIMILARITY FUNCTION -------------------------
# train = tumor ssGSEA (pathways x samples), train_cl = tumor clusters (factor),
# test = samples to classify (pathways x samples).
# Everything (pathway choice, z-scoring, cluster profiles) uses TRAIN only.
assign_cluster <- function(train, train_cl, test, exclude = character(0)) {
  if (!is.null(FIXED_PATHWAYS)) {
    # fixed list from the earlier tumor-only analysis (minus excluded pathways)
    sel <- setdiff(intersect(FIXED_PATHWAYS, rownames(train)), exclude)
  } else {
    paths <- setdiff(rownames(train), exclude)
    p     <- sapply(paths, function(pw) kruskal.test(train[pw, ] ~ train_cl)$p.value)
    fdr   <- p.adjust(p, method = "BH")
    sel   <- names(fdr)[fdr < FDR_CUT]
    if (length(sel) < MIN_PATHS) sel <- names(sort(p))[seq_len(MIN_PATHS)]
  }

  mu  <- rowMeans(train[sel, , drop = FALSE])
  sdv <- apply(train[sel, , drop = FALSE], 1, sd)
  ztr <- (train[sel, , drop = FALSE] - mu) / sdv
  zte <- (test[sel, , drop = FALSE]  - mu) / sdv

  cent <- sapply(levels(train_cl), function(k) rowMeans(ztr[, train_cl == k, drop = FALSE]))
  cor_m  <- apply(zte, 2, function(x) apply(cent, 2, function(cc) cor(x, cc, method = "spearman")))
  dist_m <- apply(zte, 2, function(x) apply(cent, 2, function(cc) sqrt(sum((x - cc)^2))))
  cor_m  <- matrix(cor_m,  nrow = ncol(cent), dimnames = list(colnames(cent), colnames(test)))
  dist_m <- matrix(dist_m, nrow = ncol(cent), dimnames = list(colnames(cent), colnames(test)))
  list(selected = sel, cor = cor_m, dist = dist_m)
}

# ----------------- 11. ASSIGN THE PDOs (TWO VERSIONS) ------------------------
run_version <- function(label, exclude) {
  r <- assign_cluster(ss_tum, clusters$cluster, ss_pdo, exclude)
  out <- data.frame(version = label, PDO = colnames(ss_pdo), row.names = NULL)
  for (k in lv) out[[paste0("spearman_", k)]]  <- r$cor[k, ]
  for (k in lv) out[[paste0("distance_", k)]]  <- r$dist[k, ]
  out$closest_by_correlation <- rownames(r$cor)[apply(r$cor, 2, which.max)]
  out$closest_by_distance    <- rownames(r$dist)[apply(r$dist, 2, which.min)]
  out$n_pathways_used        <- length(r$selected)
  list(table = out, pathways = r$selected)
}
v_all <- run_version("all_pathways", character(0))
v_noi <- run_version("without_immune_pathways", IMMUNE_PATHWAYS)

similarity <- bind_rows(v_all$table, v_noi$table)
write.csv(similarity, file.path(OUT_DIR, "PDO_cluster_similarity.csv"), row.names = FALSE)
write.csv(data.frame(version = c(rep("all_pathways", length(v_all$pathways)),
                                 rep("without_immune_pathways", length(v_noi$pathways))),
                     pathway = c(v_all$pathways, v_noi$pathways)),
          file.path(OUT_DIR, "PDO_cluster_similarity_pathways_used.csv"), row.names = FALSE)
print(similarity)

# ------------- 12. LEAVE-ONE-OUT TEST ON THE TUMORS --------------------------
# Each tumor is removed, everything is rebuilt without it, and it is assigned.
loo_run <- function(exclude) {
  pred <- sapply(seq_len(ncol(ss_tum)), function(i) {
    r <- assign_cluster(ss_tum[, -i], clusters$cluster[-i],
                        ss_tum[, i, drop = FALSE], exclude)
    rownames(r$cor)[which.max(r$cor[, 1])]
  })
  factor(pred, levels = lv)
}
loo_all <- loo_run(character(0))
loo_noi <- loo_run(IMMUNE_PATHWAYS)

loo_summary <- data.frame(
  version  = c("all_pathways", "without_immune_pathways"),
  accuracy = c(mean(loo_all == clusters$cluster), mean(loo_noi == clusters$cluster)))
for (k in lv) {
  loo_summary[[paste0("accuracy_", k)]] <- c(
    mean(loo_all[clusters$cluster == k] == k),
    mean(loo_noi[clusters$cluster == k] == k))
}
write.csv(loo_summary, file.path(OUT_DIR, "Tumor_leave_one_out_accuracy.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(table(true = clusters$cluster, predicted = loo_all)),
          file.path(OUT_DIR, "Tumor_leave_one_out_confusion_all_pathways.csv"))
write.csv(as.data.frame.matrix(table(true = clusters$cluster, predicted = loo_noi)),
          file.path(OUT_DIR, "Tumor_leave_one_out_confusion_without_immune.csv"))
print(loo_summary)

# ------------------------------ 13. FIGURES ----------------------------------
cl_cols  <- setNames(c("#1B9E77", "#D95F02", "#7570B3")[seq_along(lv)], lv)
pdo_cols <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")[seq_along(pdo_names)],
                     pdo_names)

# 13a. Boxplots: tumors by cluster; each PDO as a point at "PDO" and a dashed
#      line across the panel (to see which box its level matches)
if (length(sig_paths) > 0) {
  tum_long <- as.data.frame(t(ss_tum[sig_paths, , drop = FALSE])) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "pathway", values_to = "score") %>%
    left_join(clusters[, c("sample", "cluster")], by = "sample")
  pdo_long <- as.data.frame(t(ss_pdo[sig_paths, , drop = FALSE])) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "pathway", values_to = "score")

  p_box <- ggplot() +
    geom_boxplot(data = tum_long, aes(cluster, score, fill = cluster),
                 outlier.shape = NA, alpha = 0.6, colour = "black") +
    geom_jitter(data = tum_long, aes(cluster, score), width = 0.15, size = 1, alpha = 0.7) +
    geom_hline(data = pdo_long, aes(yintercept = score, colour = sample),
               linetype = "dashed", linewidth = 0.5) +
    geom_point(data = pdo_long, aes(x = "PDO", y = score, colour = sample),
               shape = 18, size = 3.5) +
    scale_fill_manual(values = cl_cols) +
    scale_colour_manual(values = pdo_cols, name = "PDO") +
    scale_x_discrete(limits = c(lv, "PDO")) +
    facet_wrap(~ pathway, scales = "free_y") +
    labs(x = NULL, y = "ssGSEA score",
         subtitle = "Boxes: tumors by cluster. Diamonds and dashed lines: PDOs (scored in the same run)") +
    theme_bw()
  ggsave(file.path(OUT_DIR, "plots", "Joint_boxplots_tumor_clusters_with_PDO.png"),
         p_box, width = 15, height = 11, dpi = 300)
  ggsave(file.path(OUT_DIR, "plots", "Joint_boxplots_tumor_clusters_with_PDO.pdf"),
         p_box, width = 15, height = 11)
}

# 13b. Heatmap: row z-score from the TUMOR mean and SD; tumors by cluster, PDO block
mu  <- rowMeans(ss_tum)
sdv <- apply(ss_tum, 1, sd)
z   <- (ss - mu) / sdv
col_order <- c(tumor_samples[order(clusters$cluster)], pdo_names)
z         <- z[, col_order]
grp       <- factor(group[col_order], levels = c(lv, "PDO"))

ann_cols <- c(cl_cols, PDO = "grey40")
top_ann  <- HeatmapAnnotation(Group = grp, col = list(Group = ann_cols))
ht <- Heatmap(z, name = "z-score\n(tumor-based)", top_annotation = top_ann,
              col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
              cluster_columns = FALSE, column_split = grp,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 6))
png(file.path(OUT_DIR, "plots", "Joint_ssGSEA_heatmap_tumors_and_PDO.png"),
    width = 3400, height = 2600, res = 300)
draw(ht); dev.off()
pdf(file.path(OUT_DIR, "plots", "Joint_ssGSEA_heatmap_tumors_and_PDO.pdf"), width = 11.5, height = 8.5)
draw(ht); dev.off()

# 13c. PCA of ssGSEA scores: do the PDOs sit apart from all tumors (technical/
#      biological separation) or inside the tumor cloud?
pc  <- prcomp(t(ss), scale. = TRUE)
pcd <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                  group = factor(group[rownames(pc$x)], levels = c(lv, "PDO")))
ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
p_pca <- ggplot(pcd, aes(PC1, PC2, colour = group, shape = group)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(data = subset(pcd, group == "PDO"), aes(label = sample),
                  colour = "black", size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = ann_cols) +
  scale_shape_manual(values = c(rep(16, length(lv)), 18)) +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       colour = NULL, shape = NULL) +
  theme_bw(base_size = 13)
ggsave(file.path(OUT_DIR, "plots", "Joint_ssGSEA_PCA_tumors_and_PDO.png"),
       p_pca, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(OUT_DIR, "plots", "Joint_ssGSEA_PCA_tumors_and_PDO.pdf"),
       p_pca, width = 7.5, height = 6)

message("Done. Outputs in: ", OUT_DIR)
sessionInfo()
