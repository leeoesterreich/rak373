# pILC RNA-seq: PAM50 intrinsic subtyping on the stranded (corrected) VST matrix.
# Method: genefu PAM50 (Parker 2009 centroids, Spearman correlation of each tumor to the
# 5 centroids; the tumor takes the subtype of the highest correlation).
# Outputs (no dates in names): subtype calls, subtype probabilities, subtype-vs-cluster table + tests,
# clinical features (tumor cellularity, TILs, pleomorphic component) by subtype, figures (PNG + PDF).
# Cluster labels: 50-run consensus NMF (R_NMFConsensus_50runs/Cluster_membership_k3.csv).

suppressPackageStartupMessages({
  library(genefu)          # BiocManager::install("genefu")
  library(AnnotationDbi)
  library(org.Hs.eg.db)    # symbol -> Entrez ID (PAM50 centroids are matched by Entrez ID)
  library(ggplot2)
})
SEED <- 42
set.seed(SEED)

# ---------------------------------------------------------------- 1. Paths
root    <- Sys.getenv("PILC_ROOT", "/ix1/alee/LO_LAB/Personal/Rahul/pILC")
vst_dir <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47")
vst_f   <- file.path(vst_dir, "htseq_VST_matrix_all_expression_filtered_genes.rds")
anno_f  <- file.path(vst_dir, "gene_annotation_from_GTF.tsv")
# 50-run cluster labels (same file as the PCA and correlation-heatmap figures)
clus_f  <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
clin_f  <- file.path(root, "Clinical_correlation/pILC_ClinicoPathological_details.csv")
if (!file.exists(clus_f)) stop("Cluster file not found: ", clus_f)
out_dir <- file.path(root, "RNASeq/PAM50/Stranded")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

norm_id <- function(x) {                        # TP09.9_S34 / TP09-9 -> TP09-09
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}

# ---------------------------------------------------------------- 2. Load VST matrix (genes x samples)
vst <- as.matrix(readRDS(vst_f))
if (nrow(vst) < ncol(vst)) vst <- t(vst)
message("VST matrix: ", nrow(vst), " genes x ", ncol(vst), " samples")
message("First row names: ", paste(head(rownames(vst), 3), collapse = ", "))

# ---------------------------------------------------------------- 3. Gene symbols (one row per symbol)
ids <- rownames(vst)
if (mean(grepl("^ENSG", ids)) > 0.5) {          # Ensembl IDs -> symbols via the GTF annotation
  anno <- read.delim(anno_f, check.names = FALSE, stringsAsFactors = FALSE)
  id_col  <- grep("gene_id|^Name$|ensembl", colnames(anno), ignore.case = TRUE, value = TRUE)[1]
  sym_col <- grep("gene_name|symbol", colnames(anno), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col) || is.na(sym_col))
    stop("Cannot find ID/symbol columns in the annotation. Columns: ", paste(colnames(anno), collapse = ", "))
  sym <- anno[[sym_col]][match(sub("\\..*$", "", ids), sub("\\..*$", "", anno[[id_col]]))]
} else {
  sym <- ids
}
keep <- !is.na(sym) & sym != ""
vst <- vst[keep, ]; sym <- sym[keep]
ord <- order(rowMeans(vst), decreasing = TRUE)   # for duplicated symbols keep the highest-expressed row
vst <- vst[ord, ]; sym <- sym[ord]
dup <- duplicated(sym)
vst <- vst[!dup, ]; rownames(vst) <- sym[!dup]

# ---------------------------------------------------------------- 4. Keep the PAM50 genes
data(pam50.robust, package = "genefu")
cmap <- pam50.robust$centroids.map                 # PAM50 genes with Entrez IDs
ent  <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(vst), column = "ENTREZID",
                              keytype = "SYMBOL", multiVals = "first")
use  <- !is.na(ent) & ent %in% as.character(cmap$EntrezGene.ID)
expr50 <- vst[use, , drop = FALSE]
annot  <- data.frame(probe = rownames(expr50), Gene.Symbol = rownames(expr50),
                     EntrezGene.ID = as.numeric(ent[use]), row.names = rownames(expr50))
genes_used <- data.frame(PAM50_probe = cmap$probe, EntrezGene.ID = cmap$EntrezGene.ID,
                         in_data = as.character(cmap$EntrezGene.ID) %in% as.character(annot$EntrezGene.ID))
write.csv(genes_used, file.path(out_dir, "PAM50_genes_used.csv"), row.names = FALSE)
message("PAM50 genes found in the data: ", sum(genes_used$in_data), " / ", nrow(genes_used))
if (any(!genes_used$in_data))
  message("Missing: ", paste(genes_used$PAM50_probe[!genes_used$in_data], collapse = ", "))

# ---------------------------------------------------------------- 5. PAM50 prediction
# data = samples x genes. genefu centers each gene across the cohort before correlating with
# the centroids, so calls depend on the cohort composition (see the note on ER balance).
pred <- molecular.subtyping(sbt.model = "pam50", data = t(expr50), annot = annot, do.mapping = TRUE)
sub_names <- c(Basal = "Basal-like", Her2 = "HER2-enriched", LumA = "Luminal A",
               LumB = "Luminal B", Normal = "Normal-like")
subtype <- unname(sub_names[as.character(pred$subtype)])   # readable subtype names
names(subtype) <- names(pred$subtype)                       # keep sample names

calls <- data.frame(Sample = names(subtype), Key = norm_id(names(subtype)),
                    PAM50 = factor(subtype, levels = unname(sub_names)))
if (is.matrix(pred$cor)) {                          # correlation of each tumor to each centroid
  cm <- pred$cor[calls$Sample, , drop = FALSE]
  colnames(cm) <- paste0("cor_", colnames(cm))
  calls <- cbind(calls, cm)
}

# ---------------------------------------------------------------- 5b. Subtype probabilities (genefu) and centroid correlations
# subtype.proba is taken exactly as genefu returns it. The script prints whether the rows sum to 1,
# so you can see whether they behave as probabilities. Treat them as the relative confidence of the call.
prob <- NULL
pr <- pred$subtype.proba
if (!is.null(pr) && (is.matrix(pr) || is.data.frame(pr))) {
  prob <- as.matrix(pr)
  if (is.null(rownames(prob)) && nrow(prob) == length(subtype)) rownames(prob) <- names(subtype)
  colnames(prob) <- ifelse(colnames(prob) %in% names(sub_names), sub_names[colnames(prob)], colnames(prob))
  message("subtype.proba found: row sums range ", paste(round(range(rowSums(prob)), 3), collapse = " - "))
} else {
  message("genefu did not return subtype.proba; the probability plot is skipped (centroid correlations are still plotted)")
}

# ---------------------------------------------------------------- 6. Add RNA cluster and clinical data
# "<5" style values are censored: use half the limit (same rule as the cluster-vs-clinical script)
pct <- function(x) {
  x <- trimws(as.character(x))
  out <- suppressWarnings(as.numeric(sub("^<", "", x)))
  out[grepl("^<", x) & !is.na(out)] <- out[grepl("^<", x) & !is.na(out)] / 2
  out
}
clus <- read.csv(clus_f); clus$Key <- norm_id(clus$sample)
clin <- read.csv(clin_f, check.names = FALSE, strip.white = TRUE); clin$Key <- norm_id(clin$Sample)
m <- match(calls$Key, clin$Key)
calls$Cluster <- factor(clus$cluster[match(calls$Key, clus$Key)], levels = 1:3)
calls$Tumor_cellularity <- pct(clin[["tumor_cellularity (%)"]][m])
calls$TILs              <- pct(clin[["TILs (%)"]][m])
calls$Pleomorphic_component <- pct(clin[["pleomorphic_invasive_compoent_cellularity (%)"]][m])
calls$ER  <- trimws(clin[["Tumor Marker #1 (ER)"]][m])
calls$HER2 <- trimws(clin[["BREAST: Her-2/Neu-Desc"]][m])
if (!is.null(prob)) {
  pm <- prob[calls$Sample, , drop = FALSE]
  colnames(pm) <- paste0("prob_", colnames(pm))
  calls <- cbind(calls, pm)
  write.csv(data.frame(Sample = calls$Sample, Key = calls$Key, Cluster = calls$Cluster, PAM50 = calls$PAM50, pm),
            file.path(out_dir, "PAM50_subtype_probability.csv"), row.names = FALSE)
}
write.csv(calls, file.path(out_dir, "PAM50_subtype_calls.csv"), row.names = FALSE)
print(table(calls$PAM50, useNA = "ifany"))

# ---------------------------------------------------------------- 7. Subtype vs RNA cluster, and vs clinical features
tab <- table(PAM50 = calls$PAM50, Cluster = calls$Cluster)
write.csv(as.data.frame.matrix(tab), file.path(out_dir, "PAM50_subtype_by_RNA_cluster.csv"))
tab_nz <- tab[rowSums(tab) > 0, , drop = FALSE]
p_cluster <- if (nrow(tab_nz) > 1) fisher.test(tab_nz, simulate.p.value = TRUE, B = 1e5)$p.value else NA_real_

feat <- c(Tumor_cellularity = "Tumor cellularity (%)", TILs = "TILs (%)",
          Pleomorphic_component = "Pleomorphic component cellularity (%)")
kw <- sapply(names(feat), function(v) {
  ok <- !is.na(calls[[v]])
  if (length(unique(calls$PAM50[ok])) > 1) kruskal.test(calls[[v]][ok] ~ droplevels(calls$PAM50[ok]))$p.value else NA_real_
})
n_ok <- sapply(names(feat), function(v) sum(!is.na(calls[[v]])))
stats <- data.frame(Test = c("PAM50 subtype vs RNA cluster (Fisher, simulated p, seed 42)",
                             paste(feat, "across PAM50 subtypes (Kruskal-Wallis)")),
                    p_value = c(p_cluster, kw), N = c(sum(tab), n_ok))
write.csv(stats, file.path(out_dir, "PAM50_stats.csv"), row.names = FALSE)
print(stats)

# Table paired with the clinical figure: n, median [Q1-Q3] per subtype and feature
desc <- do.call(rbind, lapply(names(feat), function(v) do.call(rbind, lapply(levels(calls$PAM50), function(s) {
  x <- calls[[v]][calls$PAM50 == s & !is.na(calls[[v]])]
  data.frame(Feature = feat[[v]], PAM50 = s, n = length(x),
             median = if (length(x)) median(x) else NA, Q1 = if (length(x)) quantile(x, .25) else NA,
             Q3 = if (length(x)) quantile(x, .75) else NA, row.names = NULL)
}))))
write.csv(desc, file.path(out_dir, "PAM50_clinical_features_by_subtype.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 8. Figures
sub_cols <- c("Luminal A" = "#1F4E9C", "Luminal B" = "#6BAED6", "HER2-enriched" = "#E377C2",
              "Basal-like" = "#D62728", "Normal-like" = "#2CA02C")
th <- theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))
save_both <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
}

# 8a. Subtype by RNA cluster
p1 <- ggplot(calls[!is.na(calls$Cluster), ], aes(Cluster, fill = PAM50)) +
  geom_bar(position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sub_cols, drop = FALSE) +
  labs(title = "PAM50 subtype by RNA cluster", y = "Number of tumors", fill = "PAM50") + th
save_both(p1, "PAM50_subtype_by_RNA_cluster", 6.5, 5)

# 8b. Tumor cellularity, TILs and pleomorphic component by subtype (one panel each, Kruskal-Wallis p in the title)
long <- do.call(rbind, lapply(names(feat), function(v)
  data.frame(PAM50 = calls$PAM50, Value = calls[[v]],
             Feature = sprintf("%s\nKruskal-Wallis p = %s", feat[[v]], ifelse(is.na(kw[[v]]), "NA", sprintf("%.3f", kw[[v]]))))))
long <- long[!is.na(long$Value), ]
long$Feature <- factor(long$Feature, levels = unique(long$Feature))
p2 <- ggplot(long, aes(PAM50, Value, fill = PAM50)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + geom_jitter(width = 0.1, size = 2) +
  facet_wrap(~Feature, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = sub_cols) +
  labs(title = "Clinical features by PAM50 subtype", x = NULL, y = "Percent (%)") +
  th + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1),
             strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11))
save_both(p2, "PAM50_clinical_features_by_subtype", 13, 5)

# 8c. Subtype probability per tumor (samples grouped by RNA cluster, then by called subtype)
ord <- calls[order(calls$Cluster, calls$PAM50, calls$Sample), ]
if (!is.null(prob)) {
  pl <- do.call(rbind, lapply(colnames(prob), function(s)
    data.frame(Key = calls$Key, Cluster = calls$Cluster, Probability = prob[calls$Sample, s], Subtype = s)))
  pl$Subtype <- factor(pl$Subtype, levels = names(sub_cols))
  pl$Key <- factor(pl$Key, levels = ord$Key)
  pl <- pl[!is.na(pl$Cluster), ]
  p3 <- ggplot(pl, aes(Key, Probability, fill = Subtype)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.2) +
    facet_grid(~paste("Cluster", Cluster), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = sub_cols, drop = FALSE) +
    labs(title = "PAM50 subtype probability per tumor", x = NULL, y = "Probability (genefu subtype.proba)") +
    th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
               strip.background = element_blank(), strip.text = element_text(face = "bold"))
  save_both(p3, "PAM50_subtype_probability_by_sample", 12, 5)
}

# 8d. Correlation of each tumor with the five PAM50 centroids (always available from genefu)
cc <- grep("^cor_", colnames(calls), value = TRUE)
if (length(cc)) {
  cl <- do.call(rbind, lapply(cc, function(cn)
    data.frame(Key = calls$Key, Cluster = calls$Cluster, Centroid = sub("^cor_", "", cn), Correlation = calls[[cn]])))
  cl$Centroid <- ifelse(cl$Centroid %in% names(sub_names), sub_names[cl$Centroid], cl$Centroid)
  cl$Centroid <- factor(cl$Centroid, levels = rev(names(sub_cols)))
  cl$Key <- factor(cl$Key, levels = ord$Key); cl <- cl[!is.na(cl$Cluster), ]
  p4 <- ggplot(cl, aes(Key, Centroid, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.3) +
    facet_grid(~paste("Cluster", Cluster), scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(title = "Correlation with PAM50 centroids", x = NULL, y = NULL) +
    th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
               strip.background = element_blank(), strip.text = element_text(face = "bold"))
  save_both(p4, "PAM50_centroid_correlation_by_sample", 12, 4)
}

message("Done. Output in ", out_dir)