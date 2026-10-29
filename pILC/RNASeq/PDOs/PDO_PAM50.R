# PDO RNA-seq: PAM50 intrinsic subtyping of the 4 PDOs, placed within the pILC tumors.
# Method: genefu PAM50 (Parker 2009 centroids, Spearman correlation to the 5 centroids;
# the sample takes the subtype of the highest correlation) - same call as the tumor script.
#
# Design (n = 4 PDOs, no RNA clusters, no clinical data):
#  * ONE joint VST is built from the raw counts of the 31 tumors + 4 PDOs (blind, intercept only),
#    because separately fitted VST matrices should not be mixed.
#  * PRIMARY PDO call = genefu run on all 35 samples together.
#  * genefu calls can depend on the cohort (see the note in the tumor script). Instead of assuming how,
#    the script measures it:
#      (a) PDO-only run (n = 4) versus the joint run;
#      (b) tumor calls: original saved calls vs tumors-only run on the joint VST vs the joint run.
# Outputs (no dates in names): calls, probabilities, centroid correlations, concordance tables, figures (PNG + PDF).
# Run first: PDO_counts_to_VST_GENCODEv47.R. The tumor PAM50 script should have been run (for the concordance check).

suppressPackageStartupMessages({
  library(genefu)          # BiocManager::install("genefu")
  library(AnnotationDbi)
  library(org.Hs.eg.db)    # symbol -> Entrez ID (PAM50 centroids are matched by Entrez ID)
  library(DESeq2)
  library(ggplot2)
})
SEED <- 42
set.seed(SEED)

# ---------------------------------------------------------------- 1. Paths
root <- Sys.getenv("PILC_ROOT", "/ix1/alee/LO_LAB/Personal/Rahul/pILC")
rna  <- file.path(root, "RNASeq")
cnt_t_f <- file.path(rna, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv")
cnt_p_f <- file.path(rna, "PDO/VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv")
clus_f  <- file.path(rna, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
tum_calls_f <- file.path(rna, "PAM50/Stranded/PAM50_subtype_calls.csv")      # saved by the tumor PAM50 script (optional)
for (f in c(cnt_t_f, cnt_p_f, clus_f)) if (!file.exists(f)) stop("Input not found: ", f)
out_dir <- file.path(rna, "PDO/PAM50")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_COUNT <- 5      # same low-expression rule as the VST scripts: count >= 5 ...
MIN_FRAC  <- 0.10   # ... in at least 10% of samples

norm_id <- function(x) {                        # TP09.9_S34 / TP09-9 -> TP09-09 (tumors only)
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}

# ---------------------------------------------------------------- 2. Raw counts of tumors and PDOs
read_counts <- function(f) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("Gene" %in% colnames(d))
  m <- as.matrix(d[, setdiff(colnames(d), "Gene"), drop = FALSE])
  mode(m) <- "numeric"; rownames(m) <- make.unique(as.character(d$Gene))
  if (anyNA(m)) stop("Missing values in ", f)
  round(m)
}
cnt_t <- read_counts(cnt_t_f); colnames(cnt_t) <- norm_id(colnames(cnt_t))
cnt_p <- read_counts(cnt_p_f)
if (ncol(cnt_t) != 31) stop("Expected 31 tumors, found ", ncol(cnt_t))
if (ncol(cnt_p) != 4)  stop("Expected 4 PDOs, found ", ncol(cnt_p))
if (any(colnames(cnt_p) %in% colnames(cnt_t))) stop("A PDO name equals a tumor name.")

common <- intersect(rownames(cnt_t), rownames(cnt_p))
message("Genes: tumors ", nrow(cnt_t), ", PDOs ", nrow(cnt_p), ", shared ", length(common))
cnt <- cbind(cnt_t[common, , drop = FALSE], cnt_p[common, , drop = FALSE])
tum_keys <- colnames(cnt_t); pdo_keys <- colnames(cnt_p)

# ---------------------------------------------------------------- 3. Joint VST (blind to group)
keep <- rowSums(cnt >= MIN_COUNT) >= ceiling(ncol(cnt) * MIN_FRAC)
message("Genes kept after the low-expression filter: ", sum(keep), " of ", nrow(cnt))
storage.mode(cnt) <- "integer"
dds <- DESeqDataSetFromMatrix(countData = cnt[keep, ],
                              colData = data.frame(sample = colnames(cnt), row.names = colnames(cnt)),
                              design = ~ 1)
vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) { message("vst() failed (", conditionMessage(e), "); using varianceStabilizingTransformation()")
                                      varianceStabilizingTransformation(dds, blind = TRUE) })
vst <- assay(vsd)
saveRDS(vst, file.path(out_dir, "joint_tumor_PDO_VST_matrix.rds"))

# ---------------------------------------------------------------- 4. Gene symbols (one row per symbol)
ids <- rownames(vst)
if (mean(grepl("^ENSG", ids)) > 0.5) stop("Row names look like Ensembl IDs; the count files should carry gene symbols.")
sym <- sub("\\|.*$", "", ids)                    # "SYMBOL|ENSG..." (shared symbols) -> "SYMBOL"
okr <- !is.na(sym) & sym != ""
vst <- vst[okr, ]; sym <- sym[okr]
ord <- order(rowMeans(vst), decreasing = TRUE)   # for duplicated symbols keep the highest-expressed row
vst <- vst[ord, ]; sym <- sym[ord]
dup <- duplicated(sym)
vst <- vst[!dup, ]; rownames(vst) <- sym[!dup]

# ---------------------------------------------------------------- 5. Keep the PAM50 genes
data(pam50.robust, package = "genefu")
cmap <- pam50.robust$centroids.map
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

# ---------------------------------------------------------------- 6. PAM50 prediction (one function, three runs)
sub_names <- c(Basal = "Basal-like", Her2 = "HER2-enriched", LumA = "Luminal A",
               LumB = "Luminal B", Normal = "Normal-like")

run_pam50 <- function(samples) {
  pred <- molecular.subtyping(sbt.model = "pam50", data = t(expr50[, samples, drop = FALSE]),
                              annot = annot, do.mapping = TRUE)
  st <- unname(sub_names[as.character(pred$subtype)]); names(st) <- names(pred$subtype)
  out <- data.frame(Sample = names(st), PAM50 = factor(st, levels = unname(sub_names)),
                    row.names = names(st), stringsAsFactors = FALSE)
  if (is.matrix(pred$cor)) {                      # correlation of each sample to each centroid
    cm <- pred$cor[out$Sample, , drop = FALSE]
    colnames(cm) <- paste0("cor_", colnames(cm))
    out <- cbind(out, cm)
  }
  pr <- pred$subtype.proba                        # taken as genefu returns it (relative confidence of the call)
  if (!is.null(pr) && (is.matrix(pr) || is.data.frame(pr))) {
    pm <- as.matrix(pr)
    if (is.null(rownames(pm)) && nrow(pm) == nrow(out)) rownames(pm) <- out$Sample
    cn <- colnames(pm)
    colnames(pm) <- paste0("prob_", ifelse(cn %in% names(sub_names), unname(sub_names[cn]), cn))
    out <- cbind(out, pm[out$Sample, , drop = FALSE])
    message("subtype.proba found: row sums range ", paste(round(range(rowSums(pm)), 3), collapse = " - "))
  }
  out
}

res_all <- run_pam50(colnames(expr50))                                   # PRIMARY: 31 tumors + 4 PDOs
res_tum <- run_pam50(tum_keys)                                           # tumors only, joint VST
res_pdo <- tryCatch(run_pam50(pdo_keys),                                 # PDOs only (n = 4), sensitivity
                    error = function(e) { message("PDO-only run failed: ", conditionMessage(e)); NULL })

# ---------------------------------------------------------------- 7. Group labels and tables
clus <- read.csv(clus_f); clus$Key <- norm_id(clus$sample)
res_all$Group <- ifelse(res_all$Sample %in% pdo_keys, "PDO",
                        paste("Cluster", clus$cluster[match(res_all$Sample, clus$Key)]))
if (any(grepl("NA$", res_all$Group))) stop("No cluster for: ", paste(res_all$Sample[grepl("NA$", res_all$Group)], collapse = ", "))
res_all$Group <- factor(res_all$Group, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "PDO"))
write.csv(res_all, file.path(out_dir, "PAM50_calls_tumors_and_PDO_joint_run.csv"), row.names = FALSE)

# PDO table: primary call + PDO-only sensitivity
pdo_tab <- res_all[res_all$Sample %in% pdo_keys, , drop = FALSE]
colnames(pdo_tab)[colnames(pdo_tab) == "PAM50"] <- "PAM50_joint_run"
if (!is.null(res_pdo)) {
  pdo_tab$PAM50_PDO_only_run <- res_pdo$PAM50[match(pdo_tab$Sample, res_pdo$Sample)]
  pdo_tab$calls_agree <- as.character(pdo_tab$PAM50_joint_run) == as.character(pdo_tab$PAM50_PDO_only_run)
}
write.csv(pdo_tab, file.path(out_dir, "PDO_PAM50_subtype_calls.csv"), row.names = FALSE)
print(pdo_tab[, c("Sample", "PAM50_joint_run", intersect(c("PAM50_PDO_only_run", "calls_agree"), colnames(pdo_tab)))], row.names = FALSE)

# Tumor concordance: does adding the PDOs (or refitting the VST) change any tumor call?
conc <- data.frame(Sample = tum_keys,
                   tumors_only_joint_VST = as.character(res_tum$PAM50[match(tum_keys, res_tum$Sample)]),
                   joint_run_with_PDOs   = as.character(res_all$PAM50[match(tum_keys, res_all$Sample)]))
if (file.exists(tum_calls_f)) {
  orig <- read.csv(tum_calls_f)
  conc$original_tumor_script <- as.character(orig$PAM50[match(tum_keys, norm_id(orig$Key))])
  conc <- conc[, c("Sample", "original_tumor_script", "tumors_only_joint_VST", "joint_run_with_PDOs")]
} else message("Saved tumor calls not found (", tum_calls_f, "); the comparison with the original calls is skipped.")
write.csv(conc, file.path(out_dir, "PAM50_tumor_call_concordance.csv"), row.names = FALSE)

cmp_pair <- function(a, b, label) {
  ok <- !is.na(a) & !is.na(b)
  message(sprintf("%s: %d / %d tumors identical", label, sum(a[ok] == b[ok]), sum(ok)))
}
if ("original_tumor_script" %in% colnames(conc)) {
  cmp_pair(conc$original_tumor_script, conc$tumors_only_joint_VST, "original calls vs tumors-only on the joint VST (effect of the VST refit)")
  cmp_pair(conc$original_tumor_script, conc$joint_run_with_PDOs,   "original calls vs joint run with PDOs (total effect)")
}
cmp_pair(conc$tumors_only_joint_VST, conc$joint_run_with_PDOs,     "tumors-only vs joint run (effect of adding the PDOs)")

message("Subtype counts, joint run:"); print(table(res_all$Group, res_all$PAM50))

# ---------------------------------------------------------------- 8. Figures
sub_cols <- c("Luminal A" = "#1F4E9C", "Luminal B" = "#6BAED6", "HER2-enriched" = "#E377C2",
              "Basal-like" = "#D62728", "Normal-like" = "#2CA02C")
th <- theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))
save_both <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
}

# 8a. Subtype counts by group (tumors by RNA cluster, then PDOs)
p1 <- ggplot(res_all, aes(Group, fill = PAM50)) +
  geom_bar(position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = sub_cols, drop = FALSE) +
  labs(title = "PAM50 subtype: tumors by RNA cluster and PDOs", x = NULL, y = "Number of samples",
       caption = "Joint genefu run on 31 tumors + 4 PDOs.") + th
save_both(p1, "PDO_PAM50_subtype_by_group", 7, 5)

ord <- res_all[order(res_all$Group, res_all$PAM50, res_all$Sample), ]

# 8b. Subtype probability per sample (tumors grouped by cluster, PDOs last)
pc <- grep("^prob_", colnames(res_all), value = TRUE)
if (length(pc)) {
  pl <- do.call(rbind, lapply(pc, function(cn)
    data.frame(Sample = res_all$Sample, Group = res_all$Group, Probability = res_all[[cn]],
               Subtype = sub("^prob_", "", cn))))
  pl$Subtype <- factor(pl$Subtype, levels = names(sub_cols))
  pl$Sample  <- factor(pl$Sample, levels = ord$Sample)
  p3 <- ggplot(pl, aes(Sample, Probability, fill = Subtype)) +
    geom_col(width = 0.85, color = "black", linewidth = 0.2) +
    facet_grid(~Group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = sub_cols, drop = FALSE) +
    labs(title = "PAM50 subtype probability per sample", x = NULL, y = "Probability (genefu subtype.proba)") +
    th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
               strip.background = element_blank(), strip.text = element_text(face = "bold"))
  save_both(p3, "PDO_PAM50_subtype_probability_by_sample", 13, 5)
} else message("No subtype probabilities returned; probability plot skipped.")

# 8c. Correlation of each sample with the five PAM50 centroids
cc <- grep("^cor_", colnames(res_all), value = TRUE)
if (length(cc)) {
  cl <- do.call(rbind, lapply(cc, function(cn) {
    nm <- sub("^cor_", "", cn)
    data.frame(Sample = res_all$Sample, Group = res_all$Group,
               Centroid = ifelse(nm %in% names(sub_names), unname(sub_names[nm]), nm), Correlation = res_all[[cn]])
  }))
  cl$Centroid <- factor(cl$Centroid, levels = rev(names(sub_cols)))
  cl$Sample   <- factor(cl$Sample, levels = ord$Sample)
  p4 <- ggplot(cl, aes(Sample, Centroid, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.3) +
    facet_grid(~Group, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(title = "Correlation with PAM50 centroids", x = NULL, y = NULL) +
    th + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
               strip.background = element_blank(), strip.text = element_text(face = "bold"))
  save_both(p4, "PDO_PAM50_centroid_correlation_by_sample", 13, 4)
}

message("Done. Output in ", out_dir)
