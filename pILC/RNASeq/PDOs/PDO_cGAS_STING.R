###############################################################################
# PDO RNA-seq | cGAS-STING signature of the 4 PDOs, placed within the pILC tumors
# Same signature and score as the tumor script (pILC cGAS-STING script):
#   7-gene panel of Chen et al., npj Breast Cancer 2024 (C6orf150 = CGAS, CCL5, CXCL10, IRF3, TBK1,
#   TMEM173 = STING1, STAT1). Score = mean of the per-gene z-scores.
#   Also split into a sensing/signalling core (CGAS, STING1, TBK1, IRF3) and a downstream output
#   (CCL5, CXCL10, STAT1; interferon / T-cell genes that rise with immune infiltration).
#
# Design (n = 4 PDOs, no RNA clusters, no TILs):
#  * One joint VST of the raw counts of 31 tumors + 4 PDOs (blind, intercept only), as in the PAM50 and
#    ssGSEA PDO scripts.
#  * z-scores use the mean and SD of the 31 TUMORS (tumor scale), so a PDO value reads as "SDs from the
#    tumor mean". Normalising 4 samples on their own would give an unstable scale.
#  * Statistics: only an EXPLORATORY PDOs-vs-tumors Wilcoxon test (4 vs 31, BH over features) and the
#    percentile of each PDO among the 31 tumors. The cluster tests, TILs-adjusted models, grade test and
#    clinical correlations of the tumor script cannot be repeated for the PDOs (no clusters, TILs or grade).
#  * PDOs contain no immune cells, so the downstream output genes are expected to be low in PDOs; the
#    sensing core is the better comparison.
# Positive value / z = higher in the PDOs than the tumor mean.
# Run first: PDO_counts_to_VST_GENCODEv47.R
#
# Reference : Chen M, Yu S, van der Sluis T, et al. cGAS-STING pathway expression correlates with genomic
#             instability and immune cell infiltration in breast cancer. npj Breast Cancer. 2024;10:1.
#             doi:10.1038/s41523-023-00609-z (author list shortened: complete it before citing)
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
ROOT    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR <- file.path(ROOT, "PDO")

CNT_T_FILE   <- file.path(ROOT, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv")
CNT_P_FILE   <- file.path(PDO_DIR, "VST_data_GENCODEv47/htseq_raw_counts_with_gene_names.csv")
CLUSTER_FILE <- file.path(ROOT, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
# optional checks (skipped if a file is missing)
SHARED     <- file.path(ROOT, "DESeq2_stranded", "Shared_inputs_for_pathway_analysis")
TUM_VST_F  <- file.path(SHARED, "VST_blind_ssGSEA_input.tsv")   # VST used by the tumor cGAS-STING script
META_F     <- file.path(SHARED, "Sample_metadata.tsv")          # cluster labels used by the tumor script
for (f in c(CNT_T_FILE, CNT_P_FILE, CLUSTER_FILE)) if (!file.exists(f)) stop("Missing input: ", f)

OUT_ROOT <- file.path(PDO_DIR, "Pathway_analysis", "cGAS_STING")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

MIN_COUNT <- 5      # same low-expression rule as the other joint-VST scripts: count >= 5 ...
MIN_FRAC  <- 0.10   # ... in at least 10% of samples
Z_CLIP    <- 3      # heatmap colour limit (SDs from the tumor mean)

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
               "Score used here: mean of per-gene z-scores (tumor mean and SD of the 31 tumors, joint tumor + PDO VST); the paper's exact formula was not checked.")

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
group_levels   <- c(cluster_levels, "PDO")
group_colors   <- c(cluster_colors, PDO = "#333333")
pdo_cols       <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")     # same PDO colours as the PCA

# ------------------------------- 2. HELPERS ----------------------------------
save_plot <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
}
save_pheatmap <- function(draw_fun, stem, w, h) {
  png(paste0(stem, ".png"), width = w * 300, height = h * 300, res = 300); draw_fun(); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = w, height = h);                        draw_fun(); invisible(dev.off())
}
fmt_p <- function(x) {
  if (is.na(x)) "NA" else if (x < 1e-4) format(x, digits = 3, scientific = TRUE)
  else format(signif(x, 3), scientific = FALSE)
}
wrap_txt <- function(x, width) paste(strwrap(x, width), collapse = "\n")   # wraps long titles/captions
norm_id <- function(x) {                        # TP09.9_S34 / TP09-9 -> TP09-09 (tumors only)
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}
read_counts <- function(f) {
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("Gene" %in% colnames(d))
  m <- as.matrix(d[, setdiff(colnames(d), "Gene"), drop = FALSE])
  mode(m) <- "numeric"; rownames(m) <- make.unique(as.character(d$Gene))
  if (anyNA(m)) stop("Missing values in ", f)
  round(m)
}
pre <- function(x) file.path(OUT_ROOT, paste0("PDO_cGAS_STING_", x))

# ------------------------------- 3. JOINT VST (tumors + PDOs) ----------------
cnt_t <- read_counts(CNT_T_FILE); colnames(cnt_t) <- norm_id(colnames(cnt_t))
cnt_p <- read_counts(CNT_P_FILE)
if (ncol(cnt_t) != 31) stop("Expected 31 tumors, found ", ncol(cnt_t))
if (ncol(cnt_p) != 4)  stop("Expected 4 PDOs, found ", ncol(cnt_p))
if (any(colnames(cnt_p) %in% colnames(cnt_t))) stop("A PDO name equals a tumor name.")
tum_keys <- colnames(cnt_t); pdo_keys <- colnames(cnt_p)

common <- intersect(rownames(cnt_t), rownames(cnt_p))
message("Genes: tumors ", nrow(cnt_t), ", PDOs ", nrow(cnt_p), ", shared ", length(common))
cnt <- cbind(cnt_t[common, , drop = FALSE], cnt_p[common, , drop = FALSE])

keep <- rowSums(cnt >= MIN_COUNT) >= ceiling(ncol(cnt) * MIN_FRAC)
message("Genes kept after the low-expression filter: ", sum(keep), " of ", nrow(cnt))
storage.mode(cnt) <- "integer"
dds <- DESeqDataSetFromMatrix(countData = cnt[keep, ],
                              colData = data.frame(sample = colnames(cnt), row.names = colnames(cnt)),
                              design = ~ 1)
vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) { message("vst() failed (", conditionMessage(e), "); using varianceStabilizingTransformation()")
                                      varianceStabilizingTransformation(dds, blind = TRUE) })
expr <- assay(vsd)

# One row per gene symbol ("SYMBOL|ENSG" -> "SYMBOL"; keep the highest-expressed row)
sym <- sub("\\|.*$", "", rownames(expr))
if (mean(grepl("^ENSG", sym)) > 0.5) stop("Row names look like Ensembl IDs; the count files should carry gene symbols.")
okr <- !is.na(sym) & sym != ""
expr <- expr[okr, ]; sym <- sym[okr]
o <- order(rowMeans(expr), decreasing = TRUE)
expr <- expr[o, ]; sym <- sym[o]
expr <- expr[!duplicated(sym), ]; rownames(expr) <- sym[!duplicated(sym)]

# ------------------------------- 4. CLUSTERS OF THE TUMORS -------------------
clus <- read.csv(CLUSTER_FILE); clus$Key <- norm_id(clus$sample)
cl <- paste("Cluster", clus$cluster[match(tum_keys, clus$Key)])
if (anyNA(clus$cluster[match(tum_keys, clus$Key)])) stop("No cluster for some tumors: check CLUSTER_FILE.")
cl <- setNames(factor(cl, levels = cluster_levels), tum_keys)
print(table(cl))

# Check: same cluster labels as the tumor cGAS-STING script (Sample_metadata.tsv)?
if (file.exists(META_F)) {
  meta <- fread(META_F)
  if (all(c("Sample", "Cluster") %in% names(meta))) {
    mk <- norm_id(as.character(meta$Sample)); mc <- as.character(meta$Cluster)
    mm <- match(tum_keys, mk)
    if (sum(!is.na(mm)) >= 10) {
      n_diff <- sum(mc[mm] != as.character(cl), na.rm = TRUE)
      message("Cluster labels vs Sample_metadata.tsv: ", n_diff, " of ", sum(!is.na(mm)), " tumors differ")
      if (n_diff > 0) warning("Cluster labels differ from Sample_metadata.tsv: 'Cluster 2' may not be the same group as in the tumor script.")
    }
  }
}

# ------------------------------- 5. SIGNATURE --------------------------------
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

all_keys  <- c(tum_keys, pdo_keys)
gene_expr <- expr[found$Symbol_used, all_keys, drop = FALSE]
rownames(gene_expr) <- found$Gene

# per-gene z-scores on the tumor scale (mean and SD of the 31 tumors)
tumor_z <- function(m) (m - rowMeans(m[, tum_keys, drop = FALSE])) / apply(m[, tum_keys, drop = FALSE], 1, sd)
zg <- tumor_z(gene_expr)
mean_z <- function(genes) {
  g <- intersect(genes, rownames(zg))
  if (length(g) == 0) rep(NA_real_, ncol(zg)) else colMeans(zg[g, , drop = FALSE])
}
scores <- rbind(cGAS_STING_score        = mean_z(rownames(zg)),
                Sensing_core_score      = mean_z(CORE),
                Downstream_output_score = mean_z(OUTPUT))
scores <- scores[!apply(scores, 1, function(x) all(is.na(x))), , drop = FALSE]
feat <- rbind(scores, gene_expr)                  # scores (z-score units) + genes (VST units)
feature_names <- rownames(feat)
main_scores   <- intersect(c("cGAS_STING_score", "Sensing_core_score", "Downstream_output_score"), feature_names)

# Check: tumor scores from the joint VST vs the tumor script's own VST (should be very similar)
if (file.exists(TUM_VST_F)) {
  tv <- fread(TUM_VST_F); tm <- as.matrix(tv[, -1, with = FALSE]); rownames(tm) <- tv[[1]]
  colnames(tm) <- norm_id(colnames(tm))
  cs <- intersect(colnames(tm), tum_keys)
  if (length(cs) >= 10 && all(found$Symbol_used %in% rownames(tm))) {
    zt <- t(scale(t(tm[found$Symbol_used, cs, drop = FALSE])))
    rho <- cor(colMeans(zt), scores["cGAS_STING_score", cs], method = "spearman")
    message(sprintf("Check: tumor cGAS-STING score, joint VST vs tumor-script VST: Spearman rho = %.3f (%d tumors)", rho, length(cs)))
  } else message("Check skipped: too few matching tumors or panel genes missing in the tumor VST file.")
}

# ------------------------------- 6. TABLES -----------------------------------
grp <- setNames(c(as.character(cl), rep("PDO", length(pdo_keys))), all_keys)
sample_tab <- data.table(Sample = all_keys, Group = grp[all_keys], t(feat[, all_keys]))
fwrite(sample_tab, pre("scores_and_gene_expression_per_sample.tsv"), sep = "\t", quote = FALSE)

fl <- melt(as.data.table(t(feat), keep.rownames = "Sample"), id.vars = "Sample",
           variable.name = "Feature", value.name = "Value")
fl[, Feature := as.character(Feature)]
fl[, Group := grp[Sample]]
tv_long <- fl[Group != "PDO"]

# each PDO against the 31 tumors
tum_stats <- tv_long[, .(Tumor_mean = mean(Value), Tumor_SD = sd(Value), Tumor_median = median(Value)), by = Feature]
cl_med <- dcast(tv_long[, .(M = median(Value)), by = .(Feature, Group)], Feature ~ Group, value.var = "M")
setnames(cl_med, cluster_levels, paste0("Median_", sub(" ", "_", cluster_levels)))
pdo_tab <- merge(fl[Group == "PDO", .(Feature, PDO = Sample, Value)], tum_stats, by = "Feature")
pdo_tab[, Z_vs_tumors := (Value - Tumor_mean) / Tumor_SD]
pdo_tab[, Percentile_in_tumors := mapply(function(f, v) 100 * mean(tv_long[Feature == f, Value] <= v), Feature, Value)]
pdo_tab <- merge(pdo_tab, cl_med, by = "Feature")
pdo_tab[, Feature := factor(Feature, levels = feature_names)]
setorder(pdo_tab, Feature, PDO)
pdo_tab[, Feature := as.character(Feature)]
fwrite(pdo_tab, pre("PDO_vs_tumors_per_PDO.tsv"), sep = "\t", quote = FALSE)
print(pdo_tab[Feature %in% main_scores, .(Feature, PDO, Value = round(Value, 2), Percentile_in_tumors = round(Percentile_in_tumors))])

# exploratory: 4 PDOs vs 31 tumors (Wilcoxon, BH over features)
wt <- fl[, {
  a <- Value[Group == "PDO"]; b <- Value[Group != "PDO"]
  .(Median_PDO = median(a), Median_tumors = median(b), Median_diff = median(a) - median(b),
    Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
}, by = Feature]
wt[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH")]
wt[, Feature := factor(Feature, levels = feature_names)]
setorder(wt, Feature)
wt[, Feature := as.character(Feature)]
fwrite(wt, pre("PDO_vs_tumors_Wilcoxon_exploratory.tsv"), sep = "\t", quote = FALSE)

# ------------------------------- 7. FIGURES ----------------------------------
# 7a. Boxplots by cluster (tumors) with the 4 PDOs as coloured diamonds
lab <- wt[, .(Feature, Label = paste0(Feature, "\nPDOs vs tumors FDR = ", vapply(Wilcoxon_FDR, fmt_p, character(1))))]
fl2 <- merge(fl, lab, by = "Feature")
fl2[, Label := factor(Label, levels = lab[match(feature_names, Feature), Label])]
fl2[, Group := factor(Group, levels = group_levels)]
pdo_pal <- setNames(pdo_cols[seq_along(pdo_keys)], pdo_keys)

box_plot <- function(dd, ncol_f, w, h, stem, title, caption) {
  p <- ggplot() +
    geom_boxplot(data = dd[Group != "PDO"], aes(Group, Value, fill = Group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(data = dd[Group != "PDO"], aes(Group, Value), width = 0.12, height = 0, size = 1.5, shape = 21, fill = "white") +
    geom_point(data = dd[Group == "PDO"], aes(Group, Value, colour = Sample),
               position = position_jitter(width = 0.12, height = 0, seed = SEED), size = 3.2, shape = 18) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    scale_colour_manual(values = pdo_pal, name = NULL) +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(~ Label, scales = "free_y", ncol = ncol_f) +
    labs(x = NULL, y = "Score (z-score units) or VST expression", title = title, caption = wrap_txt(caption, 110)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(size = 8),
          axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"),
          plot.caption = element_text(size = 8, hjust = 0), plot.title.position = "plot", plot.caption.position = "plot")
  save_plot(p, stem, w, h)
}
cap_txt <- paste0("Boxes: the 31 tumors by RNA cluster; diamonds: the 4 PDOs. Scores are mean z-scores on the tumor scale ",
                  "(tumor mean 0, SD 1). FDR: exploratory Wilcoxon, 4 PDOs vs 31 tumors. PDOs have no immune cells, so the ",
                  "downstream output genes (CCL5, CXCL10, STAT1) are expected to be low.")
box_plot(fl2, 4, 12, 3.1 * ceiling(length(feature_names) / 4) + 1.3, pre("boxplots_tumors_by_cluster_with_PDOs"),
         "cGAS-STING signature: PDOs within the tumors", cap_txt)
box_plot(fl2[Feature %in% main_scores], 3, 9, 4.6, pre("scores_tumors_by_cluster_with_PDOs"),
         "cGAS-STING scores: PDOs within the tumors", cap_txt)

# 7b. Heatmap: tumors grouped by cluster, then the PDOs (z on the tumor scale, clipped)
z_all <- tumor_z(feat)
z_all[z_all >  Z_CLIP] <-  Z_CLIP
z_all[z_all < -Z_CLIP] <- -Z_CLIP
z_all <- z_all[complete.cases(z_all), , drop = FALSE]
ord_t <- unlist(lapply(cluster_levels, function(g) {
  cols <- tum_keys[cl == g]
  if (length(cols) > 2) cols[hclust(dist(t(z_all[, cols, drop = FALSE])), "ward.D2")$order] else cols
}))
ord <- c(ord_t, pdo_keys)
gaps <- cumsum(c(table(factor(as.character(cl[ord_t]), levels = cluster_levels)), PDO = 0))[seq_along(cluster_levels)]
ann <- data.frame(Group = factor(grp[ord], levels = group_levels), row.names = ord)
draw_hm <- function() {
  pheatmap(z_all[, ord, drop = FALSE], color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           breaks = seq(-Z_CLIP, Z_CLIP, length.out = 102), annotation_col = ann,
           annotation_colors = list(Group = group_colors), cluster_rows = TRUE, cluster_cols = FALSE,
           gaps_col = as.integer(gaps), clustering_method = "ward.D2", fontsize_row = 9, fontsize_col = 7,
           border_color = NA, main = "cGAS-STING panel and scores (SDs from the tumor mean; tumors by cluster, then PDOs)")
}
save_pheatmap(draw_hm, pre("heatmap_tumors_by_cluster_then_PDOs"), 11, 5.5)

message("Done. Outputs in: ", OUT_ROOT)