############################################################
## DESeq2_ILC_vs_pILC.R
##
## PURPOSE
##   Differential expression, ILC vs pILC, on RAW HTSeq counts
##   (GENCODE v47 gene labels, all genes).
##
## SAMPLE GROUPS (preservation matters)
##   pILC_FF            31  fresh frozen
##   ILC_Bionano_FF      3  fresh frozen
##   ILC_Chelsea_FFPE   21  FFPE
##
## MODEL      design = ~ group  (reference = pILC_FF); three contrasts:
##   1. PRIMARY   ILC_Bionano_FF  vs pILC_FF   fresh frozen vs fresh frozen:
##                the only ILC vs pILC comparison NOT confounded with FFPE
##                (only 3 ILC samples -> low power; use effect sizes and
##                 concordance with contrast 2)
##   2. SUPPORT   ILC_Chelsea_FFPE vs pILC_FF  ILC vs pILC but FFPE-confounded
##   3. FFPE      ILC_Chelsea_FFPE vs ILC_Bionano_FF  ILC vs ILC: preservation
##                effect (and any Bionano-specific biology)
##   log2FC > 0 = higher in ILC (contrasts 1, 2) or higher in FFPE (contrast 3).
##   ComBat_seq-corrected counts and VST values are NOT used for the tests.
##
## GENE SCOPE (optional second argument)
##   all             all genes that pass the filter (default)
##   protein_coding  protein-coding genes only (sensitivity analysis: removes
##                   most non-coding / low-abundance genes that FFPE detects poorly)
##
## FFPE DIRECTION CHECK
##   "well_detected" = count >= 10 in at least 50% of the samples of EVERY group.
##   The FFPE ILC samples (Chelsea) are used only to check, for well-detected genes,
##   whether they change in the same direction as the fresh-frozen ILC samples.
##
## RUN
##   Rscript DESeq2_ILC_vs_pILC.R Reverse                  (or No)
##   Rscript DESeq2_ILC_vs_pILC.R Reverse protein_coding   (protein-coding only)
##   In RStudio the defaults below are used.
############################################################

SEED <- 42
set.seed(SEED)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)   # heatmap (install.packages("pheatmap") if missing)
})

############################################################
## 0. Settings
############################################################

args        <- commandArgs(trailingOnly = TRUE)
ilc_setting <- if (length(args) >= 1) args[1] else "Reverse"   # ILC counting setting: "Reverse" or "No"
stopifnot(ilc_setting %in% c("Reverse", "No"))
gene_scope  <- if (length(args) >= 2) args[2] else "all"       # "all" or "protein_coding"
stopifnot(gene_scope %in% c("all", "protein_coding"))

detect_count <- 10    # "well detected": count >= 10 ...
detect_frac  <- 0.50  # ... in at least 50% of the samples of every group

padj_cut  <- 0.05   # significance
lfc_cut   <- 1      # |log2 fold change| for the "up / down" call
min_count <- 5      # same filter rule as the VST pipeline
min_frac  <- 0.10   # count >= 5 in >= 10% of all samples
n_label   <- 10     # genes labelled per direction on the volcano
# Heatmap: ALL genes of the primary contrast (Bionano vs pILC) with
#   adjusted p-value (FDR) < heatmap_padj  and  |log2 fold change| >= heatmap_lfc
heatmap_padj <- padj_cut   # same cut-offs as the volcano, so the gene counts match
heatmap_lfc  <- lfc_cut    # (change padj_cut / lfc_cut above to change volcano and heatmap together)
max_heat     <- 200   # safety cap: if more genes pass, only the strongest by |Wald statistic| are shown

# wrap long plot text so it never runs off the figure (keeps manual line breaks)
wrap_txt <- function(x, width = 70) {
  paste(unlist(lapply(strsplit(x, "\n", fixed = TRUE)[[1]], function(l) strwrap(l, width = width))), collapse = "\n")
}

############################################################
## 1. Paths
############################################################

root     <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
ilc_dir  <- file.path(root, "4_Count_file_ILC", "2_Exploratory_strandedness", ilc_setting)
pilc_dir <- file.path(root, "4_Count_file_HTSeq_stranded", "VST_data_GENCODEv47")

out_dir <- file.path(root, "ILC_vs_pILC", "DESeq2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

need <- c(file.path(ilc_dir,  c("htseq_raw_counts_with_gene_names.rds", "sample_manifest.csv")),
          file.path(pilc_dir, c("htseq_raw_counts_with_gene_names.rds", "gene_annotation_from_GTF.tsv")))
if (any(!file.exists(need))) stop("Missing input files:\n", paste(need[!file.exists(need)], collapse = "\n"))

cat("ILC counting setting:", ilc_setting, "\n")

############################################################
## 2. Load raw counts and build the sample sheet
############################################################

ilc      <- readRDS(file.path(ilc_dir,  "htseq_raw_counts_with_gene_names.rds"))
pilc     <- readRDS(file.path(pilc_dir, "htseq_raw_counts_with_gene_names.rds"))
manifest <- fread(file.path(ilc_dir, "sample_manifest.csv"))
annot    <- fread(file.path(pilc_dir, "gene_annotation_from_GTF.tsv"))

dup <- intersect(colnames(ilc), colnames(pilc))
if (length(dup) > 0) stop("Sample names present in both cohorts: ", paste(dup, collapse = ", "))

ilc_cohort <- manifest$count_folder[match(colnames(ilc), manifest$sample)]   # Chelsea / Bionano
if (anyNA(ilc_cohort)) stop("Some ILC samples are missing from sample_manifest.csv")

genes <- intersect(rownames(ilc), rownames(pilc))
cat("Genes in ILC:", nrow(ilc), "| pILC:", nrow(pilc), "| shared:", length(genes), "\n")
if (length(genes) < 10000) stop("Too few shared genes - check that gene labels match.")

counts <- cbind(ilc[genes, , drop = FALSE], pilc[genes, , drop = FALSE])
storage.mode(counts) <- "integer"
stopifnot(!anyNA(counts))

grp <- c(ifelse(ilc_cohort == "Bionano", "ILC_Bionano_FF", "ILC_Chelsea_FFPE"), rep("pILC_FF", ncol(pilc)))
coldata <- data.frame(
  row.names = colnames(counts),
  group = factor(grp, levels = c("pILC_FF", "ILC_Bionano_FF", "ILC_Chelsea_FFPE"))
)
coldata$cohort       <- factor(ifelse(coldata$group == "pILC_FF", "pILC", "ILC"), levels = c("ILC", "pILC"))
coldata$preservation <- factor(ifelse(coldata$group == "ILC_Chelsea_FFPE", "FFPE", "Fresh frozen"),
                               levels = c("Fresh frozen", "FFPE"))
print(table(coldata$group))
stopifnot(sum(coldata$group == "ILC_Bionano_FF") == 3, sum(coldata$group == "ILC_Chelsea_FFPE") == 21,
          sum(coldata$group == "pILC_FF") == 31)

############################################################
## 3. Low-expression filter (same rule as the VST pipeline)
############################################################

if (gene_scope == "protein_coding") {
  pc_genes <- annot$Gene[annot$gene_biotype == "protein_coding"]
  counts   <- counts[rownames(counts) %in% pc_genes, , drop = FALSE]
  cat("Protein-coding genes only:", nrow(counts), "\n")
}

keep   <- rowSums(counts >= min_count) >= ceiling(ncol(counts) * min_frac)
counts <- counts[keep, , drop = FALSE]
cat("Genes after filter:", nrow(counts), "\n")

# genes that are well detected in EVERY group (needed for the FFPE direction check)
wd_vec <- Reduce(`&`, lapply(levels(coldata$group), function(g) {
  idx <- coldata$group == g
  rowSums(counts[, idx, drop = FALSE] >= detect_count) >= ceiling(detect_frac * sum(idx))
}))
names(wd_vec) <- rownames(counts)
cat("Well-detected genes in all three groups:", sum(wd_vec), "\n\n")

############################################################
## 4. DESeq2 (three contrasts from one model)
############################################################

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
dds <- DESeq(dds)
print(resultsNames(dds))

coef1 <- "group_ILC_Bionano_FF_vs_pILC_FF"
stopifnot(coef1 %in% resultsNames(dds))

res1 <- results(dds, contrast = c("group", "ILC_Bionano_FF",    "pILC_FF"),        alpha = padj_cut)  # primary
res2 <- results(dds, contrast = c("group", "ILC_Chelsea_FFPE",  "pILC_FF"),        alpha = padj_cut)  # FFPE-confounded
res3 <- results(dds, contrast = c("group", "ILC_Chelsea_FFPE",  "ILC_Bionano_FF"), alpha = padj_cut)  # FFPE effect

# NAIVE pooled comparison: all 24 ILC vs 31 pILC (ILC groups weighted by sample number).
# Shown ONLY to quantify how strongly the FFPE samples drive a pooled result.
n_ilc_b <- sum(coldata$group == "ILC_Bionano_FF")
n_ilc_c <- sum(coldata$group == "ILC_Chelsea_FFPE")
res4 <- results(dds, contrast = c(0, n_ilc_b, n_ilc_c) / (n_ilc_b + n_ilc_c), alpha = padj_cut)

# shrunken log2FC for the primary contrast (apeglm; ashr if apeglm is missing)
shr <- tryCatch(
  lfcShrink(dds, coef = coef1, type = "apeglm"),
  error = function(e) {
    message("apeglm failed (", conditionMessage(e), "); trying ashr")
    tryCatch(lfcShrink(dds, coef = coef1, res = res1, type = "ashr"),
             error = function(e2) { message("No shrinkage available: ", conditionMessage(e2)); NULL })
  }
)

############################################################
## 5. Results table (raw p and FDR-adjusted p for every gene)
############################################################

tab <- data.table(Gene = rownames(res1), as.data.frame(res1))
tab[, log2FC_shrunk := if (!is.null(shr)) shr[tab$Gene, "log2FoldChange"] else NA_real_]
tab[, log2FC_Chelsea_vs_pILC := res2[tab$Gene, "log2FoldChange"]]
tab[, padj_Chelsea_vs_pILC   := res2[tab$Gene, "padj"]]
tab[, log2FC_FFPE_effect     := res3[tab$Gene, "log2FoldChange"]]
tab[, padj_FFPE_effect       := res3[tab$Gene, "padj"]]
tab[, log2FC_pooled_24_vs_31 := res4[tab$Gene, "log2FoldChange"]]
tab[, padj_pooled_24_vs_31   := res4[tab$Gene, "padj"]]
tab <- merge(annot[, .(Gene, ensembl_id, gene_symbol, gene_biotype)], tab, by = "Gene", all.y = TRUE)

tab[, well_detected := unname(wd_vec[Gene])]
tab[, lfc_use := fifelse(!is.na(log2FC_shrunk), log2FC_shrunk, log2FoldChange)]
tab[, direction := fifelse(!is.na(padj) & padj < padj_cut & lfc_use >=  lfc_cut, "Up in ILC",
                   fifelse(!is.na(padj) & padj < padj_cut & lfc_use <= -lfc_cut, "Up in pILC", "NS"))]

# same direction in the FFPE ILC samples, and gene not strongly changed by preservation
tab[, concordant_sign := sign(log2FoldChange) == sign(log2FC_Chelsea_vs_pILC)]
tab[, preservation_sensitive := !is.na(padj_FFPE_effect) & padj_FFPE_effect < padj_cut &
                                abs(log2FC_FFPE_effect) >= lfc_cut]
tab[, robust := direction != "NS" & concordant_sign %in% TRUE & !preservation_sensitive &
                well_detected %in% TRUE]
setorder(tab, padj, na.last = TRUE)

tag <- paste0("ILC_vs_pILC_DESeq2_ILC", ilc_setting,
              if (gene_scope != "all") paste0("_", gene_scope) else "")
fwrite(tab[, !"lfc_use"], file.path(out_dir, paste0(tag, "_results.csv")))

############################################################
## 6. Console summary
############################################################

sig <- tab[direction != "NS"]
cat("\nPRIMARY contrast: ILC_Bionano_FF (n=3) vs pILC_FF (n=31)\n")
cat(sprintf("  tested genes: %d | padj < %.2f and |log2FC| >= %g: Up in ILC = %d, Up in pILC = %d\n",
            sum(!is.na(tab$padj)), padj_cut, lfc_cut,
            sum(tab$direction == "Up in ILC"), sum(tab$direction == "Up in pILC")))
cat(sprintf("  same direction in Chelsea (FFPE) vs pILC: %.1f%% of %d significant genes\n",
            100 * mean(sig$concordant_sign, na.rm = TRUE), nrow(sig)))
cat(sprintf("  flagged as preservation-sensitive (Chelsea vs Bionano): %d\n",
            sum(sig$preservation_sensitive)))
cat(sprintf("  log2FC correlation, Bionano vs pILC and Chelsea vs pILC (all tested genes): %.3f\n",
            cor(tab$log2FoldChange, tab$log2FC_Chelsea_vs_pILC, use = "complete.obs")))
cat(sprintf("  median shift (Chelsea vs pILC) - (Bionano vs pILC): %.2f log2 units (global FFPE offset)\n",
            median(tab$log2FC_Chelsea_vs_pILC - tab$log2FoldChange, na.rm = TRUE)))
cat(sprintf("SUPPORT contrast Chelsea (FFPE) vs pILC: padj < %.2f genes = %d (confounded with FFPE)\n",
            padj_cut, sum(tab$padj_Chelsea_vs_pILC < padj_cut, na.rm = TRUE)))

############################################################
## 6b. FFPE direction check (well-detected genes only)
##     Do the FFPE ILC samples move the same way as the fresh-frozen
##     ILC samples, each compared with pILC?
############################################################

cat(sprintf("\nPOOLED 24 ILC vs 31 pILC (naive, FFPE-dominated): padj < %.2f and |log2FC| >= %g: %d genes (primary Bionano vs pILC: %d)\n",
            padj_cut, lfc_cut,
            sum(tab$padj_pooled_24_vs_31 < padj_cut & abs(tab$log2FC_pooled_24_vs_31) >= lfc_cut, na.rm = TRUE),
            sum(tab$direction != "NS")))
cat(sprintf("  log2FC correlation of the pooled result with: Chelsea-only %.3f | Bionano-only %.3f\n",
            cor(tab$log2FC_pooled_24_vs_31, tab$log2FC_Chelsea_vs_pILC, use = "complete.obs"),
            cor(tab$log2FC_pooled_24_vs_31, tab$log2FoldChange, use = "complete.obs")))

wd <- tab[well_detected %in% TRUE & !is.na(log2FoldChange) & !is.na(log2FC_Chelsea_vs_pILC)]
dir_check <- function(d, label) {
  n <- nrow(d)
  k <- sum(sign(d$log2FoldChange) == sign(d$log2FC_Chelsea_vs_pILC))
  data.table(gene_set = label, n_genes = n,
             pct_same_direction = if (n > 0) 100 * k / n else NA_real_,
             binom_p_vs_50pct   = if (n > 0) binom.test(k, n, 0.5)$p.value else NA_real_,
             spearman_rho       = if (n > 2) cor(d$log2FoldChange, d$log2FC_Chelsea_vs_pILC, method = "spearman") else NA_real_)
}
sig_wd <- wd[direction != "NS"]
dc <- rbind(
  dir_check(wd,                                   "all well-detected genes"),
  dir_check(sig_wd,                               "primary-significant, well-detected"),
  dir_check(sig_wd[direction == "Up in ILC"],     "  of which Up in ILC"),
  dir_check(sig_wd[direction == "Up in pILC"],    "  of which Up in pILC"),
  dir_check(head(wd[order(-abs(stat))], 500),     "top 500 by |Wald statistic|")
)
fwrite(dc, file.path(out_dir, paste0(tag, "_direction_check_FFPE.csv")))
cat("\nFFPE direction check (well-detected genes):\n")
print(dc)

############################################################
## 7. Volcano plot of the PRIMARY contrast (PDF + PNG)
############################################################

vd <- tab[!is.na(padj)]
vd[, neglog10 := pmin(-log10(pmax(padj, 1e-300)), 300)]
lab <- rbind(head(vd[direction == "Up in ILC"],  n_label),
             head(vd[direction == "Up in pILC"], n_label))   # already sorted by padj

p <- ggplot(vd, aes(lfc_use, neglog10, colour = direction)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = lab, aes(label = Gene), size = 3, colour = "black",
                  max.overlaps = Inf, box.padding = 0.4, segment.size = 0.3) +
  scale_colour_manual(values = c("Up in ILC" = "#0072B2", "Up in pILC" = "#D55E00", "NS" = "grey75"),
                      breaks = c("Up in ILC", "Up in pILC"), name = NULL) +
  labs(x = sprintf("log2 fold change (ILC vs pILC%s)", if (!is.null(shr)) ", shrunken" else ""),
       y = "-log10 adjusted p-value",
       title = wrap_txt("ILC vs pILC, both fresh frozen", 45),
       subtitle = wrap_txt(sprintf("ILC_Bionano (n=3) vs pILC (n=31); padj < %.2f and |log2FC| >= %g:\n%d up in ILC, %d up in pILC (ILC counts: %s)",
                          padj_cut, lfc_cut, sum(vd$direction == "Up in ILC"),
                          sum(vd$direction == "Up in pILC"), ilc_setting), 70)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.subtitle = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(out_dir, paste0(tag, "_volcano.pdf")), p, width = 8, height = 7)
ggsave(file.path(out_dir, paste0(tag, "_volcano.png")), p, width = 8, height = 7, dpi = 300)

############################################################
## 8. Concordance plot: fresh frozen ILC vs FFPE ILC (both against pILC)
############################################################

sc <- wd   # well-detected genes only (see section 6b)
p_sc <- ggplot(sc, aes(log2FoldChange, log2FC_Chelsea_vs_pILC, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  scale_colour_manual(values = c("Up in ILC" = "#0072B2", "Up in pILC" = "#D55E00", "NS" = "grey75"),
                      breaks = c("Up in ILC", "Up in pILC"), name = "Primary call") +
  labs(x = "log2FC, ILC_Bionano (fresh frozen) vs pILC",
       y = "log2FC, ILC_Chelsea (FFPE) vs pILC",
       title = wrap_txt("Do FFPE and fresh-frozen ILC differ from pILC in the same way?", 45),
       subtitle = wrap_txt(sprintf("Well-detected genes (n = %d). Same direction: %.1f%% of all genes, %.1f%% of primary-significant genes; Spearman rho = %.2f\nDashed line = identical effect; offset from the line = preservation (FFPE) effect",
                          nrow(wd), dc$pct_same_direction[1], dc$pct_same_direction[2], dc$spearman_rho[1]), 68)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.subtitle = element_text(size = 10),
        legend.position = "bottom")
ggsave(file.path(out_dir, paste0(tag, "_FF_vs_FFPE_concordance.pdf")), p_sc, width = 7.5, height = 7)
ggsave(file.path(out_dir, paste0(tag, "_FF_vs_FFPE_concordance.png")), p_sc, width = 7.5, height = 7, dpi = 300)

############################################################
## 8b. PCA (PDF + PNG, scores and variance tables)
##     WHY: a quick, unbiased look at what separates the samples,
##          before trusting any DE list.
##     HOW: VST values, blind = TRUE (the design is NOT used, so the
##          groups cannot shape the plot), top n_top_pca most variable
##          genes, centred, not scaled; PC1 vs PC2 with % variance.
##     TWO PLOTS:
##       A. all samples          -> shows whether FFPE vs fresh frozen dominates
##       B. fresh frozen only    -> Bionano ILC vs pILC, the real primary comparison
##     No sample labels and no ellipses (n = 3 is too small for a group ellipse).
############################################################

n_top_pca <- 2000
vsd_blind <- assay(vst(dds, blind = TRUE))

pca_lab   <- c(pILC_FF = "pILC (fresh frozen)", ILC_Bionano_FF = "ILC Bionano (fresh frozen)",
               ILC_Chelsea_FFPE = "ILC Chelsea (FFPE)")
pca_col   <- setNames(c("#D55E00", "#009E73", "#0072B2"), pca_lab)
pca_shape <- setNames(c(15, 17, 16), pca_lab)   # square, triangle, circle

r2_of <- function(y, f) {                        # variance of a PC explained by a factor
  f <- droplevels(factor(f))
  if (nlevels(f) < 2) NA_real_ else summary(lm(y ~ f))$r.squared
}

run_pca <- function(mat, cd, stub, title, key_covariate) {
  cd  <- cd[colnames(mat), , drop = FALSE]
  top <- names(sort(apply(mat, 1, var), decreasing = TRUE))[seq_len(min(n_top_pca, nrow(mat)))]
  pc  <- prcomp(t(mat[top, ]), center = TRUE, scale. = FALSE)
  pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  k   <- 1:5

  scores <- data.frame(sample = rownames(pc$x), pc$x[, k], group = cd$group,
                       preservation = cd$preservation, row.names = NULL)
  var_tab <- data.frame(PC = paste0("PC", k), percent_variance = pct[k],
                        R2_with_group = sapply(k, function(i) r2_of(pc$x[, i], cd$group)),
                        R2_with_preservation = sapply(k, function(i) r2_of(pc$x[, i], cd$preservation)))
  fwrite(scores,  file.path(out_dir, paste0(tag, "_", stub, "_scores.csv")))
  fwrite(var_tab, file.path(out_dir, paste0(tag, "_", stub, "_variance.csv")))

  pd <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                   grp = factor(pca_lab[as.character(cd$group)], levels = pca_lab))
  r2_txt <- sprintf("PC1 R2 with %s = %.2f", key_covariate,
                    if (key_covariate == "preservation") var_tab$R2_with_preservation[1] else var_tab$R2_with_group[1])
  g <- ggplot(pd, aes(PC1, PC2, colour = grp, shape = grp)) +
    geom_point(size = 3.2, alpha = 0.85) +
    scale_colour_manual(values = pca_col, drop = TRUE, name = NULL) +
    scale_shape_manual(values = pca_shape, drop = TRUE, name = NULL) +
    labs(x = sprintf("PC1 (%.1f%% of variance)", pct[1]),
         y = sprintf("PC2 (%.1f%% of variance)", pct[2]),
         title = wrap_txt(title, 45),
         subtitle = wrap_txt(sprintf("%d samples; top %d most variable genes (VST, blind). %s",
                                     ncol(mat), length(top), r2_txt), 62)) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          plot.title.position = "plot", plot.subtitle = element_text(size = 10),
          legend.position = "bottom")
  ggsave(file.path(out_dir, paste0(tag, "_", stub, ".pdf")), g, width = 7.5, height = 6.5)
  ggsave(file.path(out_dir, paste0(tag, "_", stub, ".png")), g, width = 7.5, height = 6.5, dpi = 300)
  cat(sprintf("\nPCA %s: PC1 %.1f%%, PC2 %.1f%% | %s\n", stub, pct[1], pct[2], r2_txt))
  invisible(var_tab)
}

run_pca(vsd_blind, coldata, "PCA_all_samples",
        "PCA of all samples: what separates them?", "preservation")

ff <- rownames(coldata)[coldata$preservation == "Fresh frozen"]
run_pca(vsd_blind[, ff], coldata[ff, ], "PCA_fresh_frozen_only",
        "PCA, fresh frozen only: Bionano ILC vs pILC", "group")

############################################################
## 9. Heatmap (PDF + PNG) and its gene table
##    Genes  : ALL genes of the primary contrast (Bionano vs pILC) with
##             FDR < heatmap_padj and |log2FC| >= heatmap_lfc
##             (log2FC = shrunken value, same as the volcano plot)
##    Values : VST of all 55 samples, z-score per gene (clipped at +/- 3)
##    Columns: annotated by group (Chelsea FFPE / Bionano FF / pILC FF) and preservation
##    Rows   : side bar shows Up in ILC / Up in pILC
############################################################

vsd_mat <- assay(vst(dds, blind = FALSE))

hm <- tab[Gene %in% rownames(vsd_mat) & !is.na(padj) & padj < heatmap_padj & abs(lfc_use) >= heatmap_lfc]
cat(sprintf("\nGenes with FDR < %g and |log2FC| >= %g (Bionano vs pILC): %d\n",
            heatmap_padj, heatmap_lfc, nrow(hm)))
if (nrow(hm) < 2) stop("Fewer than 2 genes pass the heatmap cut-offs; relax heatmap_padj / heatmap_lfc.")
if (nrow(hm) > max_heat) {
  cat("More than max_heat genes pass; showing the", max_heat, "strongest by |Wald statistic|.\n")
  hm <- hm[order(-abs(stat))][seq_len(max_heat)]
}
hm[, direction_hm := fifelse(lfc_use > 0, "Up in ILC", "Up in pILC")]
top_genes <- hm$Gene
cat(sprintf("Heatmap genes: %d (Up in ILC = %d, Up in pILC = %d)\n",
            nrow(hm), sum(hm$direction_hm == "Up in ILC"), sum(hm$direction_hm == "Up in pILC")))

z <- t(scale(t(vsd_mat[top_genes, , drop = FALSE])))
z <- z[complete.cases(z), , drop = FALSE]
z[z >  3] <-  3
z[z < -3] <- -3

grp_lab <- c(ILC_Chelsea_FFPE = "ILC Chelsea (FFPE)", ILC_Bionano_FF = "ILC Bionano (fresh frozen)",
             pILC_FF = "pILC (fresh frozen)")
annotation_col <- data.frame(Group = factor(grp_lab[as.character(coldata[colnames(z), "group"])], levels = grp_lab),
                             Preservation = coldata[colnames(z), "preservation"],
                             row.names = colnames(z))
annotation_row <- data.frame(Direction = factor(hm[match(rownames(z), Gene), direction_hm],
                                                levels = c("Up in ILC", "Up in pILC")),
                             row.names = rownames(z))
annotation_colors <- list(Group = setNames(c("#0072B2", "#009E73", "#D55E00"), grp_lab),
                          Preservation = c("Fresh frozen" = "#CC79A7", FFPE = "#E69F00"),
                          Direction = c("Up in ILC" = "#56B4E9", "Up in pILC" = "#E69F00"))

ph <- pheatmap(z,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               breaks = seq(-3, 3, length.out = 101),
               annotation_col = annotation_col, annotation_row = annotation_row,
               annotation_colors = annotation_colors,
               clustering_method = "ward.D2",
               cluster_rows = TRUE, cluster_cols = TRUE,      # dendrograms: genes on the left, samples on top
               show_rownames = TRUE, show_colnames = FALSE,   # gene names on the right, no sample names
               fontsize_row = 8,
               treeheight_row = 40, treeheight_col = 40,
               main = wrap_txt(sprintf("%d genes, ILC vs pILC (FDR < %g, |log2FC| >= %g; VST, row z-score)",
                                       nrow(z), heatmap_padj, heatmap_lfc), 60),
               silent = TRUE)

hm_tag <- paste0(tag, "_heatmap_FDR", heatmap_padj, "_LFC", heatmap_lfc, "_n", nrow(z))
hm_h   <- max(6, 0.16 * nrow(z) + 4)    # figure height grows with the number of genes
pdf(file.path(out_dir, paste0(hm_tag, ".pdf")), width = 10, height = hm_h)
grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())
png(file.path(out_dir, paste0(hm_tag, ".png")), width = 10, height = hm_h, units = "in", res = 300)
grid::grid.newpage(); grid::grid.draw(ph$gtable); invisible(dev.off())

fwrite(hm[, .(Gene, ensembl_id, gene_symbol, gene_biotype, direction = direction_hm, baseMean,
              log2FoldChange, log2FC_shrunk, stat, pvalue, padj,
              log2FC_Chelsea_vs_pILC, padj_Chelsea_vs_pILC)],
       file.path(out_dir, paste0(hm_tag, "_genes.csv")))

cat("\nDone. Results in:", out_dir, "\n")