############################################################
## ILC_vs_pILC_cohort_diagnostics.R
##
## PURPOSE
##   Find out WHAT separates ILC from pILC in the PCA (PC1 ~ 50%).
##   Sample groups (preservation matters):
##     ILC_Chelsea_FFPE   21  FFPE
##     ILC_Bionano_FF      3  fresh frozen
##     pILC_FF            31  fresh frozen
##   Key question: do the 3 fresh-frozen ILC samples sit with pILC (then
##   FFPE drives the gap) or with the other ILC samples (then it is not
##   only preservation)?
##   NO batch correction here: RAW ILC counts + raw pILC counts.
##
## OUTPUT (ILC_vs_pILC/Diagnostics/, PDF + PNG + CSV, no dates)
##   1. Uncorrected PCA, 3 groups (all genes and protein-coding only)
##   2. Scree plot; bar colour = R2 with cohort and with preservation
##   3. PC vs covariate R2 heatmap
##   4. Sample QC by group (boxplots) + statistics table
##   5. PC1 top-loading genes and their biotypes
##   6. Biotype composition and log2FC distribution of the DESeq2 result
##      (only if DESeq2_ILC_vs_pILC.R was already run)
##
## RUN
##   Rscript ILC_vs_pILC_cohort_diagnostics.R Reverse     (or No)
############################################################

SEED <- 42
set.seed(SEED)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

############################################################
## 0. Settings and paths
############################################################

args        <- commandArgs(trailingOnly = TRUE)
ilc_setting <- if (length(args) >= 1) args[1] else "Reverse"
stopifnot(ilc_setting %in% c("Reverse", "No"))

n_top  <- 2000   # most variable genes used for the PCA
n_pc   <- 10     # principal components shown

root     <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
ilc_dir  <- file.path(root, "4_Count_file_ILC", "2_Exploratory_strandedness", ilc_setting)
pilc_dir <- file.path(root, "4_Count_file_HTSeq_stranded", "VST_data_GENCODEv47")
out_dir  <- file.path(root, "ILC_vs_pILC", "Diagnostics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

de_file <- file.path(root, "ILC_vs_pILC", "DESeq2",
                     paste0("ILC_vs_pILC_DESeq2_ILC", ilc_setting, "_results.csv"))
tag <- paste0("ILC_vs_pILC_ILC", ilc_setting)

grp_levels   <- c("ILC_Chelsea_FFPE", "ILC_Bionano_FF", "pILC_FF")
cols_group   <- c(ILC_Chelsea_FFPE = "#0072B2", ILC_Bionano_FF = "#009E73", pILC_FF = "#D55E00")
shapes_group <- c(ILC_Chelsea_FFPE = 16, ILC_Bionano_FF = 17, pILC_FF = 15)

save_fig <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(tag, "_", name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(tag, "_", name, ".png")), p, width = w, height = h, dpi = 300)
}

############################################################
## 1. Load RAW counts of both cohorts (no correction)
############################################################

need <- c(file.path(ilc_dir,  c("htseq_raw_counts_with_gene_names.rds", "sample_manifest.csv")),
          file.path(pilc_dir, c("htseq_raw_counts_with_gene_names.rds", "gene_annotation_from_GTF.tsv")))
if (any(!file.exists(need))) stop("Missing input files:\n", paste(need[!file.exists(need)], collapse = "\n"))

ilc      <- readRDS(file.path(ilc_dir,  "htseq_raw_counts_with_gene_names.rds"))
pilc     <- readRDS(file.path(pilc_dir, "htseq_raw_counts_with_gene_names.rds"))
manifest <- fread(file.path(ilc_dir, "sample_manifest.csv"))
annot    <- fread(file.path(pilc_dir, "gene_annotation_from_GTF.tsv"))

dup <- intersect(colnames(ilc), colnames(pilc))
if (length(dup) > 0) stop("Sample names present in both cohorts: ", paste(dup, collapse = ", "))

genes  <- intersect(rownames(ilc), rownames(pilc))
counts <- cbind(ilc[genes, , drop = FALSE], pilc[genes, , drop = FALSE])
storage.mode(counts) <- "integer"

ilc_cohort <- manifest$count_folder[match(colnames(ilc), manifest$sample)]
grp <- c(ifelse(ilc_cohort == "Bionano", "ILC_Bionano_FF", "ILC_Chelsea_FFPE"), rep("pILC_FF", ncol(pilc)))
coldata <- data.frame(row.names = colnames(counts), group = factor(grp, levels = grp_levels))
coldata$cohort       <- factor(ifelse(coldata$group == "pILC_FF", "pILC", "ILC"), levels = c("ILC", "pILC"))
coldata$preservation <- factor(ifelse(coldata$group == "ILC_Chelsea_FFPE", "FFPE", "Fresh frozen"),
                               levels = c("Fresh frozen", "FFPE"))
print(table(coldata$group))

group_labels <- setNames(sprintf("%s (n=%d)",
                         c("ILC Chelsea, FFPE", "ILC Bionano, fresh frozen", "pILC, fresh frozen"),
                         as.integer(table(coldata$group)[grp_levels])), grp_levels)
cat("Shared genes:", length(genes), "| samples:", ncol(counts), "\n")

############################################################
## 2. Gene biotype groups and per-sample QC covariates
############################################################

bt  <- annot$gene_biotype[match(rownames(counts), annot$Gene)]
sym <- annot$gene_symbol[match(rownames(counts), annot$Gene)]
bt[is.na(bt)] <- "unknown"

biotype_group <- function(b) {
  fcase(b == "protein_coding", "protein_coding",
        b == "lncRNA", "lncRNA",
        grepl("rRNA", b), "rRNA",
        grepl("pseudogene", b), "pseudogene",
        b %in% c("snoRNA", "snRNA", "scaRNA", "miRNA", "misc_RNA", "scRNA", "sRNA",
                 "vault_RNA", "Y_RNA", "ribozyme"), "small_RNA",
        default = "other")
}
bgrp    <- biotype_group(bt)
is_mito <- !is.na(sym) & grepl("^MT-", sym)

tot <- colSums(counts)
pct <- function(mask) 100 * colSums(counts[mask, , drop = FALSE]) / tot

qc <- data.table(
  sample = colnames(counts),
  group = coldata$group, cohort = coldata$cohort, preservation = coldata$preservation,
  library_size_log10  = log10(tot),
  genes_detected_ge10 = colSums(counts >= 10),
  pct_protein_coding  = pct(bgrp == "protein_coding"),
  pct_lncRNA          = pct(bgrp == "lncRNA"),
  pct_rRNA            = pct(bgrp == "rRNA"),
  pct_small_RNA       = pct(bgrp == "small_RNA"),
  pct_pseudogene      = pct(bgrp == "pseudogene"),
  pct_mito            = pct(is_mito)
)

############################################################
## 3. Filter + VST of all samples together (unsupervised)
############################################################

keep <- rowSums(counts >= 5) >= ceiling(ncol(counts) * 0.10)
dds  <- DESeqDataSetFromMatrix(counts[keep, ], coldata, design = ~ 1)
vmat <- assay(vst(dds, blind = TRUE))
bgrp_keep <- bgrp[keep]
cat("Genes after filter:", nrow(vmat), "\n")

run_pca <- function(m) {
  v   <- apply(m, 1, var)
  top <- names(sort(v, decreasing = TRUE))[1:min(n_top, length(v))]
  prcomp(t(m[top, ]), center = TRUE, scale. = FALSE)
}
pca_all <- run_pca(vmat)
pca_pc  <- run_pca(vmat[bgrp_keep == "protein_coding", ])

r2 <- function(y, x) summary(lm(y ~ x))$r.squared

############################################################
## 4. Figure 1: uncorrected PCA with the 3 groups
############################################################

plot_pca <- function(pca, title) {
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = coldata$group)
  ggplot(df, aes(PC1, PC2, colour = group, shape = group)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 3.6, alpha = 0.85) +
    scale_colour_manual(values = cols_group, labels = group_labels, name = NULL) +
    scale_shape_manual(values = shapes_group, labels = group_labels, name = NULL) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), title = title) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}
p_pca <- (plot_pca(pca_all, sprintf("All genes (top %d variable)", n_top)) +
          plot_pca(pca_pc,  sprintf("Protein-coding only (top %d variable)", n_top))) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Uncorrected PCA: where do the fresh-frozen ILC samples fall?",
                  subtitle = "If FFPE drives PC1, the 3 fresh-frozen ILC samples (green) sit with pILC (orange), not with the FFPE ILC samples (blue)") &
  theme(legend.position = "bottom")
save_fig(p_pca, "PCA_uncorrected_3groups", 13, 6.5)

pc_coords <- data.table(sample = rownames(pca_all$x), group = coldata$group,
                        PC1 = pca_all$x[, 1], PC2 = pca_all$x[, 2], PC3 = pca_all$x[, 3])
fwrite(pc_coords, file.path(out_dir, paste0(tag, "_PCA_uncorrected_coordinates.csv")))

############################################################
## 5. Figure 2: scree plot (bar colour = R2 with cohort / preservation)
############################################################

scree_one <- function(pca, label) {
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  rbind(
    data.table(geneset = label, colour_by = "R2 with cohort (ILC vs pILC)",
               PC = paste0("PC", 1:n_pc), var_pct = ve[1:n_pc],
               r2 = sapply(1:n_pc, function(k) r2(pca$x[, k], coldata$cohort))),
    data.table(geneset = label, colour_by = "R2 with preservation (FFPE vs fresh frozen)",
               PC = paste0("PC", 1:n_pc), var_pct = ve[1:n_pc],
               r2 = sapply(1:n_pc, function(k) r2(pca$x[, k], coldata$preservation)))
  )
}
scree <- rbind(scree_one(pca_all, sprintf("All genes (top %d)", n_top)),
               scree_one(pca_pc,  sprintf("Protein-coding only (top %d)", n_top)))
scree[, PC := factor(PC, levels = paste0("PC", 1:n_pc))]
fwrite(scree, file.path(out_dir, paste0(tag, "_scree_table.csv")))

p_scree <- ggplot(scree, aes(PC, var_pct, fill = r2)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sprintf("%.1f", var_pct)), vjust = -0.4, size = 3) +
  facet_grid(colour_by ~ geneset) +
  scale_fill_gradient(low = "grey85", high = "#B2182B", limits = c(0, 1), name = "R2") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% variance explained",
       title = "Scree plot (raw counts, no correction)",
       subtitle = "Red bar = that PC is mostly explained by the variable in the row label") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"), strip.background = element_rect(fill = "grey95"))
save_fig(p_scree, "scree", 12, 8)

############################################################
## 6. Figure 3: PC vs covariate association (R2)
############################################################

covs <- c("cohort", "preservation", "library_size_log10", "genes_detected_ge10", "pct_protein_coding",
          "pct_lncRNA", "pct_rRNA", "pct_small_RNA", "pct_pseudogene", "pct_mito")
r2_tab <- rbindlist(lapply(1:n_pc, function(k) {
  rbindlist(lapply(covs, function(cv) {
    data.table(PC = paste0("PC", k), covariate = cv, r2 = r2(pca_all$x[, k], qc[[cv]]))
  }))
}))
r2_tab[, PC := factor(PC, levels = paste0("PC", 1:n_pc))]
r2_tab[, covariate := factor(covariate, levels = rev(covs))]
fwrite(r2_tab, file.path(out_dir, paste0(tag, "_PC_covariate_R2.csv")))

p_r2 <- ggplot(r2_tab, aes(PC, covariate, fill = r2)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 3.2) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), name = "R2") +
  labs(x = NULL, y = NULL, title = "What explains each principal component?",
       subtitle = "R2 of PC score vs covariate (all genes). Association, not proof of cause: preservation and cohort differ in these metrics") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9))
save_fig(p_r2, "PC_covariate_R2", 11, 6)

############################################################
## 7. Figure 4: sample QC by group + statistics table
############################################################

qnum <- c("library_size_log10", "genes_detected_ge10", "pct_protein_coding", "pct_lncRNA",
          "pct_rRNA", "pct_small_RNA", "pct_pseudogene", "pct_mito")
ql <- melt(qc, id.vars = c("sample", "group"), measure.vars = qnum, variable.name = "metric")
ql[, metric := factor(metric, levels = qnum)]

pairs <- list(c("ILC_Chelsea_FFPE", "pILC_FF"),
              c("ILC_Bionano_FF",   "pILC_FF"),
              c("ILC_Chelsea_FFPE", "ILC_Bionano_FF"))
qstats <- rbindlist(lapply(qnum, function(m) {
  rbindlist(lapply(pairs, function(pr) {
    a <- qc[group == pr[1], get(m)]; b <- qc[group == pr[2], get(m)]
    data.table(metric = m, group_1 = pr[1], group_2 = pr[2],
               median_1 = median(a), median_2 = median(b), wilcox_p = wilcox.test(a, b)$p.value)
  }))
}))
qstats[, wilcox_FDR := p.adjust(wilcox_p, method = "BH")]
fwrite(qstats, file.path(out_dir, paste0(tag, "_sample_QC_by_group_stats.csv")))
fwrite(qc,     file.path(out_dir, paste0(tag, "_sample_QC_values.csv")))
print(qstats[group_1 == "ILC_Chelsea_FFPE" & group_2 == "pILC_FF"])

p_qc <- ggplot(ql, aes(group, value, colour = group)) +
  geom_boxplot(outlier.shape = NA, fill = NA, colour = "grey30", width = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = cols_group, guide = "none") +
  scale_x_discrete(labels = c("Chelsea\nFFPE", "Bionano\nFF", "pILC\nFF")) +
  labs(x = NULL, y = NULL, title = "Sample-level QC by group (raw counts)",
       subtitle = "FFPE typically shows lower gene detection and different RNA-biotype composition") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey95"))
save_fig(p_qc, "sample_QC_by_group", 13, 6.5)

############################################################
## 8. PC1 top-loading genes and their biotypes
############################################################

ld <- pca_all$rotation[, 1]
ld_tab <- data.table(Gene = names(ld), PC1_loading = as.numeric(ld))
ld_tab[, biotype := bt[match(Gene, rownames(counts))]]
ld_tab[, biotype_group := bgrp[match(Gene, rownames(counts))]]
ld_tab <- ld_tab[order(-abs(PC1_loading))]
fwrite(head(ld_tab, 200), file.path(out_dir, paste0(tag, "_PC1_top200_loadings.csv")))
cat("\nBiotype groups of the 200 strongest PC1 genes:\n")
print(table(head(ld_tab$biotype_group, 200)))

############################################################
## 9. DESeq2 result: log2FC distribution and biotype composition
############################################################

if (file.exists(de_file)) {
  de <- fread(de_file)[!is.na(padj)]
  de[, biotype_group := biotype_group(gene_biotype)]
  med_lfc <- median(de$log2FoldChange, na.rm = TRUE)
  cat(sprintf("\nMedian log2FC (primary contrast) over tested genes: %.2f (0 = no global shift)\n", med_lfc))

  comp <- de[, .N, by = .(direction, biotype_group)][, pct := 100 * N / sum(N), by = direction]
  fwrite(comp, file.path(out_dir, paste0(tag, "_DE_biotype_composition.csv")))
  comp[, direction := factor(direction, levels = c("Up in ILC", "Up in pILC", "NS"))]

  p_h <- ggplot(de, aes(log2FoldChange)) +
    geom_histogram(bins = 120, fill = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = med_lfc, colour = "#B2182B") +
    labs(x = "log2 fold change (primary contrast)", y = "Genes",
         title = "Global shift", subtitle = sprintf("Red line = median (%.2f)", med_lfc)) +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  p_b <- ggplot(comp, aes(direction, pct, fill = biotype_group)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2", name = "Biotype") +
    labs(x = NULL, y = "% of genes", title = "Biotype of DE genes") +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  save_fig(p_h + p_b + plot_layout(widths = c(1.3, 1)), "DE_logFC_and_biotype", 12, 5)
} else {
  message("DESeq2 results not found (", de_file, "); skipped section 9. Run DESeq2_ILC_vs_pILC.R first.")
}

############################################################
## 10. Console summary
############################################################

cat("\nR2 with preservation, all genes:   ",
    paste(sprintf("PC%d=%.2f", 1:3, sapply(1:3, function(k) r2(pca_all$x[, k], coldata$preservation))), collapse = "  "), "\n")
cat("R2 with cohort, all genes:         ",
    paste(sprintf("PC%d=%.2f", 1:3, sapply(1:3, function(k) r2(pca_all$x[, k], coldata$cohort))), collapse = "  "), "\n")
cat("R2 with cohort, protein-coding only:",
    paste(sprintf("PC%d=%.2f", 1:3, sapply(1:3, function(k) r2(pca_pc$x[, k], coldata$cohort))), collapse = "  "), "\n")
cat("\nDone. Results in:", out_dir, "\n")
