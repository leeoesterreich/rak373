###############################################################################
# PDO RNA-seq | cGAS-STING panel in the CANCER-CELL part of the tumors vs the 4 PDOs
#
# Question: PDOs are almost pure cancer cells, bulk tumors are mixtures (immune, stromal, fat).
#           Is a PDO "low" for a gene because its cancer cells express less, or because it lacks
#           the non-cancer cells that contribute to the bulk tumor value?
#
# Approach: InstaPrism (PBCA reference) already splits each bulk sample into cell types. Its
#           cell-type-specific expression (Z) for "Cancer Epithelial Cells" is read from the SAVED
#           deconvolution results (nothing is rerun) and compared with the PDOs' own cancer-cell Z.
#           The same panel (CGAS, STING1, TBK1, IRF3, STAT1, CCL5, CXCL10) and the same score are
#           used as in PDO_cGAS_STING.R. Both views are made side by side:
#             "Bulk"         = log2(TPM + 1) of the whole sample
#             "Cancer cells" = log2(InstaPrism cancer-epithelial Z + 1)
#           z-scores use the mean and SD of the 31 TUMORS (tumor scale), so a PDO value reads as
#           "SDs from the tumor mean". PDO percentile = % of the 31 tumors with a value <= the PDO.
#
# Caution : Z is a model ESTIMATE, as good as the PBCA reference; it is noisier for genes with low
#           cancer-cell expression and for tumors with a low cancer fraction (the fraction is shown).
#           Compare positions within the tumor range rather than exact differences.
#           This is mRNA, not pathway activity. n = 4 PDOs, tests are exploratory.
# Reads saved files only; it does NOT rerun InstaPrism and does not overwrite existing files.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(InstaPrism)    # reconstruct_Z_ct_initial()
  library(ggplot2)
  library(pheatmap)
})

SEED <- 42
set.seed(SEED)

# ------------------------------- 1. SETTINGS ---------------------------------
RNA_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
PDO_DIR <- file.path(RNA_DIR, "PDO")
INSTA_T <- file.path(RNA_DIR, "Immune_deconvolution_BC_Atlas/InstaPrism_PBCA_reference/Stranded")
INSTA_P <- file.path(PDO_DIR, "InstaPrism_PBCA_reference")

# Saved InstaPrism results. TUM_RDS = NULL: searched in INSTA_T (must find exactly one "*deconv*.rds");
# otherwise set the full path of the tumor result here.
TUM_RDS <- NULL
PDO_RDS <- file.path(INSTA_P, "PDO_InstaPrism_PBCA_deconv_result.rds")
TPM_T     <- file.path(RNA_DIR, "4_Count_file_HTSeq_stranded/Max_TPM/tpm_matrix_pILC.txt")
TPM_P     <- file.path(PDO_DIR, "Max_TPM/tpm_matrix_PDO.txt")
CLUSTER_F <- file.path(RNA_DIR, "NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
RANK_FILE <- file.path(INSTA_T, "Downstream_stats", "Cancer_epithelial_genes",
                       "pILC_InstaPrism_PBCA_cancer_epithelial_gene_ranking.csv")   # optional lookup

OUT_DIR <- file.path(INSTA_P, "Downstream", "cGAS_STING_cancer_cell")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CANCER_PATTERN <- "cancer.*epithelial"      # cell-type name searched in the InstaPrism result
Z_CLIP <- 3

PANEL <- list(CGAS    = c("CGAS", "MB21D1", "C6orf150"),
              STING1  = c("STING1", "TMEM173"),
              TBK1    = "TBK1",
              IRF3    = "IRF3",
              STAT1   = "STAT1",
              CCL5    = "CCL5",
              CXCL10  = "CXCL10")
CORE   <- c("CGAS", "STING1", "TBK1", "IRF3")
OUTPUT <- c("CCL5", "CXCL10", "STAT1")

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
group_levels   <- c(cluster_levels, "PDO")
group_colors   <- c(cluster_colors, PDO = "#333333")
pdo_cols       <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
SRC_LEVELS     <- c("Bulk (whole sample)", "Cancer cells (InstaPrism)")

# ------------------------------- 2. HELPERS ----------------------------------
need <- function(path, what) if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)
make_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x <- sub("([0-9])T$", "\\1", x)
  gsub("_0+([0-9])", "_\\1", x)
}
wrap_txt <- function(x, width) paste(strwrap(x, width), collapse = "\n")
fmt_p <- function(x) {
  if (is.na(x)) "NA" else if (x < 1e-4) format(x, digits = 3, scientific = TRUE)
  else format(signif(x, 3), scientific = FALSE)
}
save_plot <- function(p, stem, w, h) {
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
}
save_pheatmap <- function(draw_fun, stem, w, h) {
  png(paste0(stem, ".png"), width = w * 300, height = h * 300, res = 300); draw_fun(); invisible(dev.off())
  pdf(paste0(stem, ".pdf"), width = w, height = h);                        draw_fun(); invisible(dev.off())
}
pre <- function(x) file.path(OUT_DIR, paste0("PDO_cGAS_STING_cancer_cell_", x))
read_tpm <- function(path) {
  d <- fread(path, data.table = FALSE, check.names = FALSE)
  m <- as.matrix(d[, -1]); rownames(m) <- d[[1]]; storage.mode(m) <- "numeric"
  colnames(m) <- make_key(colnames(m))
  if (anyDuplicated(colnames(m))) stop("Duplicated sample keys in ", path)
  m
}

# cancer-cell Z (genes x samples) and cancer fraction from a saved InstaPrism result
cancer_Z <- function(rds) {
  res <- readRDS(rds)
  th  <- res@Post.ini.ct@theta                          # cell types x samples (or the reverse)
  ct_names <- if (any(grepl(CANCER_PATTERN, rownames(th), ignore.case = TRUE))) rownames(th) else colnames(th)
  ct <- grep(CANCER_PATTERN, ct_names, ignore.case = TRUE, value = TRUE)
  if (length(ct) != 1) stop("Could not find exactly one cancer-epithelial cell type in ", rds, ": ", paste(ct, collapse = ", "))
  frac <- if (ct %in% rownames(th)) th[ct, ] else th[, ct]
  Z <- reconstruct_Z_ct_initial(InstaPrism_obj = res, cell.type.of.interest = ct)
  smp <- names(frac)
  if (all(smp %in% rownames(Z))) Z <- t(Z)              # tutorial output: samples x genes -> genes x samples
  if (!all(smp %in% colnames(Z))) stop("Sample names of Z and theta do not match in ", rds)
  Z <- Z[, smp, drop = FALSE]
  colnames(Z) <- make_key(colnames(Z)); names(frac) <- colnames(Z)
  message("Cancer-cell Z from ", basename(rds), ": ", nrow(Z), " genes x ", ncol(Z), " samples (cell type: ", ct, ")")
  list(Z = pmax(Z, 0), frac = frac, ct = ct)
}

# ------------------------------- 3. LOAD -------------------------------------
if (is.null(TUM_RDS)) {
  cand <- list.files(INSTA_T, pattern = "deconv.*\\.rds$", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (length(cand) != 1) stop("Set TUM_RDS: found ", length(cand), " candidate tumor results in ", INSTA_T, "\n", paste(cand, collapse = "\n"))
  TUM_RDS <- cand
}
for (f in c(TUM_RDS, PDO_RDS, TPM_T, TPM_P, CLUSTER_F)) need(f, "Input file")
if (!exists("reconstruct_Z_ct_initial")) stop("InstaPrism is too old: reconstruct_Z_ct_initial() not found.")

ct_t <- cancer_Z(TUM_RDS)
ct_p <- cancer_Z(PDO_RDS)
tpm_t <- read_tpm(TPM_T); tpm_p <- read_tpm(TPM_P)
pdo_keys <- colnames(ct_p$Z)
if (length(pdo_keys) != 4) stop("Expected 4 PDOs, found ", length(pdo_keys))
pdo_label <- setNames(c("BO-046", "BO-114", "BO-189", "BO-268")[match(pdo_keys, make_key(c("BO-046", "BO-114", "BO-189", "BO-268")))], pdo_keys)
if (anyNA(pdo_label)) pdo_label <- setNames(pdo_keys, pdo_keys)

cl_raw <- fread(CLUSTER_F, data.table = FALSE)
s_hit <- grep("^sample$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
s_col <- if (length(s_hit) > 0) s_hit[1] else colnames(cl_raw)[1]
c_hit <- grep("^cluster$", colnames(cl_raw), ignore.case = TRUE, value = TRUE)
c_col <- if (length(c_hit) > 0) c_hit[1] else colnames(cl_raw)[ncol(cl_raw)]
cl <- data.frame(key = make_key(cl_raw[[s_col]]),
                 num = suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cl_raw[[c_col]])))))
if (anyNA(cl$num) || !all(cl$num %in% 1:3)) stop("Cluster values must be 1, 2 or 3 in ", CLUSTER_F)
tum_keys <- Reduce(intersect, list(colnames(ct_t$Z), colnames(tpm_t), cl$key))
message("Tumors available in Z, TPM and the cluster file: ", length(tum_keys))
if (length(tum_keys) < 25) stop("Too few matched tumors. Check the sample IDs.")
tum_cl <- setNames(paste("Cluster", cl$num[match(tum_keys, cl$key)]), tum_keys)

# panel genes present in all four matrices
gene_map <- rbindlist(lapply(names(PANEL), function(g) {
  ok <- PANEL[[g]][PANEL[[g]] %in% rownames(ct_t$Z) & PANEL[[g]] %in% rownames(ct_p$Z) &
                   PANEL[[g]] %in% rownames(tpm_t) & PANEL[[g]] %in% rownames(tpm_p)]
  data.table(Gene = g, Accepted_symbols = paste(PANEL[[g]], collapse = " / "),
             Symbol_used = ifelse(length(ok) == 0, NA_character_, ok[1]), Found = length(ok) > 0)
}))
fwrite(gene_map, pre("gene_map.tsv"), sep = "\t", quote = FALSE, na = "")
print(gene_map)
if (sum(gene_map$Found) < 5) stop("Fewer than 5 of the 7 panel genes are present in all matrices: check gene_map.tsv.")
found <- gene_map[Found == TRUE]

# optional: does the tumor cancer-cell specificity ranking contain these genes?
if (file.exists(RANK_FILE)) {
  rk <- fread(RANK_FILE, data.table = FALSE)
  if (all(c("gene", "rank") %in% colnames(rk))) {
    hit <- rk[match(found$Symbol_used, rk$gene), , drop = FALSE]
    rk_out <- data.table(Gene = found$Gene, Symbol_used = found$Symbol_used,
                         Rank_in_cancer_epithelial_ranking = hit$rank)
    fwrite(rk_out, pre("rank_in_cancer_epithelial_gene_ranking.tsv"), sep = "\t", quote = FALSE, na = "")
    message("Panel genes in the tumor cancer-epithelial ranking file (NA = not in the file):"); print(rk_out)
  }
}

# ------------------------------- 4. FEATURES (bulk and cancer cells) ---------
build <- function(src) {
  if (src == SRC_LEVELS[1]) {
    et <- log2(tpm_t[found$Symbol_used, tum_keys, drop = FALSE] + 1)
    ep <- log2(tpm_p[found$Symbol_used, pdo_keys, drop = FALSE] + 1)
  } else {
    et <- log2(ct_t$Z[found$Symbol_used, tum_keys, drop = FALSE] + 1)
    ep <- log2(ct_p$Z[found$Symbol_used, pdo_keys, drop = FALSE] + 1)
  }
  rownames(et) <- rownames(ep) <- found$Gene
  m   <- cbind(et, ep)
  mu  <- rowMeans(et); sdv <- apply(et, 1, sd); sdv[sdv == 0] <- NA
  zg  <- (m - mu) / sdv; zg[!is.finite(zg)] <- NA                    # a constant or missing gene cannot be z-scored
  mean_z <- function(g) { g <- intersect(g, rownames(zg)); if (!length(g)) rep(NA_real_, ncol(zg)) else colMeans(zg[g, , drop = FALSE], na.rm = TRUE) }
  scores <- rbind(cGAS_STING_score = mean_z(rownames(zg)), Sensing_core_score = mean_z(CORE), Downstream_output_score = mean_z(OUTPUT))
  feat <- rbind(scores, m)                                           # scores in z units, genes in log2 units
  fl <- melt(as.data.table(t(feat), keep.rownames = "Sample"), id.vars = "Sample", variable.name = "Feature", value.name = "Value")
  fl[, Feature := as.character(Feature)]
  fl[, Source := src]
  fl
}
fl <- rbindlist(lapply(SRC_LEVELS, build))
fl[, Group := ifelse(Sample %in% pdo_keys, "PDO", tum_cl[Sample])]
fl[, Source := factor(Source, levels = SRC_LEVELS)]
feature_names <- c("cGAS_STING_score", "Sensing_core_score", "Downstream_output_score", found$Gene)
feature_names <- feature_names[feature_names %in% fl$Feature]

# ------------------------------- 5. TABLES -----------------------------------
fwrite(dcast(fl, Sample + Group ~ Source + Feature, value.var = "Value"),
       pre("values_per_sample.tsv"), sep = "\t", quote = FALSE)

tv <- fl[Group != "PDO"]
pdo_tab <- fl[Group == "PDO", .(Source, Feature, PDO = pdo_label[Sample], Value)]
st <- tv[, .(Tumor_mean = mean(Value, na.rm = TRUE), Tumor_SD = sd(Value, na.rm = TRUE), Tumor_median = median(Value, na.rm = TRUE)), by = .(Source, Feature)]
pdo_tab <- merge(pdo_tab, st, by = c("Source", "Feature"))
pdo_tab[, Z_vs_tumors := (Value - Tumor_mean) / Tumor_SD]
pdo_tab[, Percentile_in_tumors := mapply(function(s, f, v) 100 * mean(tv[as.character(Source) == s & Feature == f, Value] <= v, na.rm = TRUE),
                                         as.character(Source), Feature, Value)]
pdo_tab[, Feature := factor(Feature, levels = feature_names)]
setorder(pdo_tab, Source, Feature, PDO)
pdo_tab[, Feature := as.character(Feature)]
fwrite(pdo_tab, pre("PDO_vs_tumors_per_PDO.tsv"), sep = "\t", quote = FALSE)
print(pdo_tab[Feature %in% c("cGAS_STING_score", "Sensing_core_score", "Downstream_output_score"),
              .(Source, Feature, PDO, Value = round(Value, 2), Percentile_in_tumors = round(Percentile_in_tumors))])

wt <- fl[, {
  a <- Value[Group == "PDO"]; b <- Value[Group != "PDO"]
  .(Median_PDO = median(a, na.rm = TRUE), Median_tumors = median(b, na.rm = TRUE),
    Wilcoxon_p = suppressWarnings(wilcox.test(a, b)$p.value))
}, by = .(Source, Feature)]
wt[, Median_diff := Median_PDO - Median_tumors]
wt[, Wilcoxon_FDR := p.adjust(Wilcoxon_p, method = "BH"), by = Source]
wt[, Feature := factor(Feature, levels = feature_names)]
setorder(wt, Source, Feature)
wt[, Feature := as.character(Feature)]
fwrite(wt, pre("PDO_vs_tumors_Wilcoxon_exploratory.tsv"), sep = "\t", quote = FALSE)

# cancer fraction per sample (how reliable the cancer-cell estimate is)
cf <- rbind(data.table(Sample = tum_keys, Group = tum_cl[tum_keys], Cancer_fraction = ct_t$frac[tum_keys]),
            data.table(Sample = pdo_label[pdo_keys], Group = "PDO", Cancer_fraction = ct_p$frac[pdo_keys]))
fwrite(cf, pre("cancer_fraction_per_sample.tsv"), sep = "\t", quote = FALSE)
if (any(ct_p$frac < 0.8)) message("Check: PDO cancer fraction below 0.8 in ", paste(pdo_label[ct_p$frac < 0.8], collapse = ", "))

# ------------------------------- 6. FIGURES ----------------------------------
pdo_pal <- setNames(pdo_cols[seq_along(pdo_keys)], pdo_label[pdo_keys])

# 6a. Where each PDO sits among the 31 tumors: bulk vs cancer cells
pp <- copy(pdo_tab); pp[, Feature := factor(Feature, levels = feature_names)]
p_pct <- ggplot(pp, aes(Feature, Percentile_in_tumors, colour = PDO)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = 90, fill = "grey92") +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey50") +
  geom_point(shape = 18, size = 3.6, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = pdo_pal, name = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 10, 25, 50, 75, 90, 100)) +
  facet_wrap(~ Source, ncol = 1) +
  labs(x = NULL, y = "PDO percentile among the 31 tumors",
       title = "cGAS-STING: where each PDO sits among the tumors",
       caption = wrap_txt(paste0("Percentile = % of the 31 tumors with a value at or below the PDO. Grey band: 10th to 90th percentile; dashed line: median. ",
                                 "Top: whole-sample log2(TPM + 1). Bottom: InstaPrism cancer-epithelial expression (model estimate). ",
                                 "A gene low in the top panel but not in the bottom panel is low mainly because the PDO lacks non-cancer cells."), 105)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"), plot.caption = element_text(size = 8, hjust = 0),
        plot.title.position = "plot", plot.caption.position = "plot")
save_plot(p_pct, pre("PDO_percentile_bulk_vs_cancer_cells"), 9, 8)

# 6b. Boxplots of the cancer-cell values (tumors by cluster, PDOs as diamonds)
lab <- wt[Source == SRC_LEVELS[2], .(Feature, Label = paste0(Feature, "\nPDOs vs tumors FDR = ", vapply(Wilcoxon_FDR, fmt_p, character(1))))]
d2 <- merge(fl[Source == SRC_LEVELS[2]], lab, by = "Feature")
d2[, Label := factor(Label, levels = lab[match(feature_names, Feature), Label])]
d2[, Group := factor(Group, levels = group_levels)]
d2[, PDOname := pdo_label[Sample]]
p_box <- ggplot() +
  geom_boxplot(data = d2[Group != "PDO"], aes(Group, Value, fill = Group), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(data = d2[Group != "PDO"], aes(Group, Value), width = 0.12, height = 0, size = 1.5, shape = 21, fill = "white") +
  geom_point(data = d2[Group == "PDO"], aes(Group, Value, colour = PDOname),
             position = position_jitter(width = 0.12, height = 0, seed = SEED), size = 3.2, shape = 18) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_colour_manual(values = pdo_pal, name = NULL) +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ Label, scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "Score (z-score units) or log2(cancer-cell Z + 1)",
       title = "cGAS-STING in the cancer-cell part of the tumors vs the PDOs",
       caption = wrap_txt("Cancer-cell expression is the InstaPrism (PBCA) cancer-epithelial estimate for tumors and PDOs. Scores are mean z-scores on the tumor scale. FDR: exploratory Wilcoxon, 4 PDOs vs 31 tumors.", 110)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 8, hjust = 0), plot.title.position = "plot", plot.caption.position = "plot")
save_plot(p_box, pre("boxplots_tumors_by_cluster_with_PDOs"), 12, 3.1 * ceiling(length(feature_names) / 4) + 1.3)

# 6c. Heatmap (cancer-cell values on the tumor scale): tumors by cluster, then PDOs; strip = cancer fraction
hm <- dcast(fl[Source == SRC_LEVELS[2]], Feature ~ Sample, value.var = "Value")
hm_m <- as.matrix(hm[, -1, with = FALSE]); rownames(hm_m) <- hm$Feature
hm_m <- hm_m[feature_names, , drop = FALSE]
mu <- rowMeans(hm_m[, tum_keys, drop = FALSE]); sdv <- apply(hm_m[, tum_keys, drop = FALSE], 1, sd)
z_all <- (hm_m - mu) / sdv
z_all[z_all >  Z_CLIP] <-  Z_CLIP
z_all[z_all < -Z_CLIP] <- -Z_CLIP
z_all <- z_all[complete.cases(z_all), , drop = FALSE]
ord_t <- unlist(lapply(cluster_levels, function(g) {
  cols <- tum_keys[tum_cl == g]
  if (length(cols) > 2) cols[hclust(dist(t(z_all[, cols, drop = FALSE])), "ward.D2")$order] else cols
}))
ord <- c(ord_t, pdo_keys)
gaps <- cumsum(table(factor(tum_cl[ord_t], levels = cluster_levels)))
frac_all <- c(ct_t$frac, ct_p$frac)
ann <- data.frame(Group = factor(c(tum_cl[ord_t], rep("PDO", length(pdo_keys))), levels = group_levels),
                  Cancer_fraction = frac_all[ord], row.names = ord)
colnames(z_all) <- ifelse(colnames(z_all) %in% pdo_keys, pdo_label[colnames(z_all)], colnames(z_all))
rownames(ann) <- ifelse(rownames(ann) %in% pdo_keys, pdo_label[rownames(ann)], rownames(ann))
ord_lab <- ifelse(ord %in% pdo_keys, pdo_label[ord], ord)
draw_hm <- function() {
  pheatmap(z_all[, ord_lab, drop = FALSE], color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           breaks = seq(-Z_CLIP, Z_CLIP, length.out = 102), annotation_col = ann,
           annotation_colors = list(Group = group_colors, Cancer_fraction = c("white", "#1B7837")),
           cluster_rows = TRUE, cluster_cols = FALSE, gaps_col = as.integer(gaps),
           clustering_method = "ward.D2", fontsize_row = 9, fontsize_col = 7, border_color = NA,
           main = "cGAS-STING, cancer-cell expression (SDs from the tumor mean; tumors by cluster, then PDOs)")
}
save_pheatmap(draw_hm, pre("heatmap_tumors_by_cluster_then_PDOs"), 11, 5.5)

message("Done. Outputs in: ", OUT_DIR)
