###############################################################################
# PDO RNA-seq | ssGSEA (MSigDB Hallmark) of the 4 PDOs, placed within the pILC tumors
# Same call as the tumor script: GSVA::ssgseaParam(..., normalize = TRUE) + gsva().
#
# Design (n = 4 PDOs, no RNA clusters, no grade):
#  * ONE joint run on 31 tumors + 4 PDOs, so PDO and tumor scores share one scale
#    (normalising a 4-sample matrix on its own would give a different scale).
#  * The joint VST is built from the raw counts of both sets (blind, intercept only).
#  * Check: tumor scores from the joint run are compared with the saved tumor ssGSEA scores
#    (Spearman per pathway), to show whether adding the PDOs or refitting the VST changes them.
#  * Statistics: only an EXPLORATORY PDOs-vs-tumors Wilcoxon test (4 vs 31, BH over pathways).
#    PDOs also differ from tumors in culture conditions and in having no stromal/immune cells.
# Positive difference / z = higher in the PDOs.
# Run first: PDO_counts_to_VST_GENCODEv47.R   Output names carry no dates.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)          # needs >= 1.50 (ssgseaParam)
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
HALLMARK_GMT <- file.path(ROOT, "Pathway_analysis", "Dataset", "h.all.v2025.1.Hs.symbols.gmt")
TUM_SCORES_F <- file.path(ROOT, "Pathway_analysis/ssGSEA/ssGSEA_Hallmark_scores.csv")   # saved by the tumor script (optional check)
for (f in c(CNT_T_FILE, CNT_P_FILE, CLUSTER_FILE, HALLMARK_GMT))
  if (!file.exists(f)) stop("Missing input: ", f)
if (packageVersion("GSVA") < "1.50.0") stop("GSVA >= 1.50 is needed (ssgseaParam).")

OUT_SS <- file.path(PDO_DIR, "Pathway_analysis", "ssGSEA")
dir.create(OUT_SS, recursive = TRUE, showWarnings = FALSE)

FDR_CUT   <- 0.05
MIN_COUNT <- 5      # same low-expression rule as the VST scripts: count >= 5 ...
MIN_FRAC  <- 0.10   # ... in at least 10% of samples
N_SHOW    <- 12     # bar plot: max pathways per direction (significant ones only)
PER_PAGE  <- 12     # boxplots: panels per figure; more significant pathways -> extra figures (_page1, _page2, ...)
KEY_PATHWAYS <- c("ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "ANDROGEN_RESPONSE")
ADD_KEY_PATHWAYS <- FALSE   # TRUE: also add the key pathways to the boxplots even if not significant

cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")
cluster_levels <- names(cluster_colors)
group_levels   <- c(cluster_levels, "PDO")
group_colors   <- c(cluster_colors, PDO = "#333333")

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

# ------------------------------- 3. HALLMARK SETS ----------------------------
gmt      <- strsplit(readLines(HALLMARK_GMT), "\t")
hallmark <- setNames(lapply(gmt, function(x) unique(x[-(1:2)])),
                     sub("^HALLMARK_", "", vapply(gmt, `[`, character(1), 1)))
message("Hallmark sets: ", length(hallmark))

# ------------------------------- 4. JOINT VST (tumors + PDOs) ----------------
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
saveRDS(expr, file.path(OUT_SS, "joint_tumor_PDO_VST_matrix.rds"))

coverage <- mean(unique(unlist(hallmark)) %in% rownames(expr))
message(sprintf("ssGSEA input: %d genes x %d samples; %.1f%% of Hallmark genes present", nrow(expr), ncol(expr), 100 * coverage))
if (coverage < 0.8) warning("Less than 80% of Hallmark genes are in the ssGSEA input.")

# ------------------------------- 5. ssGSEA (joint run) -----------------------
ss_all <- gsva(ssgseaParam(exprData = expr, geneSets = hallmark, normalize = TRUE))
ss_t <- ss_all[, tum_keys, drop = FALSE]
ss_p <- ss_all[, pdo_keys, drop = FALSE]

fwrite(as.data.table(ss_p,   keep.rownames = "Pathway"), file.path(OUT_SS, "ssGSEA_Hallmark_scores_PDO.csv"))
fwrite(as.data.table(ss_all, keep.rownames = "Pathway"), file.path(OUT_SS, "ssGSEA_Hallmark_scores_tumors_and_PDO_joint_run.csv"))

# ------------------------------- 6. Cluster labels of the tumors -------------
clus <- read.csv(CLUSTER_FILE); clus$Key <- norm_id(clus$sample)
cl <- factor(paste("Cluster", clus$cluster[match(tum_keys, clus$Key)]), levels = cluster_levels)
if (anyNA(cl)) stop("No cluster for: ", paste(tum_keys[is.na(cl)], collapse = ", "))
print(table(cl))

# ------------------------------- 7. Check: joint-run tumor scores vs the saved tumor scores
if (file.exists(TUM_SCORES_F)) {
  orig <- fread(TUM_SCORES_F, data.table = FALSE, check.names = FALSE)
  rownames(orig) <- orig[[1]]; orig <- as.matrix(orig[, -1, drop = FALSE])
  colnames(orig) <- norm_id(colnames(orig))
  cs <- intersect(colnames(orig), tum_keys); ps <- intersect(rownames(orig), rownames(ss_t))
  if (length(cs) >= 10 && length(ps) >= 10) {
    rho <- sapply(ps, function(p) cor(orig[p, cs], ss_t[p, cs], method = "spearman"))
    chk <- data.table(Pathway = ps, Spearman_rho_saved_vs_joint_run = unname(rho))
    setorder(chk, Spearman_rho_saved_vs_joint_run)
    fwrite(chk, file.path(OUT_SS, "ssGSEA_tumor_scores_saved_vs_joint_run.csv"))
    message(sprintf("Tumor scores, saved vs joint run (%d tumors, %d pathways): median rho = %.3f, minimum = %.3f (%s)",
                    length(cs), length(ps), median(rho), min(rho), names(which.min(rho))))
  } else message("Saved tumor scores could not be matched to the tumors; check skipped.")
} else message("Saved tumor scores not found (", TUM_SCORES_F, "); check skipped.")

# ------------------------------- 8. PDOs versus tumors -----------------------
tum_mean <- rowMeans(ss_t); tum_sd <- apply(ss_t, 1, sd)
tum_med  <- apply(ss_t, 1, median)
q        <- t(apply(ss_t, 1, quantile, probs = c(0.25, 0.75)))
cl_med   <- sapply(cluster_levels, function(k) apply(ss_t[, cl == k, drop = FALSE], 1, median))
colnames(cl_med) <- paste0("tumor_median_", sub(" ", "_", cluster_levels))
pdo_mean <- rowMeans(ss_p)

p_w <- sapply(rownames(ss_t), function(p) suppressWarnings(wilcox.test(ss_p[p, ], ss_t[p, ], exact = FALSE)$p.value))
cmp <- data.frame(Pathway = rownames(ss_t), ss_p,
                  PDO_mean = pdo_mean, tumor_median_all = tum_med, tumor_Q1 = q[, 1], tumor_Q3 = q[, 2],
                  cl_med, check.names = FALSE, row.names = NULL)
cmp$PDO_mean_z_vs_tumors      <- (pdo_mean - tum_mean) / tum_sd
cmp$PDO_percentile_in_tumors  <- sapply(rownames(ss_t), function(p) 100 * mean(ss_t[p, ] <= pdo_mean[p]))
cmp$PDO_mean_minus_tumor_median <- pdo_mean - tum_med
cmp$Wilcoxon_p   <- p_w
cmp$Wilcoxon_FDR <- p.adjust(p_w, method = "BH")
cmp$Higher_in    <- ifelse(cmp$PDO_mean_z_vs_tumors > 0, "PDO", "tumors")
cmp <- cmp[order(-abs(cmp$PDO_mean_z_vs_tumors)), ]
write.csv(cmp, file.path(OUT_SS, "ssGSEA_Hallmark_PDO_vs_tumors.csv"), row.names = FALSE)
message("Pathways with Wilcoxon FDR < ", FDR_CUT, " (PDOs vs tumors, exploratory): ", sum(cmp$Wilcoxon_FDR < FDR_CUT))

# ------------------------------- 9. Heatmap: tumors by cluster, then PDOs ----
ord  <- order(cl)
mat  <- cbind(ss_t[, ord, drop = FALSE], ss_p)
z    <- (mat - tum_mean) / tum_sd                       # SDs from the tumor mean
z[z >  3] <-  3
z[z < -3] <- -3
grp  <- factor(c(as.character(cl[ord]), rep("PDO", ncol(ss_p))), levels = group_levels)
ann  <- data.frame(Group = grp, row.names = colnames(mat))
n_cl <- as.integer(table(cl))
draw_ss_heatmap <- function() {
  pheatmap(z, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           annotation_col = ann, annotation_colors = list(Group = group_colors),
           cluster_rows = TRUE, cluster_cols = FALSE,
           gaps_col = cumsum(n_cl),                            # gaps after every cluster, incl. before the PDOs
           fontsize_row = 8, fontsize_col = 7, border_color = NA,
           main = "ssGSEA Hallmark scores (SDs from the tumor mean; tumors by cluster, then PDOs)")
}
save_pheatmap(draw_ss_heatmap, file.path(OUT_SS, "ssGSEA_Hallmark_heatmap_PDO_vs_tumors"), 11, 9)

# ------------------------------- 9b. Top / bottom pathways by mean score across the 4 PDOs
# Ranking of ABSOLUTE ssGSEA scores of the PDOs (descriptive, no test, no tumor comparison).
N_TOPBOT <- 5
pdo_cols <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[seq_along(pdo_keys)], pdo_keys)
ord_p    <- names(sort(rowMeans(ss_p), decreasing = TRUE))
sel_top  <- head(ord_p, N_TOPBOT)
sel_bot  <- tail(ord_p, N_TOPBOT)                              # also in descending order
panel_lv <- c(paste("Top", N_TOPBOT), paste("Bottom", N_TOPBOT))
tb_long <- rbindlist(lapply(c(sel_top, sel_bot), function(p) data.table(
  Pathway = p, PDO = pdo_keys, Score = as.numeric(ss_p[p, pdo_keys]),
  Panel = if (p %in% sel_top) panel_lv[1] else panel_lv[2])))
tb_long[, Panel := factor(Panel, levels = panel_lv)]
tb_long[, Pathway_lab := factor(gsub("_", " ", Pathway), levels = rev(gsub("_", " ", c(sel_top, sel_bot))))]
tb_sum <- tb_long[, .(Mean = mean(Score), Min = min(Score), Max = max(Score)), by = .(Panel, Pathway, Pathway_lab)]

p_tb <- ggplot(tb_long, aes(Score, Pathway_lab)) +
  geom_segment(data = tb_sum, aes(x = Min, xend = Max, y = Pathway_lab, yend = Pathway_lab),
               colour = "grey65", linewidth = 1, inherit.aes = FALSE) +
  geom_point(aes(colour = PDO), size = 3.5, alpha = 0.9) +
  geom_point(data = tb_sum, aes(x = Mean, y = Pathway_lab), shape = 18, size = 5, colour = "black", inherit.aes = FALSE) +
  facet_grid(Panel ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = pdo_cols, name = "PDO") +
  labs(x = "ssGSEA score", y = NULL,
       title = paste0("Top ", N_TOPBOT, " and bottom ", N_TOPBOT, " Hallmark pathways\nby mean ssGSEA score"),
       subtitle = paste0("Dots: individual PDOs. Black diamond: mean of the ", length(pdo_keys), " PDOs.\nGrey line: min to max."),
       caption = paste0("Descriptive ranking only: one sample per PDO (n = ", length(pdo_keys), ", no replicates), no statistical test.")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(), plot.caption = element_text(size = 9, hjust = 0))
save_plot(p_tb, file.path(OUT_SS, "ssGSEA_Hallmark_PDO_top_bottom_pathways"), 9, 6.5)

tb_tab <- merge(tb_sum[, .(Panel, Pathway, Mean, Min, Max)],
                dcast(tb_long, Panel + Pathway ~ PDO, value.var = "Score"), by = c("Panel", "Pathway"))
setorder(tb_tab, -Mean)
fwrite(tb_tab, file.path(OUT_SS, "ssGSEA_Hallmark_PDO_top_bottom_pathways.csv"))

# ------------------------------- 9c. Where the 4 PDOs differ most (ranked by range across the PDOs)
# Range = max - min of the ssGSEA score over the 4 PDOs. Descriptive only (one sample per PDO, no test).
N_DIFF <- 10
rng_all <- apply(ss_p, 1, function(x) max(x) - min(x))
sel_d   <- head(names(sort(rng_all, decreasing = TRUE)), N_DIFF)        # largest range first
d_long <- rbindlist(lapply(sel_d, function(p) data.table(Pathway = p, PDO = pdo_keys, Score = as.numeric(ss_p[p, pdo_keys]))))
d_long[, Pathway_lab := factor(gsub("_", " ", Pathway), levels = rev(gsub("_", " ", sel_d)))]   # largest range at the top
d_sum <- d_long[, .(Mean = mean(Score), Min = min(Score), Max = max(Score), Range = max(Score) - min(Score),
                    Highest_PDO = PDO[which.max(Score)], Lowest_PDO = PDO[which.min(Score)]),
                by = .(Pathway, Pathway_lab)]

p_diff <- ggplot(d_long, aes(Score, Pathway_lab)) +
  geom_segment(data = d_sum, aes(x = Min, xend = Max, y = Pathway_lab, yend = Pathway_lab),
               colour = "grey65", linewidth = 1, inherit.aes = FALSE) +
  geom_point(aes(colour = PDO), size = 3.5, alpha = 0.9) +
  geom_point(data = d_sum, aes(x = Mean, y = Pathway_lab), shape = 18, size = 5, colour = "black", inherit.aes = FALSE) +
  scale_colour_manual(values = pdo_cols, name = "PDO") +
  labs(x = "ssGSEA score", y = NULL,
       title = "Where the four PDOs differ most",
       caption = paste0("Top ", N_DIFF, " Hallmark pathways ranked by range (max minus min) across the ", length(pdo_keys),
                        " PDOs. Dots: PDOs; diamond: mean. Descriptive only (one sample per PDO, no test).")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 9, hjust = 0))
save_plot(p_diff, file.path(OUT_SS, "ssGSEA_Hallmark_PDO_largest_range_pathways"), 8.5, 5.5)

d_tab <- merge(d_sum[, .(Pathway, Range, Mean, Min, Max, Highest_PDO, Lowest_PDO)],
               dcast(d_long, Pathway ~ PDO, value.var = "Score"), by = "Pathway")
setorder(d_tab, -Range)
fwrite(d_tab, file.path(OUT_SS, "ssGSEA_Hallmark_PDO_largest_range_pathways.csv"))
# All 50 pathways, for reference
fwrite(data.table(Pathway = names(rng_all), Range = rng_all, ss_p)[order(-Range)],
       file.path(OUT_SS, "ssGSEA_Hallmark_PDO_range_all_pathways.csv"))

# ------------------------------- 10. Bar plot: pathways most different in the PDOs
sig  <- cmp[cmp$Wilcoxon_FDR < FDR_CUT, ]                 # significant pathways only (already sorted by |z|)
if (nrow(sig) == 0) stop("No pathway has Wilcoxon FDR < ", FDR_CUT, "; nothing to plot (see the CSV).")
up   <- head(sig$Pathway[order(-sig$PDO_mean_z_vs_tumors)], N_SHOW)
down <- head(sig$Pathway[order( sig$PDO_mean_z_vs_tumors)], N_SHOW)
bar_df <- sig[sig$Pathway %in% c(up, down), ]
bar_df$Direction <- ifelse(bar_df$PDO_mean_z_vs_tumors > 0, "Higher in PDOs", "Lower in PDOs")
bar_df$Pathway <- factor(bar_df$Pathway, levels = bar_df$Pathway[order(bar_df$PDO_mean_z_vs_tumors)])
p_bar <- ggplot(bar_df, aes(PDO_mean_z_vs_tumors, Pathway, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = c("Higher in PDOs" = "#B2182B", "Lower in PDOs" = "#2166AC"), name = NULL) +
  labs(x = "Mean PDO ssGSEA score (SD from the tumor mean)", y = NULL,
       title = "Hallmark pathways significantly different in the PDOs",
       caption = paste0("Only pathways with FDR < ", FDR_CUT, " (PDOs vs tumors, Wilcoxon, BH); up to ", N_SHOW,
                        " per direction. Exploratory; PDOs also differ from tumors in culture conditions and cell composition.")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
save_plot(p_bar, file.path(OUT_SS, "ssGSEA_Hallmark_PDO_top_pathways_barplot"), 8, 7)

# ------------------------------- 11. Boxplots: tumors by cluster, PDOs overlaid
# Only significant pathways (FDR < FDR_CUT), ordered by distance from the tumor mean (|z|).
# Many pathways can tie at the lowest possible FDR with 4 vs 31 samples, so |z| is used for ordering.
show_all <- sig$Pathway
if (ADD_KEY_PATHWAYS) show_all <- unique(c(show_all, intersect(KEY_PATHWAYS, cmp$Pathway)))
n_pages  <- ceiling(length(show_all) / PER_PAGE)
message(length(show_all), " pathways in the boxplots (", n_pages, " figure(s) of up to ", PER_PAGE, " panels).")
pdo_cols <- setNames(c("#D55E00", "#0072B2", "#009E73", "#CC79A7")[seq_along(pdo_keys)], pdo_keys)

for (pg in seq_len(n_pages)) {
  show <- show_all[((pg - 1) * PER_PAGE + 1):min(pg * PER_PAGE, length(show_all))]
  lab  <- setNames(paste0(show, "\nPDO vs tumors FDR = ",
                          vapply(cmp$Wilcoxon_FDR[match(show, cmp$Pathway)], fmt_p, character(1))), show)

  tum_long <- data.table(Pathway = rep(show, times = length(tum_keys)),
                         Sample  = rep(tum_keys, each = length(show)),
                         Score   = as.vector(ss_t[show, , drop = FALSE]))
  tum_long[, Group := as.character(cl)[match(Sample, tum_keys)]]
  pdo_long <- data.table(Pathway = rep(show, times = length(pdo_keys)),
                         Sample  = rep(pdo_keys, each = length(show)),
                         Score   = as.vector(ss_p[show, , drop = FALSE]), Group = "PDO")
  tum_long[, Panel := factor(lab[Pathway], levels = unname(lab))]
  pdo_long[, Panel := factor(lab[Pathway], levels = unname(lab))]

  p_box <- ggplot() +
    geom_boxplot(data = tum_long, aes(Group, Score, fill = Group), outlier.shape = NA, alpha = 0.5, colour = "black") +
    geom_point(data = tum_long, aes(Group, Score), colour = "grey30", size = 1.2,
               position = position_jitter(width = 0.15, height = 0, seed = SEED)) +
    geom_point(data = pdo_long, aes(Group, Score, colour = Sample), size = 3,
               position = position_jitter(width = 0.1, height = 0, seed = SEED)) +
    facet_wrap(~ Panel, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = cluster_colors, guide = "none") +
    scale_colour_manual(values = pdo_cols, name = "PDO") +
    scale_x_discrete(limits = group_levels) +
    labs(x = NULL, y = "ssGSEA Enrichment Score",
         title = paste0("Hallmark ssGSEA, significant pathways: tumors by cluster and PDOs",
                        if (n_pages > 1) paste0(" (page ", pg, " of ", n_pages, ")") else ""),
         caption = paste0("Only pathways with FDR < ", FDR_CUT, " (PDOs vs all tumors, Wilcoxon, BH over ", nrow(cmp),
                          " pathways), exploratory. Boxes: tumors by RNA cluster; grey points: tumors; coloured points: PDOs.")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
          plot.caption = element_text(size = 9, hjust = 0))
  stem <- if (n_pages > 1) paste0("ssGSEA_Hallmark_boxplots_PDO_vs_tumors_significant_page", pg)
          else "ssGSEA_Hallmark_boxplots_PDO_vs_tumors_significant"
  save_plot(p_box, file.path(OUT_SS, stem), 12, 3.4 * ceiling(length(show) / 3) + 1)
}

# Matching table for the figures
write.csv(cmp[match(show_all, cmp$Pathway), c("Pathway", "PDO_mean", "tumor_median_all", "PDO_mean_z_vs_tumors",
                                              "PDO_percentile_in_tumors", "Wilcoxon_p", "Wilcoxon_FDR", "Higher_in")],
          file.path(OUT_SS, "ssGSEA_Hallmark_PDO_vs_tumors_significant_pathways.csv"), row.names = FALSE)

message("Done. Outputs in: ", OUT_SS)