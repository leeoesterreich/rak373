############################################################
## PURPOSE
##   Do the leading-edge genes from the estrogen-vs-chromatin correlation
##   analysis carry more WES mutations in one RNA cluster, and are mutated
##   genes expressed differently? Run separately for the chromatin and the
##   estrogen leading-edge gene sets.
##
## QUESTIONS (for each gene set)
##   1. Is an excess of leading-edge gene mutations in a cluster real, or does
##      it depend on a few low-VAF calls? (VAF >= 3%, >= 5%, >= 3% without indels)
##   2. Is it specific to this gene set? (random gene sets with the same
##      per-gene mutation counts; total WES mutation burden by cluster)
##   3. Is a mutated gene expressed differently in the tumor that carries the
##      mutation (same gene, same tumor), also after accounting for cluster?
##   4. Is the gene-set expression score different in tumors with a
##      leading-edge mutation, after accounting for cluster?
##
## INPUT
##   - leading-edge gene tables written by pILC_chromatin_estrogen_correlation.R
##     (run that script first)
##   - saved VST matrix, cluster membership CSV, WES calls (PCGR audit table)
##
## NOTE
##   All tests are exploratory (few mutated tumors). The clusters came from the
##   same expression data, and the leading-edge genes were chosen from the same
##   tumors' expression, so the expression results are descriptive.
############################################################

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ----------------------------
# Paths (fixed names, no date or time stamps in any output)
# ----------------------------
root       <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"
vst_dir    <- file.path(root, "RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47")
vst_file   <- file.path(vst_dir, "htseq_VST_matrix_all_expression_filtered_genes.rds")
annot_file <- file.path(vst_dir, "gene_annotation_from_GTF.tsv")
cluster_csv_file <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
wes_audit  <- file.path(root, paste0("WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/",
                                     "Combined_MAF/Full_MAF_nonsilent_variant_PCGR_VAF_match_audit_2026-09-23.tsv"))

corr_dir   <- file.path(root, "RNASeq/Downstream/Chromatin_Estrogen_correlation")
gene_files <- c(Chromatin = file.path(corr_dir, "tables", "chromatin_leading_edge_genes.tsv"),
                Estrogen  = file.path(corr_dir, "tables", "estrogen_leading_edge_genes.tsv"))
out_root   <- file.path(corr_dir, "mutation_vs_expression")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (f in c(vst_file, annot_file, cluster_csv_file, wes_audit, gene_files)) {
  if (!file.exists(f)) stop("Input file not found: ", f,
                            "\n(The leading-edge gene tables come from pILC_chromatin_estrogen_correlation.R.)")
}

# ----------------------------
# Settings
# ----------------------------
N_GENES <- 60          # leading-edge genes per set: top N by |rho| with the other family's score (NULL = all)
N_PERM  <- 5000        # random gene sets for the matched null

# Scenarios for the mutation calls
scenarios <- list(
  `VAF >= 3%`            = list(thr = 0.03, drop_indels = FALSE),
  `VAF >= 5%`            = list(thr = 0.05, drop_indels = FALSE),
  `VAF >= 3%, no indels` = list(thr = 0.03, drop_indels = TRUE)
)
FIG1_SCENARIOS <- c("VAF >= 3%", "VAF >= 5%")   # drawn in Fig 1

nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Splice_Region", "Translation_Start_Site"
)
indel_classes <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")
cluster_cols  <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")

# ----------------------------
# Helpers
# ----------------------------
std_id <- function(x) {                       # same sample-ID rule as the oncoplot script
  x <- toupper(trimws(as.character(x)))
  x <- sub("^X(?=TP)", "", x, perl = TRUE)
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("\\.", "-", x)
  x <- sub("^JC287-", "", x)
  x <- sub("(_T|-T)$", "", x)
  x <- sub("(?<=[0-9])T$", "", x, perl = TRUE)
  sub("-0+(?=[0-9])", "-", x, perl = TRUE)    # TP09-09 = TP09-9
}
save_plot <- function(file_base, width, height, draw_fn, res = 300) {
  pdf(paste0(file_base, ".pdf"), width = width, height = height); draw_fn(); dev.off()
  png(paste0(file_base, ".png"), width = width, height = height, units = "in", res = res); draw_fn(); dev.off()
  cat("Saved:", paste0(file_base, c(".pdf", ".png")), sep = "\n  "); cat("\n")
}
p_txt <- function(p) ifelse(is.na(p), "P = NA", ifelse(p < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p))))
theme_nat <- function() theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
        axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.caption = element_text(hjust = 0, size = 6, colour = "grey30"),
        legend.key.size = unit(3, "mm"))
safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)

# ----------------------------
# 1. Expression: saved VST matrix with gene symbols (genes x tumors)
# ----------------------------
vst <- as.matrix(readRDS(vst_file))
if (is.null(rownames(vst)) || is.null(colnames(vst))) stop("VST matrix needs gene row names and sample column names.")
storage.mode(vst) <- "double"
raw_ids <- rownames(vst)
sym <- toupper(raw_ids)
is_e <- grepl("^ENSG[0-9]+", sym)
if (any(is_e)) {                              # Ensembl IDs -> symbols through the annotation file
  annot <- fread(annot_file)
  key   <- sub("\\.[0-9]+$", "", annot$Gene)
  sym[is_e] <- toupper(annot$gene_symbol[match(sub("\\.[0-9]+$", "", raw_ids[is_e]), key)])
}
keep <- !is.na(sym) & sym != "" & !grepl("^ENSG[0-9]+", sym)
vst <- vst[keep, , drop = FALSE]; sym <- sym[keep]
o <- order(-rowMeans(vst)); vst <- vst[o, , drop = FALSE]; sym <- sym[o]
dup <- duplicated(sym)                        # shared symbol: keep the most expressed gene
mat <- vst[!dup, , drop = FALSE]; rownames(mat) <- sym[!dup]
colnames(mat) <- std_id(colnames(mat))
if (anyDuplicated(colnames(mat))) stop("Two VST columns map to the same tumor ID.")
cat("Expression matrix:", nrow(mat), "genes x", ncol(mat), "tumors\n")

# ----------------------------
# 2. RNA clusters (same file as the correlation script)
# ----------------------------
cc <- fread(cluster_csv_file)
samp_col <- intersect(c("sample", "Sample", "Name", "name"), names(cc))
samp_col <- if (length(samp_col)) samp_col[1] else names(cc)[1]
clus_col <- setdiff(grep("cluster|membership", names(cc), ignore.case = TRUE, value = TRUE), samp_col)[1]
if (is.na(clus_col)) stop("No cluster column found in ", cluster_csv_file)
message("Cluster CSV: sample column '", samp_col, "', cluster column '", clus_col, "'")
cl <- data.table(Tumor = std_id(cc[[samp_col]]),
                 Cluster = paste("Cluster", as.integer(gsub("[^0-9]", "", as.character(cc[[clus_col]])))))
if (anyNA(cl$Cluster) || any(!cl$Cluster %in% names(cluster_cols))) stop("Cluster membership must be 1, 2 or 3.")

# ----------------------------
# 3. WES calls (all genes, non-silent)
# ----------------------------
wes <- fread(wes_audit, sep = "\t", na.strings = c("", "NA", "."),
             select = c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "VAF_TUMOR_fraction"))
wes[, Tumor := std_id(Tumor_Sample_Barcode)]
wes_tumors <- unique(wes$Tumor)
calls_all <- wes[Variant_Classification %in% nonsilent_classes & !is.na(VAF_TUMOR_fraction),
                 .(Tumor, Gene = toupper(Hugo_Symbol), Class = Variant_Classification, VAF = VAF_TUMOR_fraction)]
rm(wes)

# Tumors used: WES + RNA cluster + expression
tt <- Reduce(intersect, list(wes_tumors, cl$Tumor, colnames(mat)))
tumor_cl <- cl$Cluster[match(tt, cl$Tumor)]
cl_levels <- names(cluster_cols)[names(cluster_cols) %in% tumor_cl]
cat("Tumors with WES and RNA-seq:", length(tt), "\n"); print(table(tumor_cl))
if (length(tt) < 10) stop("Fewer than 10 tumors matched across WES, clusters and expression; check the sample IDs.")
idx_cl <- lapply(cl_levels, function(k) which(tumor_cl == k)); names(idx_cl) <- cl_levels

# Statistic used for the cluster question: largest minus smallest % of mutated tumors across clusters
range_stat <- function(a) {
  p <- vapply(idx_cl, function(i) mean(a[i]), numeric(1))
  max(p) - min(p)
}

# ----------------------------
# 4. Analysis for one leading-edge gene set
# ----------------------------
run_family <- function(family, gene_file) {
  fam_dir <- file.path(out_root, tolower(family))
  dir.create(fam_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- genes
  g <- fread(gene_file)
  if (!("Gene" %in% names(g))) stop("No 'Gene' column in ", gene_file)
  g[, Gene := toupper(trimws(Gene))]
  if ("Rho_with_other_score" %in% names(g)) g <- g[order(-abs(Rho_with_other_score))]
  genes <- unique(g$Gene)
  if (!is.null(N_GENES)) genes <- head(genes, N_GENES)
  le <- intersect(genes, rownames(mat))
  cat("\n=====", family, "leading-edge genes:", length(genes), "| found in the VST matrix:", length(le), "=====\n")
  if (length(le) < 3) stop("Fewer than 3 leading-edge genes in the expression matrix for ", family)
  fwrite(data.table(Gene = genes, In_VST_matrix = genes %in% rownames(mat)),
         file.path(fam_dir, "genes_used.tsv"), sep = "\t")

  # z-scores across the tumors with WES + RNA, and the mean-z score of the gene set
  zmat <- t(scale(t(mat[le, tt, drop = FALSE])))
  zmat[!is.finite(zmat)] <- NA
  sig <- colMeans(zmat, na.rm = TRUE)

  res_tests <- list(); res_perm <- list(); res_pairs <- list(); res_pct <- list()
  null_list <- list(); tumor_tabs <- list()
  add_test <- function(sc, test, N, est, p)
    res_tests[[length(res_tests) + 1]] <<- data.table(Scenario = sc, Test = test, N = N, Estimate = est, P = p)

  for (sc in names(scenarios)) {
    s <- scenarios[[sc]]
    calls <- calls_all[VAF >= s$thr]
    if (s$drop_indels) calls <- calls[!(Class %in% indel_classes)]
    calls <- calls[Tumor %in% tt]

    # ---- tumor level: >= 1 leading-edge mutation
    le_calls <- calls[Gene %in% genes]
    n_genes <- setNames(rep(0L, length(tt)), tt)
    ng <- le_calls[, .(n = uniqueN(Gene)), by = Tumor]
    n_genes[ng$Tumor] <- ng$n
    any_mut <- n_genes > 0

    pct <- data.table(Scenario = sc, Cluster = cl_levels)
    pct[, `:=`(N = vapply(Cluster, function(k) sum(tumor_cl == k), integer(1)),
               Mutated = vapply(Cluster, function(k) sum(any_mut[tumor_cl == k]), integer(1)))]
    res_pct[[sc]] <- pct

    fp <- safe(fisher.test(table(factor(any_mut, c(FALSE, TRUE)), factor(tumor_cl, cl_levels)))$p.value)
    obs <- range_stat(any_mut)
    top_cluster <- pct$Cluster[which.max(pct$Mutated / pct$N)]

    # ---- matched random gene sets (same mutation count per gene, +/- 1 tumor)
    cnt <- calls[, .(n = uniqueN(Tumor)), by = Gene]
    universe <- union(rownames(mat), cnt$Gene)
    n_u <- setNames(rep(0L, length(universe)), universe); n_u[cnt$Gene] <- cnt$n
    mut_genes <- names(n_u)[n_u > 0 & !(names(n_u) %in% genes)]
    M <- matrix(FALSE, nrow = length(tt), ncol = length(mut_genes), dimnames = list(tt, mut_genes))
    cm <- calls[Gene %in% mut_genes]
    M[cbind(cm$Tumor, cm$Gene)] <- TRUE
    n_le <- vapply(genes, function(gg) if (gg %in% names(n_u)) as.integer(n_u[[gg]]) else 0L, integer(1))
    cand <- lapply(n_le, function(k) if (k == 0L) NULL else which(abs(n_u[mut_genes] - k) <= 1))
    set.seed(42)
    null_stat <- replicate(N_PERM, {
      pick <- unlist(lapply(cand, function(cx) if (is.null(cx) || !length(cx)) NULL else cx[sample.int(length(cx), 1)]))
      a <- if (length(pick)) rowSums(M[, pick, drop = FALSE]) > 0 else rep(FALSE, length(tt))
      range_stat(a)
    })
    perm_p <- (sum(null_stat >= obs) + 1) / (N_PERM + 1)
    null_list[[sc]] <- data.table(Scenario = sc, Null = null_stat, Observed = obs)

    # ---- total WES mutation burden by cluster
    burden <- setNames(rep(0L, length(tt)), tt)
    b <- calls[, .N, by = Tumor]; burden[b$Tumor] <- b$N
    kw_burden <- safe(kruskal.test(burden ~ factor(tumor_cl))$p.value)

    perm_row <- data.table(Scenario = sc, Cluster_with_highest_percent = top_cluster, Observed_range = obs,
                           Null_mean = mean(null_stat), Empirical_p_one_sided = perm_p,
                           Fisher_p_across_clusters = fp, Burden_KW_p = kw_burden)
    for (k in cl_levels) perm_row[[paste0("Median_burden_", sub(" ", "", k))]] <- median(burden[tumor_cl == k])
    res_perm[[sc]] <- perm_row

    # ---- expression of the mutated gene in the mutated tumor
    cell <- le_calls[Gene %in% le, .(Class = if (.N > 1) "Multi_Hit" else Class[1], VAF = max(VAF)), by = .(Tumor, Gene)]
    pairs <- cell[Tumor %in% tt]
    if (nrow(pairs) > 0) {
      pairs[, `:=`(z = zmat[cbind(Gene, Tumor)],
                   pct_rank = mapply(function(gg, t) mean(zmat[gg, ] <= zmat[gg, t], na.rm = TRUE), Gene, Tumor),
                   Cluster = tumor_cl[match(Tumor, tt)], Scenario = sc)]
      res_pairs[[sc]] <- pairs

      long <- as.data.table(as.table(zmat)); setnames(long, c("Gene", "Tumor", "z"))
      long[, Cluster := tumor_cl[match(Tumor, tt)]]
      long[, Mutated := paste(Tumor, Gene) %in% paste(pairs$Tumor, pairs$Gene)]
      if (sum(long$Mutated) >= 3 && sum(!long$Mutated) >= 3) {
        fit <- tryCatch(summary(lm(z ~ Mutated + Cluster + Gene, data = long))$coefficients, error = function(e) NULL)
        if (!is.null(fit) && "MutatedTRUE" %in% rownames(fit))
          add_test(sc, "Mutated gene expression (z) vs not, adjusted for cluster and gene (lm)",
                   sum(long$Mutated), fit["MutatedTRUE", "Estimate"], fit["MutatedTRUE", "Pr(>|t|)"])
      }
    }

    # ---- gene-set score vs leading-edge mutation
    if (sum(any_mut) > 0 && sum(!any_mut) > 0) {
      add_test(sc, "Gene-set score: tumors with >= 1 mutation vs none (Wilcoxon)", length(tt),
               median(sig[any_mut], na.rm = TRUE) - median(sig[!any_mut], na.rm = TRUE),
               safe(wilcox.test(sig[any_mut], sig[!any_mut], exact = FALSE)$p.value))
      fit2 <- tryCatch(summary(lm(sig ~ factor(tumor_cl) + any_mut))$coefficients, error = function(e) NULL)
      if (!is.null(fit2) && "any_mutTRUE" %in% rownames(fit2))
        add_test(sc, "Gene-set score ~ cluster + any mutation (lm, mutation term)", length(tt),
                 fit2["any_mutTRUE", "Estimate"], fit2["any_mutTRUE", "Pr(>|t|)"])
      sp <- suppressWarnings(tryCatch(cor.test(sig, n_genes, method = "spearman", exact = FALSE), error = function(e) NULL))
      if (!is.null(sp))
        add_test(sc, "Gene-set score vs number of mutated leading-edge genes (Spearman rho)", length(tt),
                 unname(sp$estimate), sp$p.value)
    }
    add_test(sc, "Tumors with >= 1 mutation across clusters (Fisher; estimate = largest minus smallest %)",
             length(tt), obs, fp)

    tumor_tabs[[sc]] <- data.table(Scenario = sc, Tumor = tt, Cluster = tumor_cl, n_mutated_genes = n_genes,
                                   gene_set_score = sig, total_burden = burden)
  }

  pct_dt  <- rbindlist(res_pct)[, Percent := round(100 * Mutated / N, 1)]
  perm_dt <- rbindlist(res_perm, fill = TRUE)
  tests_dt <- rbindlist(res_tests, fill = TRUE)
  pairs_dt <- rbindlist(res_pairs, fill = TRUE)
  if (nrow(tests_dt)) tests_dt[, P := signif(P, 3)]
  fwrite(pct_dt,  file.path(fam_dir, "percent_tumors_mutated_by_cluster.tsv"), sep = "\t")
  fwrite(perm_dt, file.path(fam_dir, "matched_random_geneset_test.tsv"), sep = "\t")
  fwrite(tests_dt, file.path(fam_dir, "mutation_vs_expression_tests.tsv"), sep = "\t")
  fwrite(pairs_dt, file.path(fam_dir, "mutated_gene_tumor_expression.tsv"), sep = "\t")
  fwrite(rbindlist(tumor_tabs), file.path(fam_dir, "per_tumor_summary.tsv"), sep = "\t")

  cat("\nTumors with >= 1", family, "leading-edge mutation, by cluster:\n"); print(pct_dt)
  cat("\nMatched random gene sets and mutation burden:\n"); print(perm_dt)
  cat("\nTests:\n"); print(tests_dt)
  if (nrow(pairs_dt)) {
    cat("\nExpression of each mutated gene in its mutated tumor (z across tumors, percentile rank):\n")
    print(pairs_dt[Scenario == names(scenarios)[1],
                   .(Tumor, Gene, Class, VAF = round(100 * VAF, 1), Cluster, z = round(z, 2), pct_rank = round(pct_rank, 2))][order(Gene, Tumor)])
  }

  # ---------------- Fig 1: % of tumors with >= 1 mutation, by cluster
  pct_plot <- pct_dt[Scenario %in% FIG1_SCENARIOS]
  pct_plot[, Scenario := factor(Scenario, FIG1_SCENARIOS)]
  fisher_lab <- perm_dt[Scenario %in% FIG1_SCENARIOS,
                        .(Scenario = factor(Scenario, FIG1_SCENARIOS),
                          label = paste0("Across clusters: Fisher ", p_txt(Fisher_p_across_clusters)))]
  p1 <- ggplot(pct_plot, aes(Cluster, Percent, fill = Cluster)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_text(aes(label = paste0(Mutated, "/", N)), vjust = -0.4, size = 2.6) +
    geom_text(data = fisher_lab, aes(x = 0.5, y = Inf, label = label), inherit.aes = FALSE,
              hjust = 0, vjust = 1.4, size = 2.4) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    facet_wrap(~Scenario) +
    labs(x = NULL, y = paste0("Tumors with >= 1 ", tolower(family), " leading-edge mutation (%)"),
         title = paste0(family, " leading-edge genes (", length(genes), ")"),
         caption = "WES, non-silent calls. Numbers above bars: mutated tumors / tumors in the cluster.") +
    theme_nat() + theme(plot.title = element_text(face = "bold", size = 9))
  save_plot(file.path(fam_dir, "Fig1_percent_mutated_by_cluster"), width = 5.2, height = 3.2, draw_fn = function() print(p1))

  # ---------------- Fig 2: matched random gene sets
  null_dt <- rbindlist(null_list)
  null_dt[, Scenario := factor(Scenario, names(scenarios))]
  obs_dt <- perm_dt[, .(Scenario = factor(Scenario, names(scenarios)), Observed = Observed_range,
                        label = paste0("Empirical ", p_txt(Empirical_p_one_sided)))]
  p2 <- ggplot(null_dt, aes(100 * Null)) +
    geom_histogram(bins = 30, fill = "grey75", colour = "white", linewidth = 0.2) +
    geom_vline(data = obs_dt, aes(xintercept = 100 * Observed), colour = "#D62728", linewidth = 0.6) +
    geom_text(data = obs_dt, aes(x = 100 * Observed, y = Inf, label = label), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.6, size = 2.4, colour = "#D62728") +
    facet_wrap(~Scenario) +
    labs(x = "Largest minus smallest cluster: % tumors with >= 1 mutation", y = "Random gene sets",
         title = paste0(family, " leading-edge genes versus matched random gene sets"),
         caption = paste0("Grey: ", N_PERM, " random gene sets with the same per-gene mutation counts. ",
                          "Red line: observed leading-edge value.")) +
    theme_nat() + theme(plot.title = element_text(face = "bold", size = 9))
  save_plot(file.path(fam_dir, "Fig2_matched_random_gene_sets"), width = 7.2, height = 3.2, draw_fn = function() print(p2))

  # ---------------- Fig 3: expression of each mutated gene
  sc1 <- names(scenarios)[1]
  pr <- if (nrow(pairs_dt) > 0) copy(pairs_dt[Scenario == sc1]) else data.table()
  if (nrow(pr) > 0) {
    mg <- pr[, .N, by = Gene][order(-N, Gene), Gene]
    bg <- as.data.table(as.table(zmat[mg, , drop = FALSE])); setnames(bg, c("Gene", "Tumor", "z"))
    bg[, `:=`(Gene = factor(Gene, mg), Cluster = tumor_cl[match(Tumor, tt)])]
    place <- function(d) {                    # keep mutated points and labels from overlapping within a gene
      d <- d[order(z)]
      zz <- d$z; n <- length(zz)
      dx <- numeric(n); ly <- zz
      for (i in seq_len(n)[-1]) {
        if (zz[i] - zz[i - 1] < 0.12 && dx[i - 1] == 0) dx[i] <- 0.10
        ly[i] <- max(ly[i], ly[i - 1] + 0.30)
      }
      d[, `:=`(dx = dx, ly = ly)]
      d
    }
    pr <- pr[, place(.SD), by = Gene]
    pr[, x_pt := match(Gene, mg) + dx]
    pr[, x_lab := x_pt + 0.16]
    p3 <- ggplot(bg, aes(Gene, z)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey50", fill = "grey95", linewidth = 0.3) +
      geom_point(aes(colour = Cluster), size = 0.9, alpha = 0.7,
                 position = position_jitter(width = 0.15, height = 0, seed = 42)) +
      scale_colour_manual(values = cluster_cols, guide = "none") +
      geom_segment(data = pr, aes(x = x_pt, xend = x_lab - 0.02, y = z, yend = ly), inherit.aes = FALSE,
                   colour = "grey40", linewidth = 0.25) +
      geom_point(data = pr, aes(x = x_pt, y = z, fill = Cluster), inherit.aes = FALSE, shape = 21, size = 2.8,
                 colour = "black", stroke = 0.4) +
      geom_text(data = pr, aes(x = x_lab, y = ly, label = Tumor), inherit.aes = FALSE, hjust = 0, size = 2.1) +
      scale_fill_manual(values = cluster_cols, name = "Mutated tumor (cluster)") +
      labs(x = NULL, y = "Gene expression (z-score across tumors)",
           title = paste0("Expression of each mutated ", tolower(family), " leading-edge gene (", sc1, ")"),
           caption = paste0("Small points: all tumors with WES and RNA-seq, coloured by cluster. ",
                            "Large outlined points: tumors carrying a mutation in that gene.")) +
      theme_nat() + theme(plot.title = element_text(face = "bold", size = 9),
                          axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(file.path(fam_dir, "Fig3_expression_of_mutated_genes"),
              width = max(7.2, 0.45 * length(mg) + 2.5), height = 4.2, draw_fn = function() print(p3))
  }

  # ---------------- Fig 4: gene-set score by cluster and mutation
  st <- rbindlist(tumor_tabs)[Scenario == sc1]
  lev <- c(">= 1 leading-edge mutation", "No mutation")
  st[, Mutation := factor(ifelse(n_mutated_genes > 0, lev[1], lev[2]), lev)]
  lm_p <- tests_dt[Scenario == sc1 & grepl("score ~ cluster", Test), P]
  lm_txt <- if (length(lm_p)) p_txt(lm_p[1]) else "not estimable"
  p4 <- ggplot(st, aes(Cluster, gene_set_score)) +
    geom_boxplot(aes(fill = Cluster), alpha = 0.35, outlier.shape = NA, width = 0.55, colour = "grey25", linewidth = 0.35) +
    geom_point(aes(colour = Mutation, shape = Mutation, size = Mutation),
               position = position_jitter(width = 0.12, height = 0, seed = 42)) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_colour_manual(values = setNames(c("#D62728", "#52606D"), lev), name = NULL) +
    scale_shape_manual(values = setNames(c(17, 16), lev), name = NULL) +
    scale_size_manual(values = setNames(c(2.6, 1.8), lev), name = NULL) +
    labs(x = NULL, y = paste0(family, " leading-edge score (mean z)"),
         title = paste0("Expression score by cluster and mutation (", sc1, ")"),
         caption = paste(strwrap(paste0("Mutation term in a model with cluster: ", lm_txt,
                                        " (lm). Triangles: tumors with >= 1 WES mutation in the leading-edge genes."), 75),
                         collapse = "\n")) +
    theme_nat() + theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
  save_plot(file.path(fam_dir, "Fig4_score_by_cluster_and_mutation"), width = 4.2, height = 3.8, draw_fn = function() print(p4))
  invisible(NULL)
}

# ----------------------------
# 5. Run for both gene sets
# ----------------------------
for (fam in names(gene_files)) run_family(fam, gene_files[[fam]])

writeLines(
  c("pILC leading-edge genes: WES mutations versus expression",
    "Seed: 42",
    paste0("Gene sets: leading-edge genes of the top chromatin and top estrogen pathways (top ", N_GENES,
           " genes per set by |rho| with the other family's score), from pILC_chromatin_estrogen_correlation.R."),
    "Mutations: WES non-silent calls; scenarios VAF >= 3%, VAF >= 5%, VAF >= 3% without indels.",
    "Cluster question: Fisher exact across the three clusters; statistic = largest minus smallest % mutated tumors.",
    paste0("Matched null: ", N_PERM, " random gene sets with the same per-gene mutation counts (+/- 1 tumor); one-sided empirical P."),
    "Expression: saved VST matrix, z-scored across tumors with WES and RNA-seq; score = mean z of the gene set.",
    "All tests are exploratory. Clusters and leading-edge genes were derived from the same expression data."),
  file.path(out_root, "README_mutation_vs_expression.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
cat("All outputs are in:", out_root, "\n")