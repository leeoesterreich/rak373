# ==============================================================================
# pILC RNA-seq, NMF k = 3 clusters: are they real, stable and reproducible?
#
#   Part 1  Is k = 3 the best solution?  k = 2-4: cophenetic correlation,
#           consensus CDF, PAC, silhouette widths (k = 3). The saved consensus
#           matrices (Consensus_matrix_k2/3/4.csv) are reused; nothing is refitted.
#   Part 2  Is the separation stronger than chance?  SigClust, each cluster pair
#   Part 3  Are the clusters robust?  80% sample subsets, gene subsets,
#           other clustering methods (hierarchical, k-means)
#   Part 4  Can a classifier reproduce them?  Nearest centroid, leave-one-out CV
#           (genes selected INSIDE each fold) + label-permutation test
#
# Cluster labels of record = the GenePattern k = 3 membership (numbers 1, 2, 3).
# Matrix = the stranded GENCODE v47 VST matrix used by the R_NMFConsensus script, and
# NMF = that script's own KL-divergence implementation (nmf_kl), so the validation tests
# the same method that produced the saved consensus matrices. No NMF package is needed.
# Seed 42 everywhere. No dates in any file name.
#
# Outputs (in out_dir):
#   Table1_k_selection.csv          Fig1_k_selection.pdf/.png
#   Table2_sigclust_pairs.csv
#   Table3_stability_summary.csv    Fig3_stability.pdf/.png
#   Table4_loocv_summary.csv        Fig4_loocv_permutation.pdf/.png
#   Table5_sample_level.csv         Fig2_silhouette_k3.pdf/.png
#
# Set QUICK <- TRUE first to test the script in a few minutes with small settings.
# ==============================================================================

suppressPackageStartupMessages({
  library(cluster)    # silhouette()
  library(mclust)     # adjustedRandIndex()
  library(sigclust)   # sigclust()
  library(parallel)   # mclapply()
})

SEED <- 42
set.seed(SEED)

# ------------------------------------------------------------------------------
# 0. Settings (edit here only)
# ------------------------------------------------------------------------------
# QUICK = small settings, only to test that the script runs. Switch it on without editing
# the file:  QUICK=1 Rscript Cluster_validation_k3.R   (results then go to a separate folder)
options(warn = 1)   # print warnings immediately
QUICK <- identical(Sys.getenv("QUICK"), "1")

root     <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
vst_file <- file.path(root, "4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/vst_transformed_data.gct")  # NMF input (GCT)
gp_file  <- file.path(root, "Archive/NMF/GenePattern_NMFConsensus/vst_transformed_data.consensus.k.3.gct")
out_dir  <- file.path(root, "NMF/Cluster_validation")
nmf_dir  <- file.path(root, "NMF/R_NMFConsensus_50runs")      # saved Consensus_matrix_k2/3/4.csv

K          <- 3                    # number of clusters of record
KS         <- 2:4                  # Part 1: k values compared (n = 31 limits useful k)
N_CLUST    <- 50                   # Part 1: NMF random starts per k if a k must be recomputed
MAX_ITER   <- 2000                 # NMF: max iterations (as in the R_NMFConsensus script)
STOP_FREQ  <- 10                   # NMF: convergence check every 10 iterations
STOP_CONV  <- 40                   # NMF: stop after 40 unchanged checks
EPS        <- 1e-9                 # NMF: avoids division by zero / log(0)
PAC_LOW    <- 0.10                 # Part 1: PAC window (ambiguous consensus values)
PAC_HIGH   <- 0.90
NSIM_SC    <- 1000                 # Part 2: SigClust simulations
SC_GENES   <- 5000                 # Part 2: SigClust uses the most variable genes (Inf = all; slow, memory-heavy)
FRAC       <- 0.80                 # Part 3: fraction of samples / genes kept per repeat
B_SAMPLE   <- 100                  # Part 3: number of 80% sample subsets
B_GENE     <- 100                  # Part 3: number of 80% gene subsets
NRUN_STAB  <- 10                   # Part 3: NMF random starts per perturbed run
TOP_GENES  <- c(1000, 2000, 5000, 10000)   # Part 3: most-variable-gene sets
STAB_THR   <- 0.80                 # Part 3: samples below this are called borderline (exploratory)
N_SELECT   <- 200                  # Part 4: genes chosen per fold (primary analysis)
N_SELECT_SENS <- c(50, 100, 200, 500)      # Part 4: sensitivity (observed accuracy only)
NPERM      <- 1000                 # Part 4: label permutations

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
dc <- detectCores(); if (is.na(dc)) dc <- 1L
N_CORES <- max(1L, min(N_CORES, dc))

if (QUICK) {
  N_CLUST <- 5; NSIM_SC <- 100; B_SAMPLE <- 6; B_GENE <- 6; NRUN_STAB <- 3
  TOP_GENES <- c(1000, 2000); NPERM <- 50
  out_dir <- paste0(out_dir, "_quick_test")      # test results never overwrite the real ones
}

# Create the output folder now and stop at once if that is not possible
# (a failure here used to stay hidden until the first file was written).
dir.create(out_dir, showWarnings = TRUE, recursive = TRUE)
if (!dir.exists(out_dir) || file.access(out_dir, 2) != 0)
  stop("Cannot create or write to the output folder: ", out_dir,
       "\nCheck the path and your write permission, or change out_dir above.")

# ------------------------------------------------------------------------------
# 1. Data and cluster labels of record
# ------------------------------------------------------------------------------
# "TP13.1150_S40" -> "TP13-1150" (remove _S## suffix, dots -> hyphens)
clean_id <- function(x) gsub("\\.", "-", sub("_S[0-9]+$", "", x))

# GCT: 2 header lines, then Name, Description, one column per sample (as in the NMF script)
g <- read.delim(vst_file, skip = 2, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(g$Name))
vst <- as.matrix(g[, -(1:2)]); rownames(vst) <- g$Name
storage.mode(vst) <- "double"
colnames(vst) <- clean_id(colnames(vst))
stopifnot(ncol(vst) == 31, !anyNA(vst), min(vst) >= 0)    # NMF needs complete, non-negative data
vst <- pmax(vst, EPS)

# GenePattern membership list: Name, Description, membership.ordered
gp <- read.table(gp_file, skip = 2, header = TRUE, sep = "\t",
                 check.names = FALSE, quote = "", comment.char = "")
gp_lab <- setNames(as.integer(gp$membership.ordered), clean_id(gp$Name))
stopifnot(setequal(names(gp_lab), colnames(vst)))
gp_lab <- gp_lab[colnames(vst)]                    # same sample order as the matrix

samples <- colnames(vst)
n <- ncol(vst); p <- nrow(vst)
y <- unname(gp_lab)
cat("Matrix:", p, "genes x", n, "samples | cluster sizes:", table(y), "\n")

gene_rank <- order(apply(vst, 1, var), decreasing = TRUE)   # most variable first (no labels used)

# Colours (Okabe-Ito, colour-blind safe)
ccol <- c("#0072B2", "#E69F00", "#009E73")        # clusters 1, 2, 3
kcol <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")

# ---- helpers -----------------------------------------------------------------
ari <- function(a, b) adjustedRandIndex(a, b)

all_perms <- function(v) {                          # all orderings of v (k is small)
  if (length(v) <= 1) return(list(v))
  out <- list()
  for (i in seq_along(v)) for (pp in all_perms(v[-i])) out[[length(out) + 1]] <- c(v[i], pp)
  out
}
# Recode the labels in `lab` so they agree with `ref` as much as possible.
align_labels <- function(lab, ref, k) {
  best <- lab; best_agree <- -1
  for (pp in all_perms(seq_len(k))) {
    new <- pp[lab]
    a <- sum(new == ref)
    if (a > best_agree) { best_agree <- a; best <- new }
  }
  best
}

# One NMF run, KL-divergence multiplicative updates (Lee & Seung). Copied unchanged
# from the R_NMFConsensus script. X = genes x samples. Returns the cluster of each sample.
nmf_kl <- function(X, k) {
  n <- nrow(X); m <- ncol(X)
  W <- matrix(runif(n * k), n, k) + EPS
  H <- matrix(runif(k * m), k, m) + EPS
  prev <- NULL; stable <- 0L; it <- 0L
  repeat {
    it <- it + 1L
    WH <- W %*% H + EPS
    H  <- H * sweep(crossprod(W, X / WH), 1, colSums(W), "/")
    WH <- W %*% H + EPS
    W  <- W * sweep(tcrossprod(X / WH, H), 2, rowSums(H), "/")
    W[W < EPS] <- EPS; H[H < EPS] <- EPS
    if (it %% STOP_FREQ == 0L) {
      cl <- max.col(t(H), ties.method = "first")
      if (!is.null(prev) && identical(cl, prev)) stable <- stable + 1L else stable <- 0L
      prev <- cl
      if (stable >= STOP_CONV) break
    }
    if (it >= MAX_ITER) break
  }
  max.col(t(H), ties.method = "first")
}

# Consensus matrix over `nrun` random starts (same logic and seeding as consensus_k in
# the R_NMFConsensus script: set.seed(seed), then nrun runs in a row).
consensus_matrix <- function(X, k, nrun, seed) {
  set.seed(seed)
  m <- ncol(X)
  C <- matrix(0, m, m)
  for (r in seq_len(nrun)) {
    cl <- nmf_kl(X, k)
    C  <- C + outer(cl, cl, "==")
  }
  C <- C / nrun
  dimnames(C) <- list(colnames(X), colnames(X))
  C
}

# Save one figure as PDF and PNG
save_fig <- function(name, w, h, plot_fun) {
  pdf(file.path(out_dir, paste0(name, ".pdf")), width = w, height = h); plot_fun(); dev.off()
  # cairo works without a display (batch jobs); plain png() can fail on compute nodes
  png(file.path(out_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300,
      type = if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib")
  plot_fun(); dev.off()
}

# ==============================================================================
# PART 1. Is k = 3 the best solution?
# ==============================================================================
cat("\n== Part 1: k selection (k =", min(KS), "to", max(KS), ") ==\n")
stopifnot(K %in% KS)

# Fallback if a saved matrix is missing or unreadable: recompute it with the same method,
# starts and seed as the R_NMFConsensus script.
fit_consensus <- function(k) consensus_matrix(vst, k, N_CLUST, SEED)

# Reuse the consensus matrices already saved by the 50-run R-NMF (Consensus_matrix_k#.csv).
# A file is accepted only if it is a square matrix with the same 31 samples and values in [0, 1];
# otherwise that k is recomputed, so a bad file costs time, not correctness.
read_consensus <- function(k) {
  f <- file.path(nmf_dir, paste0("Consensus_matrix_k", k, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.numeric(d[[1]])) { rn <- d[[1]]; d <- d[, -1, drop = FALSE] } else rn <- colnames(d)
  C <- as.matrix(d)
  rownames(C) <- clean_id(rn); colnames(C) <- clean_id(colnames(d))
  ok <- nrow(C) == n && ncol(C) == n && setequal(rownames(C), samples) &&
        setequal(colnames(C), samples) && all(C >= 0 & C <= 1)
  if (!ok) { warning("Saved consensus matrix for k = ", k, " has an unexpected format; recomputing"); return(NULL) }
  C[samples, samples]
}
cons <- lapply(KS, read_consensus)
names(cons) <- paste0("k", KS)
todo <- KS[vapply(cons, is.null, logical(1))]
cat("Loaded from saved files: k =", setdiff(KS, todo), "| computed here: k =", todo, "\n")
if (length(todo) > 0) {
  new <- mclapply(todo, fit_consensus, mc.cores = min(length(todo), N_CORES), mc.preschedule = FALSE)
  stopifnot(!any(vapply(new, function(x) inherits(x, "try-error") || is.null(x), logical(1))))
  cons[paste0("k", todo)] <- new
}

upper_vals <- function(M) M[upper.tri(M)]

k_tbl <- do.call(rbind, lapply(seq_along(KS), function(i) {
  k <- KS[i]; C <- cons[[i]]
  d  <- as.dist(1 - C)
  hc <- hclust(d, method = "average")
  Fn <- ecdf(upper_vals(C))
  cl <- cutree(hc, k)
  data.frame(k = k,
             source           = if (k %in% todo) "computed here" else "saved file",
             cophenetic       = cor(d, cophenetic(hc)),               # higher = more stable
             PAC              = Fn(PAC_HIGH) - Fn(PAC_LOW),           # lower = fewer ambiguous pairs
             mean_silhouette  = mean(silhouette(cl, d)[, "sil_width"]),
             min_cluster_size = min(table(cl)))
}))
print(k_tbl, row.names = FALSE)
write.csv(k_tbl, file.path(out_dir, "Table1_k_selection.csv"), row.names = FALSE)

# Does R-NMF at k = 3 reproduce the GenePattern clusters? (information only)
C3  <- cons[[paste0("k", K)]]
cl3 <- cutree(hclust(as.dist(1 - C3), "average"), K)
cat("ARI, R-NMF consensus clusters vs GenePattern (k =", K, "):", round(ari(cl3, y), 3), "\n")

# Silhouette widths for the clusters of record (k = 3), two distances:
#   consensus (1 - consensus) and expression (1 - Pearson correlation, all genes)
sil_c <- silhouette(y, as.dist(1 - C3))
sil_e <- silhouette(y, as.dist(1 - cor(vst)))
sil_tbl <- data.frame(sample = samples, cluster = y,
                      nearest_other_cluster = sil_c[, "neighbor"],
                      silhouette_consensus  = sil_c[, "sil_width"],
                      silhouette_expression = sil_e[, "sil_width"])
cat("Mean silhouette at k =", K, "| consensus:", round(mean(sil_tbl$silhouette_consensus), 3),
    "| expression:", round(mean(sil_tbl$silhouette_expression), 3), "\n")
cat("Samples with negative consensus silhouette:",
    paste(sil_tbl$sample[sil_tbl$silhouette_consensus < 0], collapse = ", "), "\n")

save_fig("Fig1_k_selection", 11, 4, function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Consensus value", ylab = "CDF",
       main = "Consensus CDF")
  for (i in seq_along(KS)) {
    v <- sort(upper_vals(cons[[i]]))
    lines(v, seq_along(v) / length(v), col = kcol[i], lwd = 2)
  }
  abline(v = c(PAC_LOW, PAC_HIGH), lty = 3, col = "grey50")
  legend("bottomright", legend = paste("k =", KS), col = kcol[seq_along(KS)], lwd = 2, bty = "n")
  plot(k_tbl$k, k_tbl$cophenetic, type = "b", pch = 19, xaxt = "n", xlab = "k",
       ylab = "Cophenetic correlation", main = "Cophenetic (higher = better)")
  axis(1, KS)
  plot(k_tbl$k, k_tbl$PAC, type = "b", pch = 19, xaxt = "n", xlab = "k",
       ylab = sprintf("PAC (%.2f-%.2f)", PAC_LOW, PAC_HIGH), main = "PAC (lower = better)")
  axis(1, KS)
})

save_fig("Fig2_silhouette_k3", 10, 6, function() {
  par(mfrow = c(1, 2), mar = c(4.5, 6, 3, 1))
  for (nm in c("silhouette_consensus", "silhouette_expression")) {
    o <- order(sil_tbl$cluster, -sil_tbl[[nm]])
    barplot(rev(sil_tbl[[nm]][o]), horiz = TRUE, names.arg = rev(sil_tbl$sample[o]),
            las = 1, cex.names = 0.6, col = ccol[rev(sil_tbl$cluster[o])], border = NA,
            xlim = c(-0.5, 1), xlab = "Silhouette width",
            main = sprintf("%s (mean %.2f)", sub("silhouette_", "", nm), mean(sil_tbl[[nm]])))
    abline(v = 0)
  }
})

# ==============================================================================
# PART 2. Is the separation stronger than chance?  SigClust (Liu et al. 2008)
# ==============================================================================
# For each pair of clusters: null = ONE Gaussian cloud with the same covariance
# eigenvalues. The statistic is the 2-means cluster index (within-SS / total-SS);
# a small index = strong separation. Our NMF labels are used as the split.
cat("\n== Part 2: SigClust ==\n")
cl_pairs <- combn(sort(unique(y)), 2, simplify = FALSE)
# SigClust on ~20,000+ genes is very slow and needs a lot of memory. It is normally run on the
# most variable genes. The ranking uses all samples and no labels (same ranking as Part 3).
sc_genes <- gene_rank[seq_len(min(SC_GENES, p))]
cat("SigClust uses the", length(sc_genes), "most variable genes\n")

# Cluster index computed by hand: within-cluster SS / total SS
cluster_index <- function(x, lab) {
  ctr <- function(m) sum(sweep(m, 2, colMeans(m))^2)
  sum(sapply(unique(lab), function(l) ctr(x[lab == l, , drop = FALSE]))) / ctr(x)
}

sc_tbl <- do.call(rbind, lapply(cl_pairs, function(pr) {
  cat("  SigClust: cluster", pr[1], "vs", pr[2], "\n")
  keep <- samples[y %in% pr]
  x    <- t(vst[sc_genes, keep, drop = FALSE])              # samples x genes
  lab  <- ifelse(y[samples %in% keep] == pr[1], 1, 2)
  ci_hand <- cluster_index(x, lab)

  run_sc <- function(use_labels) {
    set.seed(SEED)
    tryCatch(
      if (use_labels) sigclust(x, nsim = NSIM_SC, nrep = 1, labflag = 1, label = lab, icovest = 1)
      else            sigclust(x, nsim = NSIM_SC, nrep = 1, labflag = 0, icovest = 1),
      error = function(e) { message("SigClust failed for pair ", pr[1], "-", pr[2], ": ",
                                    conditionMessage(e)); NULL })
  }
  get <- function(res, s) tryCatch(as.numeric(slot(res, s)), error = function(e) NA_real_)

  res <- run_sc(TRUE); split_used <- "NMF labels"
  # Safety check: the index SigClust reports must equal the index of OUR split.
  # If not, the label coding was not understood: use SigClust's own 2-means split instead.
  if (is.null(res) || !isTRUE(abs(get(res, "xcindex") - ci_hand) < 1e-6)) {
    warning("SigClust did not use the NMF labels for pair ", pr[1], "-", pr[2],
            "; using its own 2-means split instead")
    res <- run_sc(FALSE); split_used <- "2-means"
  }
  data.frame(cluster_a = pr[1], cluster_b = pr[2], n_a = sum(lab == 1), n_b = sum(lab == 2),
             n_genes = length(sc_genes), split_used = split_used,
             cluster_index_NMF_split = ci_hand,
             cluster_index_sigclust  = get(res, "xcindex"),
             p_gaussian_fit = get(res, "pvalnorm"),   # smooth p-value (can be below 1/nsim)
             p_empirical    = get(res, "pval"))          # limited to >= 1/nsim
}))
sc_tbl$p_holm <- p.adjust(sc_tbl$p_gaussian_fit, method = "holm")
sc_tbl$p_BH   <- p.adjust(sc_tbl$p_gaussian_fit, method = "BH")
print(sc_tbl, row.names = FALSE)
write.csv(sc_tbl, file.path(out_dir, "Table2_sigclust_pairs.csv"), row.names = FALSE)

# ==============================================================================
# PART 3. Are the clusters robust?
# ==============================================================================
cat("\n== Part 3: stability ==\n")
set.seed(SEED)
sub_idx  <- replicate(B_SAMPLE, sort(sample.int(n, round(FRAC * n))), simplify = FALSE)
gene_idx <- replicate(B_GENE,   sort(sample.int(p, round(FRAC * p))), simplify = FALSE)

jobs <- c(
  lapply(seq_len(B_SAMPLE), function(b)
    list(type = "samples_80pct", samples = sub_idx[[b]], genes = seq_len(p), seed = SEED + b)),
  lapply(seq_len(B_GENE), function(b)
    list(type = "genes_80pct", samples = seq_len(n), genes = gene_idx[[b]], seed = SEED + 1000 + b)),
  lapply(seq_along(TOP_GENES), function(i)
    list(type = paste0("top", TOP_GENES[i], "_variable_genes"), samples = seq_len(n),
         genes = gene_rank[seq_len(min(TOP_GENES[i], p))], seed = SEED + 2000 + i))
)
cat("Running", length(jobs), "perturbed NMF fits (", NRUN_STAB, "starts each) on",
    N_CORES, "cores\n")

run_job <- function(j) {                              # consensus clusters of one perturbed run
  C <- consensus_matrix(vst[j$genes, j$samples, drop = FALSE], K, NRUN_STAB, j$seed)
  cutree(hclust(as.dist(1 - C), "average"), K)
}
res_list <- mclapply(jobs, function(j) tryCatch(run_job(j), error = function(e) NULL),
                     mc.cores = N_CORES, mc.preschedule = FALSE)
failed <- vapply(res_list, is.null, logical(1))
if (any(failed)) warning(sum(failed), " perturbed runs failed and were skipped")

run_tbl <- kept_tbl <- NULL
for (i in which(!failed)) {
  sm   <- samples[jobs[[i]]$samples]
  ref  <- gp_lab[sm]
  lab  <- res_list[[i]]
  kept <- align_labels(lab, ref, K) == ref              # did the sample keep its cluster?
  run_tbl  <- rbind(run_tbl,  data.frame(run = i, type = jobs[[i]]$type, ARI = ari(lab, ref)))
  kept_tbl <- rbind(kept_tbl, data.frame(run = i, type = jobs[[i]]$type, sample = sm, kept = kept))
}
nmf_types <- unique(run_tbl$type)

# Other methods on the full cohort (all genes and the 5,000 most variable genes)
alt <- list()
for (gs in c("all_genes", "top5000")) {
  genes_use <- if (gs == "all_genes") seq_len(p) else gene_rank[seq_len(min(5000, p))]
  X <- t(vst[genes_use, , drop = FALSE])                # samples x genes
  alt[[paste0("hclust_ward_euclid | ", gs)]]  <- cutree(hclust(dist(X), "ward.D2"), K)
  alt[[paste0("hclust_average_corr | ", gs)]] <- cutree(hclust(as.dist(1 - cor(t(X))), "average"), K)
  set.seed(SEED)
  alt[[paste0("kmeans | ", gs)]] <- kmeans(X, centers = K, nstart = 50, iter.max = 100)$cluster
}
alt_kept <- sapply(alt, function(l) align_labels(l, y, K) == y)     # samples x methods

summ_nmf <- do.call(rbind, lapply(nmf_types, function(tp) {
  a <- run_tbl$ARI[run_tbl$type == tp]
  data.frame(analysis = tp, n_runs = length(a), mean_ARI = mean(a), median_ARI = median(a),
             p05_ARI = unname(quantile(a, 0.05)), min_ARI = min(a), frac_identical = mean(a > 0.999))
}))
summ_alt <- do.call(rbind, lapply(names(alt), function(m) {
  a <- ari(alt[[m]], y)
  data.frame(analysis = paste0("method: ", m), n_runs = 1, mean_ARI = a, median_ARI = a,
             p05_ARI = a, min_ARI = a, frac_identical = as.numeric(a > 0.999))
}))
stab_summary <- rbind(summ_nmf, summ_alt)
print(stab_summary, row.names = FALSE, digits = 3)
write.csv(stab_summary, file.path(out_dir, "Table3_stability_summary.csv"), row.names = FALSE)

# Fraction of runs in which each sample keeps its cluster
per_sample <- function(types) {
  d <- kept_tbl[kept_tbl$type %in% types, ]
  as.numeric(tapply(d$kept, factor(d$sample, levels = samples), mean))
}
stab_samples <- per_sample("samples_80pct")
stab_genes   <- per_sample(setdiff(nmf_types, "samples_80pct"))
stab_overall <- per_sample(nmf_types)
stab_other   <- rowMeans(alt_kept)

save_fig("Fig3_stability", 12, 6, function() {
  par(mfrow = c(1, 2), mar = c(11, 4.5, 3, 1))
  boxplot(ARI ~ type, data = run_tbl, las = 2, ylim = c(0, 1), col = "grey90",
          ylab = "ARI vs GenePattern clusters", xlab = "", main = "Cluster stability (k = 3)")
  o <- order(y, stab_overall)
  barplot(stab_overall[o], names.arg = samples[o], las = 2, cex.names = 0.7, col = ccol[y[o]],
          border = NA, ylim = c(0, 1), ylab = "Fraction of runs keeping the cluster",
          main = "Per-sample stability")
  abline(h = STAB_THR, lty = 2)
})

# ==============================================================================
# PART 4. Can a nearest-centroid classifier reproduce the clusters?
# ==============================================================================
cat("\n== Part 4: nearest-centroid LOOCV + label permutation ==\n")
X <- t(vst)                                              # samples x genes

# One-way ANOVA F statistic for every gene (vectorised)
f_stat <- function(Xm, lab) {
  kk <- K; nn <- nrow(Xm)
  ng <- tabulate(lab, kk)
  M  <- rowsum(Xm, lab) / ng                             # kk x genes: class means
  grand <- colMeans(Xm)
  bss <- colSums(ng * sweep(M, 2, grand)^2)
  wss <- colSums((Xm - M[lab, , drop = FALSE])^2)
  (bss / (kk - 1)) / (wss / (nn - kk))
}

# Leave-one-out CV. EVERYTHING that uses labels (gene selection, centroids,
# scaling) is computed on the 30 training samples only.
loocv_nc <- function(Xm, lab, G) {
  nn <- nrow(Xm); pred <- integer(nn)
  for (i in seq_len(nn)) {
    tr  <- setdiff(seq_len(nn), i)
    Xtr <- Xm[tr, , drop = FALSE]; ytr <- lab[tr]
    sel <- order(f_stat(Xtr, ytr), decreasing = TRUE)[seq_len(G)]   # top-G genes, training only
    Xs  <- Xtr[, sel, drop = FALSE]
    mu  <- colMeans(Xs)
    sdv <- sqrt(colSums(sweep(Xs, 2, mu)^2) / (length(tr) - 1)); sdv[sdv == 0] <- 1
    Z    <- sweep(sweep(Xs, 2, mu), 2, sdv, "/")
    cent <- rowsum(Z, ytr) / tabulate(ytr, K)                        # K x G centroids
    z    <- (Xm[i, sel] - mu) / sdv
    pred[i] <- which.min(rowSums(sweep(cent, 2, z)^2))               # nearest centroid
  }
  pred
}

pred_obs <- loocv_nc(X, y, N_SELECT)
acc_obs  <- mean(pred_obs == y)
conf     <- table(true = y, predicted = pred_obs)
recall   <- tapply(pred_obs == y, y, mean)
cat("LOOCV accuracy (", N_SELECT, "genes per fold):", round(acc_obs, 3), "\n")
print(conf)

sens <- data.frame(genes_per_fold = N_SELECT_SENS,
                   accuracy = sapply(N_SELECT_SENS, function(g) mean(loocv_nc(X, y, g) == y)))
print(sens, row.names = FALSE)

# Permutation test: shuffle the labels, repeat the WHOLE procedure (incl. gene selection)
set.seed(SEED)
perm_labs <- replicate(NPERM, sample(y), simplify = FALSE)
acc_perm  <- unlist(mclapply(perm_labs, function(yp) mean(loocv_nc(X, yp, N_SELECT) == yp),
                             mc.cores = N_CORES, mc.preschedule = TRUE))
p_perm <- (1 + sum(acc_perm >= acc_obs)) / (NPERM + 1)              # empirical p-value
cat("Permutation p-value:", signif(p_perm, 3), "(smallest possible =",
    signif(1 / (NPERM + 1), 3), ")\n")

loocv_summary <- data.frame(
  genes_per_fold = N_SELECT, accuracy = acc_obs,
  recall_cluster1 = recall[["1"]], recall_cluster2 = recall[["2"]], recall_cluster3 = recall[["3"]],
  null_mean = mean(acc_perm), null_p95 = unname(quantile(acc_perm, 0.95)),
  n_permutations = NPERM, p_permutation = p_perm)
write.csv(loocv_summary, file.path(out_dir, "Table4_loocv_summary.csv"), row.names = FALSE)

save_fig("Fig4_loocv_permutation", 6, 4.5, function() {
  par(mar = c(4.5, 4.5, 3, 1))
  hist(acc_perm, breaks = 20, xlim = c(0, 1), col = "grey85", border = "white",
       xlab = "LOOCV accuracy", main = "Shuffled-label null")
  abline(v = acc_obs, col = "firebrick", lwd = 2)
  legend("topright", bty = "n", legend = c(sprintf("observed = %.2f", acc_obs),
                                           sprintf("p = %.3g", p_perm)))
})

# ==============================================================================
# One table with everything per sample
# ==============================================================================
sample_tbl <- data.frame(
  sample = samples, cluster = y,
  nearest_other_cluster    = sil_tbl$nearest_other_cluster,
  silhouette_consensus     = sil_tbl$silhouette_consensus,
  silhouette_expression    = sil_tbl$silhouette_expression,
  stability_sample_subsets = stab_samples,
  stability_gene_subsets   = stab_genes,
  stability_overall        = stab_overall,
  stability_other_methods  = stab_other,
  loocv_predicted          = pred_obs,
  loocv_correct            = pred_obs == y,
  low_stability            = stab_overall < STAB_THR,
  negative_silhouette      = sil_tbl$silhouette_consensus < 0)
write.csv(sample_tbl, file.path(out_dir, "Table5_sample_level.csv"), row.names = FALSE)

cat("\nBorderline samples (stability <", STAB_THR, "or negative silhouette or LOOCV error):\n")
print(sample_tbl[sample_tbl$low_stability | sample_tbl$negative_silhouette | !sample_tbl$loocv_correct,
                 c("sample", "cluster", "silhouette_consensus", "stability_overall", "loocv_predicted")],
      row.names = FALSE)
cat("\nDone. Outputs in:", out_dir, "\n")
print(sessionInfo())