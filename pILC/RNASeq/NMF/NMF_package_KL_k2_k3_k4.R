###############################################################################
# Consensus NMF with the R package "NMF", method = "KL", k = 2, 3, 4, 50 runs
#
# WHY THIS SCRIPT
#   To Run the NMF package 
#     - input (vst_transformed_data.gct), same seed (42), 50 random starts
#     - consensus matrix = fraction of runs in which two tumors share a cluster
#     - final clusters = average-linkage tree on 1 - consensus, cut at k, and the
#       clusters are numbered in the order they appear in the tree
#
# METHOD
#   "KL" = Kullback-Leibler divergence updates (Brunet et al. 2004) with a stopping
#   rule based on the stationarity of the objective value. ("brunet" instead stops
#   when the cluster assignment stops changing.)
#
# OUTPUT (folder OUT_DIR, no dates in file names; same names as the hand-written run)
#   Cophenetic_correlation.csv       cophenetic correlation and PAC for k = 2, 3, 4
#   Consensus_matrix_k{k}.csv        tumor x tumor consensus
#   Cluster_membership_k{k}.csv      cluster of each tumor ("sample","cluster")
#   NMF_best_run_k3_W.csv            gene weights W (genes x 3) of the best run (lowest objective)
#   NMF_best_run_k3_H.csv            tumor weights H (3 x tumors) of the best run
#   Run_settings_and_session_info.txt
# Each k is saved as soon as it is finished, so a time-out does not lose earlier k.
###############################################################################

suppressPackageStartupMessages(library(NMF))     # nmf(), consensus()

# ------------------------------ 1. SETTINGS ----------------------------------
ROOT     <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
GCT_FILE <- file.path(ROOT, "4_Count_file_HTSeq_stranded", "VST_data_GENCODEv47",
                      "vst_transformed_data.gct")
OUT_DIR  <- file.path(ROOT, "NMF", "R_NMFConsensus_50runs")     # <- folder you use for the NMF-package run
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(OUT_DIR))

K_RANGE <- 2:4
K_WH    <- 3            # W and H of the best run are saved for this k only
N_CLUST <- 50           # number of random starts per k
SEED    <- 42
METHOD  <- "KL"
EPS     <- 1e-9
PAC_LO  <- 0.10; PAC_HI <- 0.90     # PAC = fraction of tumor pairs with consensus in this range

# ------------------------------ 2. READ THE VST MATRIX -----------------------
g <- read.delim(GCT_FILE, skip = 2, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(g$Name))
X <- as.matrix(g[, -(1:2)]); rownames(X) <- g$Name
storage.mode(X) <- "double"
cat("Matrix:", nrow(X), "genes x", ncol(X), "tumors | method:", METHOD,
    "| NMF package", as.character(packageVersion("NMF")), "\n")
stopifnot(ncol(X) == 31, !anyNA(X), min(X) >= 0)    # NMF needs non-negative values
X <- pmax(X, EPS)                                   # same small floor as the hand-written script

# ------------------------------ 3. NMF FOR k = 2, 3, 4 -----------------------
summary_tab <- data.frame()
for (k in K_RANGE) {
  cat("\n=== k =", k, "===\n")
  t0 <- Sys.time()
  fit <- nmf(X, rank = k, method = METHOD, nrun = N_CLUST, seed = SEED)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat("  finished in", mins, "minutes\n")

  # consensus matrix from the package (fraction of runs in which two tumors cluster together)
  C <- as.matrix(consensus(fit))
  dimnames(C) <- list(colnames(X), colnames(X))

  # final clusters, done as in the hand-written script
  d  <- as.dist(1 - C)
  hc <- hclust(d, method = "average")
  coph <- cor(as.vector(d), as.vector(cophenetic(hc)))
  cl <- cutree(hc, k)
  cl_ord <- match(cl, unique(cl[hc$order]))         # 1, 2, 3 ... in tree order
  names(cl_ord) <- colnames(X)

  # ---- best run, k = K_WH only: W (genes x k) and H (k x tumors), for projecting new samples (e.g. PDOs) ----
  if (k == K_WH) {
    best <- NMF::fit(fit)                           # run with the lowest objective among the N_CLUST runs
    W <- basis(best); H <- coef(best)
    rownames(W) <- rownames(X); colnames(W) <- paste0("P", seq_len(k))
    rownames(H) <- paste0("P", seq_len(k)); colnames(H) <- colnames(X)
    write.csv(W, file.path(OUT_DIR, sprintf("NMF_best_run_k%d_W.csv", k)))
    write.csv(H, file.path(OUT_DIR, sprintf("NMF_best_run_k%d_H.csv", k)))
  }

  ut  <- upper.tri(C)
  pac <- mean(C[ut] > PAC_LO & C[ut] < PAC_HI)

  write.csv(C, file.path(OUT_DIR, sprintf("Consensus_matrix_k%d.csv", k)))
  write.csv(data.frame(sample = names(cl_ord), cluster = cl_ord),
            file.path(OUT_DIR, sprintf("Cluster_membership_k%d.csv", k)), row.names = FALSE)

  summary_tab <- rbind(summary_tab,
                       data.frame(k = k, cophenetic_correlation = round(coph, 4),
                                  PAC = round(pac, 4), minutes = mins,
                                  cluster_sizes = paste(table(cl_ord), collapse = "/")))
  write.csv(summary_tab, file.path(OUT_DIR, "Cophenetic_correlation.csv"), row.names = FALSE)
  print(summary_tab)
}

# ------------------------------ 4. SAVE SETTINGS -----------------------------
sink(file.path(OUT_DIR, "Run_settings_and_session_info.txt"))
cat("method:", METHOD, " genes:", nrow(X), " tumors:", ncol(X), "\n")
cat("k:", K_RANGE, " runs per k:", N_CLUST, " seed:", SEED, "\n\n")
print(sessionInfo())
sink()
cat("\nDone. Results in:", OUT_DIR, "\n")