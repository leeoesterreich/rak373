###############################################################################
# Plots and GenePattern-numbered clusters for the NMF-package KL run
#
# WHY THIS SCRIPT
#   The NMF run itself (NMF_package_KL_k2_k3_k4.R) is slow and is already finished.
#   This script only READS its saved files, so it runs in seconds and needs no
#   NMF package. It adds:
#     1. Cophenetic_coefficient_plot   cophenetic correlation for k = 2, 3, 4
#     2. Consensus_heatmap_k2/3/4      consensus matrix with tree and cluster bar
#     3. Cluster_membership_k3_GPlabels.csv
#        The package numbers its clusters arbitrarily, so the k = 3 clusters are
#        renumbered to the GenePattern numbers: each package cluster gets the
#        GenePattern cluster that shares the MOST samples with it (one-to-one).
#        Output columns: "sample","cluster".
#
# INPUT  (folder IN_DIR, written by NMF_package_KL_k2_k3_k4.R)
#   Cophenetic_correlation.csv, Consensus_matrix_k{k}.csv, Cluster_membership_k{k}.csv
#   GenePattern consensus.k.3.gct (tab-separated membership list)
# OUTPUT (same folder, PNG + PDF for the figures, no dates in file names)
###############################################################################

library(mclust)      # adjustedRandIndex()

# ------------------------------ 1. SETTINGS ----------------------------------
ROOT   <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
IN_DIR <- file.path(ROOT, "NMF", "NMFpackage_KL_50runs")
GP_FILE <- file.path(ROOT, "Archive/NMF/GenePattern_NMFConsensus",
                     "vst_transformed_data.consensus.k.3.gct")
K_RANGE <- 2:4
K_GP    <- 3                                   # the k that is compared with GenePattern
PAL     <- c("#4E79A7", "#F28E2B", "#59A14F", "#B07AA1")   # cluster colours (as in the hand-written run)

# "TP13.1150_S40" -> "TP13-1150": makes sample names comparable between files
clean_id <- function(x) gsub("\\.", "-", sub("_S[0-9]+$", "", x))

# Save one figure as PNG and PDF
save_fig <- function(stub, w, h, draw) {
  png(paste0(stub, ".png"), width = w * 250, height = h * 250, res = 250); draw(); dev.off()
  pdf(paste0(stub, ".pdf"), width = w, height = h);                         draw(); dev.off()
}

# ------------------------------ 2. READ THE SAVED RESULTS --------------------
coph_tab <- read.csv(file.path(IN_DIR, "Cophenetic_correlation.csv"))
print(coph_tab)

res <- list()
for (k in K_RANGE) {
  C <- as.matrix(read.csv(file.path(IN_DIR, sprintf("Consensus_matrix_k%d.csv", k)),
                          row.names = 1, check.names = FALSE))
  stopifnot(nrow(C) == ncol(C), all(rownames(C) == colnames(C)))
  # tree and clusters rebuilt exactly as in the run script (average linkage on 1 - consensus)
  hc <- hclust(as.dist(1 - C), method = "average")
  cl <- cutree(hc, k)
  cl_ord <- match(cl, unique(cl[hc$order]))          # 1, 2, 3 ... in tree order
  names(cl_ord) <- colnames(C)
  # safety check against the membership file written by the run
  saved <- read.csv(file.path(IN_DIR, sprintf("Cluster_membership_k%d.csv", k)),
                    stringsAsFactors = FALSE)
  if (!all(saved$cluster == cl_ord[saved$sample]))
    warning("k = ", k, ": rebuilt clusters differ from Cluster_membership_k", k, ".csv")
  res[[as.character(k)]] <- list(C = C, hc = hc, membership = cl_ord,
                                 coph = coph_tab$cophenetic_correlation[coph_tab$k == k])
}

# ------------------------------ 3. COPHENETIC PLOT ---------------------------
draw_coph <- function() {
  par(mar = c(4.5, 4.5, 3, 1))
  y <- coph_tab$cophenetic_correlation
  plot(coph_tab$k, y, type = "b", pch = 19, lwd = 2, xaxt = "n",
       ylim = range(y) + c(-0.01, 0.01),
       xlab = "Number of clusters (k)", ylab = "Cophenetic correlation",
       main = "NMF package (KL), cophenetic correlation")
  axis(1, at = coph_tab$k)
  text(coph_tab$k, y, labels = sprintf("%.3f", y), pos = 3, cex = 0.9)
  points(coph_tab$k[which.max(y)], max(y), pch = 21, cex = 2.2, bg = NA, col = "firebrick", lwd = 2)
}
save_fig(file.path(IN_DIR, "Cophenetic_coefficient_plot"), 5, 4.5, draw_coph)

# ------------------------------ 4. CONSENSUS HEATMAPS ------------------------
for (k in K_RANGE) {
  r <- res[[as.character(k)]]
  side <- PAL[r$membership]                          # cluster bar, same order as the columns of C
  draw_heat <- function() {
    heatmap(r$C, Rowv = as.dendrogram(r$hc), Colv = as.dendrogram(r$hc), symm = TRUE,
            col = colorRampPalette(c("white", "#B2182B"))(100), scale = "none",
            ColSideColors = side, labRow = NA,
            main = sprintf("Consensus matrix (NMF package, KL), k = %d (cophenetic = %.3f)",
                           k, r$coph))
  }
  save_fig(file.path(IN_DIR, sprintf("Consensus_heatmap_k%d", k)), 7, 7, draw_heat)
}

# ------------------------------ 5. RENUMBER k = 3 TO THE GENEPATTERN LABELS --
pkg <- res[[as.character(K_GP)]]$membership         # package clusters, original sample names

gp <- read.table(GP_FILE, skip = 2, header = TRUE, sep = "\t",
                 check.names = FALSE, quote = "", comment.char = "")
stopifnot("membership.ordered" %in% colnames(gp))
gp_lab <- setNames(as.integer(gp$membership.ordered), clean_id(gp$Name))

stopifnot(setequal(clean_id(names(pkg)), names(gp_lab)))
gp_lab <- gp_lab[clean_id(names(pkg))]               # same sample order as the package result

# For each package cluster, the GenePattern cluster that shares the most samples
overlap <- table(package = pkg, GenePattern = gp_lab)
print(overlap)
map <- apply(overlap, 1, function(r) as.integer(colnames(overlap)[which.max(r)]))
stopifnot(!anyDuplicated(map))                       # one-to-one, otherwise stop and look at the table
cat("Mapping (package cluster -> GenePattern cluster):\n"); print(map)

new_lab <- unname(map[as.character(pkg)])
cat("ARI vs GenePattern:", round(adjustedRandIndex(gp_lab, new_lab), 3), "\n")
cat("Samples in a different group than GenePattern:",
    if (any(new_lab != gp_lab)) paste(names(pkg)[new_lab != gp_lab], collapse = ", ") else "none", "\n")
print(table(new_lab))

out <- data.frame(sample = names(pkg), cluster = new_lab, stringsAsFactors = FALSE)
write.csv(out, file.path(IN_DIR, "Cluster_membership_k3_GPlabels.csv"), row.names = FALSE)
stopifnot(identical(colnames(out), c("sample", "cluster")))
cat("Saved to:", IN_DIR, "\n")
