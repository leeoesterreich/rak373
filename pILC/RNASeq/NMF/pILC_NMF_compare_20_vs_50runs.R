###############################################################################
# Compare the consensus NMF results of 20 runs and 50 runs (reads existing files only; nothing is re-run)
#
# Inputs  (same file names in both folders):
#   Cophenetic_correlation.csv, Cluster_membership_k{2,3,4}.csv, Consensus_matrix_k{2,3,4}.csv,
#   NMF_best_run_k3_H.csv
# Output  (folder Comparison_20_vs_50runs, no dates in names):
#   Cophenetic_20_vs_50.csv                 stability for each k
#   Cophenetic_coefficient_20_vs_50_plot.png/.pdf   plot of the table above (k = 2, 3, 4)
#   Cluster_comparison_by_k.csv             adjusted Rand index and number of changed samples per k
#   Cluster_membership_20_vs_50_k3.csv      per sample: cluster in both results, changed or not
#   Sample_consensus_change_20_vs_50_k3.csv own vs other-cluster support in both results
#   Consensus_20_vs_50_k3_scatter.png/.pdf  + Consensus_20_vs_50_k3_correlation.csv
#   Best_run_program_20_vs_50_k3.csv        best-run program of each sample in both results
###############################################################################
suppressPackageStartupMessages(library(ggplot2))
SEED <- 42; set.seed(SEED)

# ------------------------------ 1. Paths -------------------------------------
root  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/NMF"
d20   <- file.path(root, "R_NMFConsensus_20runs")
d50   <- file.path(root, "R_NMFConsensus_50runs")
out   <- file.path(root, "R_NMFConsensus_50runs", "Comparison_20_vs_50runs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cluster_cols <- c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3")


# Read a CSV and say which file is the problem if it is missing or empty
rd <- function(f, ...) {
  if (!file.exists(f)) stop("File not found: ", f)
  if (file.info(f)$size == 0) stop("File is empty (0 bytes): ", f)
  tryCatch(read.csv(f, ...), error = function(e) stop("Cannot read ", f, " -> ", conditionMessage(e)))
}

norm_id <- function(x) {                      # TP09.9_S34 / TP09-9 -> TP09-09
  x <- gsub("\\.", "-", sub("_S[0-9]+$", "", x))
  sub("^(TP[0-9]+)-([0-9])$", "\\1-0\\2", x)
}
# Adjusted Rand index: 1 = identical groups, ~0 = chance; ignores cluster numbering
adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); c2 <- function(x) x * (x - 1) / 2
  sij <- sum(c2(tab)); sa <- sum(c2(rowSums(tab))); sb <- sum(c2(colSums(tab)))
  e <- sa * sb / c2(n)
  if (isTRUE(all.equal(0.5 * (sa + sb), e))) return(1)
  (sij - e) / (0.5 * (sa + sb) - e)
}
read_mem <- function(dir, k) {
  m <- rd(file.path(dir, sprintf("Cluster_membership_k%d.csv", k)))
  data.frame(Key = norm_id(m$sample), Cluster = m$cluster)
}
read_cons <- function(dir, k) {
  C <- as.matrix(rd(file.path(dir, sprintf("Consensus_matrix_k%d.csv", k)), row.names = 1, check.names = FALSE))
  rownames(C) <- colnames(C) <- norm_id(colnames(C)); C
}

# ------------------------------ 2. Cophenetic correlation --------------------
# Recomputed from the consensus matrices with the same formula as the NMF script
# (average-linkage tree on 1 - consensus), so the Cophenetic_correlation.csv files are not needed.
coph_from <- function(dir, k) {
  d <- as.dist(1 - read_cons(dir, k))
  round(cor(as.vector(d), as.vector(cophenetic(hclust(d, method = "average")))), 4)
}
coph <- data.frame(k = 2:4,
                   cophenetic_20runs = sapply(2:4, function(k) coph_from(d20, k)),
                   cophenetic_50runs = sapply(2:4, function(k) coph_from(d50, k)))
coph$difference <- coph$cophenetic_50runs - coph$cophenetic_20runs
write.csv(coph, file.path(out, "Cophenetic_20_vs_50.csv"), row.names = FALSE)
print(coph)

# Plot: cophenetic coefficient for k = 2, 3, 4, 20 runs vs 50 runs (table above is the data behind it)
coph_long <- rbind(data.frame(k = coph$k, Runs = "50 runs", Coefficient = coph$cophenetic_50runs),
                   data.frame(k = coph$k, Runs = "20 runs", Coefficient = coph$cophenetic_20runs))
coph_long$Runs <- factor(coph_long$Runs, levels = c("50 runs", "20 runs"))
p_coph <- ggplot(coph_long, aes(factor(k), Coefficient, group = Runs, color = Runs, linetype = Runs)) +
  geom_line(linewidth = 1) + geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.3f", Coefficient)), vjust = -1.2, size = 4.5, show.legend = FALSE) +
  scale_color_manual(values = c("50 runs" = "#B2182B", "20 runs" = "#7F7F7F")) +
  scale_linetype_manual(values = c("50 runs" = "solid", "20 runs" = "dashed")) +
  scale_y_continuous(limits = c(min(coph_long$Coefficient) - 0.03, 1.005), expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "Cophenetic coefficient: 20 vs 50 runs", x = "Number of clusters (k)", y = "Cophenetic coefficient") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18), legend.title = element_blank(),
        axis.line = element_line(linewidth = 0.8))
ggsave(file.path(out, "Cophenetic_coefficient_20_vs_50_plot.png"), p_coph, width = 6.5, height = 5, dpi = 300)
ggsave(file.path(out, "Cophenetic_coefficient_20_vs_50_plot.pdf"), p_coph, width = 6.5, height = 5)

# ------------------------------ 3. Cluster labels for every k ----------------
by_k <- do.call(rbind, lapply(2:4, function(k) {
  a <- read_mem(d20, k); b <- read_mem(d50, k)
  b <- b[match(a$Key, b$Key), ]
  stopifnot(!anyNA(b$Cluster))
  tab <- table(New = b$Cluster, Old = a$Cluster)
  map <- apply(tab, 1, which.max)                          # new label -> old label with most overlap
  changed <- as.character(a$Cluster) != colnames(tab)[map][b$Cluster]
  data.frame(k = k, adjusted_Rand_index = round(adj_rand(a$Cluster, b$Cluster), 3),
             samples_changed = sum(changed), n_samples = nrow(a))
}))
write.csv(by_k, file.path(out, "Cluster_comparison_by_k.csv"), row.names = FALSE)
print(by_k)

# ------------------------------ 4. k = 3: labels per sample ------------------
m20 <- read_mem(d20, 3); m50 <- read_mem(d50, 3); m50 <- m50[match(m20$Key, m50$Key), ]
tab3 <- table(New = m50$Cluster, Old = m20$Cluster); print(tab3)
map3 <- apply(tab3, 1, which.max)
lab <- data.frame(Key = m20$Key, Cluster_20runs = m20$Cluster, Cluster_50runs = m50$Cluster,
                  Cluster_50runs_as_20run_label = as.integer(colnames(tab3)[map3][m50$Cluster]))
lab$Changed <- lab$Cluster_50runs_as_20run_label != lab$Cluster_20runs
write.csv(lab, file.path(out, "Cluster_membership_20_vs_50_k3.csv"), row.names = FALSE)
if (any(lab$Changed)) { message("Samples that changed cluster:"); print(lab[lab$Changed, ]) } else
  message("No sample changed cluster between 20 and 50 runs (k = 3).")

# ------------------------------ 5. k = 3: consensus support per sample -------
C20 <- read_cons(d20, 3); C50 <- read_cons(d50, 3)
C50 <- C50[rownames(C20), colnames(C20)]
cl  <- m20$Cluster[match(rownames(C20), m20$Key)]            # old labels, so the comparison is like for like
own_other <- function(Cm, cl) {
  diag(Cm) <- NA; ks <- sort(unique(cl))
  mc <- sapply(ks, function(k) rowMeans(Cm[, cl == k, drop = FALSE], na.rm = TRUE))
  idx <- cbind(seq_along(cl), match(cl, ks)); own <- mc[idx]; mc[idx] <- -Inf
  data.frame(own = own, other = apply(mc, 1, max), other_cluster = ks[max.col(mc, ties.method = "first")])
}
o20 <- own_other(C20, cl); o50 <- own_other(C50, cl)
chg <- data.frame(Key = rownames(C20), Cluster = cl,
                  own_consensus_20 = o20$own, own_consensus_50 = o50$own,
                  strongest_other_cluster_20 = o20$other_cluster, strongest_other_cluster_50 = o50$other_cluster,
                  other_consensus_20 = o20$other, other_consensus_50 = o50$other,
                  other_minus_own_20 = o20$other - o20$own, other_minus_own_50 = o50$other - o50$own,
                  mean_abs_change_in_consensus = rowMeans(abs(C50 - C20)))
chg <- chg[order(-chg$other_minus_own_50), ]                  # most ambiguous samples first
write.csv(chg, file.path(out, "Sample_consensus_change_20_vs_50_k3.csv"), row.names = FALSE)
print(round(head(chg[, c("own_consensus_20", "own_consensus_50", "other_minus_own_20", "other_minus_own_50")], 8), 3))
print(head(chg$Key, 8))

# ------------------------------ 6. Scatter: consensus values, 20 vs 50 -------
up <- upper.tri(C20)
sc <- data.frame(c20 = C20[up], c50 = C50[up])
cor_tab <- data.frame(Measure = c("Pearson correlation of consensus values (sample pairs)",
                                  "Mean absolute difference", "Largest absolute difference"),
                      Value = round(c(cor(sc$c20, sc$c50), mean(abs(sc$c20 - sc$c50)), max(abs(sc$c20 - sc$c50))), 3))
write.csv(cor_tab, file.path(out, "Consensus_20_vs_50_k3_correlation.csv"), row.names = FALSE)
print(cor_tab)
p <- ggplot(sc, aes(c20, c50)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(alpha = 0.35, size = 1.5) +
  labs(title = "Consensus values, k = 3", x = "Consensus (20 runs)", y = "Consensus (50 runs)",
       subtitle = sprintf("Pearson r = %.3f", cor(sc$c20, sc$c50))) +
  theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out, "Consensus_20_vs_50_k3_scatter.png"), p, width = 6, height = 6, dpi = 300)
ggsave(file.path(out, "Consensus_20_vs_50_k3_scatter.pdf"), p, width = 6, height = 6)

# ------------------------------ 7. Best-run program of each sample -----------
# (the W/H used for the marker-gene heatmap comes from the best run, so check it as well)
bestH <- function(dir) {
  H <- as.matrix(rd(file.path(dir, "NMF_best_run_k3_H.csv"), row.names = 1, check.names = FALSE))
  colnames(H) <- norm_id(colnames(H)); max.col(t(H), ties.method = "first") -> pr; names(pr) <- colnames(H); pr
}
p20 <- bestH(d20); p50 <- bestH(d50)[names(p20)]
tabp <- table(New = p50, Old = p20); print(tabp)
mapp <- apply(tabp, 1, which.max)
bp <- data.frame(Key = names(p20), Program_20runs = p20, Program_50runs = p50,
                 Program_50runs_as_20run_label = as.integer(colnames(tabp)[mapp][p50]), row.names = NULL)
bp$Changed <- bp$Program_50runs_as_20run_label != bp$Program_20runs
bp$Final_cluster_20runs <- m20$Cluster[match(bp$Key, m20$Key)]
write.csv(bp, file.path(out, "Best_run_program_20_vs_50_k3.csv"), row.names = FALSE)

message("Done. Output in ", out)
