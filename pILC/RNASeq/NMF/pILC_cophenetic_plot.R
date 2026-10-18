###############################################################################
# Cophenetic coefficient plot for k = 2, 3, 4 -- 50-run consensus NMF only
# (the 20 vs 50 comparison plot is made by pILC_NMF_compare_20_vs_50runs.R)
# Values are recomputed from the consensus matrices with the same formula as the NMF script:
#   d = 1 - consensus;  tree = hclust(d, "average");  coefficient = cor(d, cophenetic(tree))
# Output (folder R_NMFConsensus_50runs, no dates in names):
#   Cophenetic_coefficient_plot.png / .pdf
#   Cophenetic_coefficient_by_k.csv        (the numbers behind the plot)
###############################################################################
suppressPackageStartupMessages(library(ggplot2))
SEED <- 42; set.seed(SEED)

# ------------------------------ 1. Settings ----------------------------------
root    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/NMF"
d50     <- file.path(root, "R_NMFConsensus_50runs")
out_dir <- d50
K_RANGE <- 2:4

# ------------------------------ 2. Cophenetic coefficient per k --------------
coph_from <- function(dir, k) {
  f <- file.path(dir, sprintf("Consensus_matrix_k%d.csv", k))
  if (!file.exists(f) || file.info(f)$size == 0) stop("Missing or empty file: ", f)
  C <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  d <- as.dist(1 - C)
  round(cor(as.vector(d), as.vector(cophenetic(hclust(d, method = "average")))), 4)
}
tab <- data.frame(k = K_RANGE, cophenetic_50runs = sapply(K_RANGE, function(k) coph_from(d50, k)))
write.csv(tab, file.path(out_dir, "Cophenetic_coefficient_by_k.csv"), row.names = FALSE)
print(tab)

# ------------------------------ 3. Plot --------------------------------------
p <- ggplot(tab, aes(factor(k), cophenetic_50runs, group = 1)) +
  geom_line(linewidth = 1, color = "#B2182B") +
  geom_point(size = 4, color = "#B2182B") +
  geom_text(aes(label = sprintf("%.3f", cophenetic_50runs)), vjust = -1.2, size = 4.5) +
  scale_y_continuous(limits = c(min(tab$cophenetic_50runs) - 0.03, 1.005),
                     expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "Cophenetic coefficient (50 runs)", x = "Number of clusters (k)", y = "Cophenetic coefficient") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20), axis.line = element_line(linewidth = 0.8))

ggsave(file.path(out_dir, "Cophenetic_coefficient_plot.png"), p, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "Cophenetic_coefficient_plot.pdf"), p, width = 6, height = 5)
message("Done. Output in ", out_dir)
