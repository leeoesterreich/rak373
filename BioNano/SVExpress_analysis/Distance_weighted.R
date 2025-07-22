```r
# Generate a Distance-Weighted Breakpoint Matrix, emphasizes closer breakpoints
# Consider a ±1Mb region around the gene’s TSS/start coordinate.
# If a breakpoint exists within the ±1Mb window, score =log2(distance in bp+1)
# If no breakpoint exists in window, score=log2(1,000,000+1)

library(GenomicRanges)
library(data.table)
library(dplyr)

# -------------------------
# Load SV breakpoints
# -------------------------
sv_raw <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/SV_df_unique.csv", stringsAsFactors = FALSE)
sv_raw$sample <- gsub("_", "-", sv_raw$sample)

# Create two sets of breakpoints (one from each end)
sv_ends1 <- GRanges(seqnames = sv_raw$RefcontigID1,
                    ranges = IRanges(start = sv_raw$RefStartPos, width = 1),
                    sample = sv_raw$sample)

sv_ends2 <- GRanges(seqnames = sv_raw$RefcontigID2,
                    ranges = IRanges(start = sv_raw$RefEndPos, width = 1),
                    sample = sv_raw$sample)

# Combine both ends
sv_breakpoints <- c(sv_ends1, sv_ends2)

# -------------------------
# Load gene coordinates
# -------------------------
gene_coords <- read.table("gene_coords.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_coords) <- c("chr", "start", "end", "gene")

# Warn if duplicates are found
dup_genes <- gene_coords$gene[duplicated(gene_coords$gene)]
if (length(dup_genes) > 0) {
  message("Warning: Duplicate gene symbols found. They will be collapsed by mean distance.")
}

# Use gene TSS
gene_tss <- GRanges(seqnames = gene_coords$chr,
                    ranges = IRanges(start = gene_coords$start, width = 1),
                    gene = gene_coords$gene)

# -------------------------
# Initialize Distance Matrix
# -------------------------
genes <- gene_coords$gene
samples <- unique(sv_raw$sample)
log_default <- log2(1e6 + 1)

dist_matrix <- matrix(log_default, nrow = length(genes), ncol = length(samples),
                      dimnames = list(genes, samples))

# -------------------------
# Compute Distance-weighted Scores
# -------------------------
for (i in seq_along(gene_tss)) {
  gene_gr <- gene_tss[i]
  gene_name <- gene_gr$gene

  # ±1Mb window around TSS
  gene_window <- promoters(gene_gr, upstream = 1e6, downstream = 1e6)

  # Find overlapping SVs
  overlaps <- findOverlaps(gene_window, sv_breakpoints, ignore.strand = TRUE)
  if (length(overlaps) == 0) next

  matched_svs <- sv_breakpoints[subjectHits(overlaps)]

  dt <- data.table(
    sample = matched_svs$sample,
    dist = abs(start(matched_svs) - start(gene_gr))
  )
  closest <- dt[, .(min_dist = min(dist)), by = sample]
  closest$score <- log2(closest$min_dist + 1)

  dist_matrix[gene_name, closest$sample] <- closest$score

  if (i %% 100 == 0) cat("Processed", i, "of", length(genes), "genes\n")
}

# -------------------------
# Merge duplicate genes (average distances)
# -------------------------
dist_df <- as.data.frame(dist_matrix)
dist_df$Gene <- rownames(dist_df)

dist_df_collapsed <- dist_df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE), .groups = "drop") # More robust

# If want to change the aggregation logic from median to min/ max/ mean
# summarise(across(where(is.numeric), min, na.rm = TRUE), .groups = "drop") # For strongest SV effect
# summarise(across(where(is.numeric), max, na.rm = TRUE), .groups = "drop") # For farthest effect
# summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") # Computing the mean distance score for each sample for that gene

# Final matrix
final_matrix <- as.matrix(dist_df_collapsed[, -1])
rownames(final_matrix) <- dist_df_collapsed$Gene

# -------------------------
# Save to file
# -------------------------
write.table(final_matrix, file = "distance_weighted_breakpoint_matrix.txt", sep = "\t",
            quote = FALSE, col.names = NA)

head(final_matrix)

#####################################################################################################################
# Generate Expression and CNA matrix as earlier
#####################################################################################################################

# Load SV, expression & CNA Matrices followed SV-expression association analysis

SV <- as.matrix(read.table("distance_weighted_breakpoint_matrix.txt", header = TRUE, row.names = 1, check.names = FALSE))
copy <- as.matrix(read.table("CNA.txt", header = TRUE, row.names = 1, check.names = FALSE))
log2cpm <- as.matrix(read.table("expression_matrix_cleaned.txt", header = TRUE, row.names = 1, check.names = FALSE))

# ------------------------
# Align All Matrices by Shared Genes and Samples
# ------------------------
shared_samples <- Reduce(intersect, list(colnames(log2cpm), colnames(SV), colnames(copy)))
shared_genes <- Reduce(intersect, list(rownames(log2cpm), rownames(SV), rownames(copy)))

expression <- log2cpm[shared_genes, shared_samples, drop = FALSE]
SV <- SV[shared_genes, shared_samples, drop = FALSE]
copy <- copy[shared_genes, shared_samples, drop = FALSE]

cat("Final matrix dimensions (genes x samples):\n")
cat("Expression:", dim(expression), "\n")
cat("SV:", dim(SV), "\n")
cat("CNA:", dim(copy), "\n")

# ------------------------
# Linear Model Integration
# ------------------------
SV_integration <- function(expression, copy, SV) {
  cat("Starting time:", date(), "
")
  J <- ncol(expression)
  I <- nrow(expression)
  beta <- matrix(NA, nrow = I, ncol = 4)

  for (i in 1:I) {
    if (sd(SV[i, ], na.rm = TRUE) > 0.001) {
      tryCatch({
        z1 <- lm(as.numeric(expression[i, ]) ~ as.numeric(SV[i, ]))
        z1_coef <- summary(z1)$coefficients
        if (nrow(z1_coef) >= 2 && !anyNA(z1_coef[2, 3:4])) {
          beta[i, 1] <- z1_coef[2, 3]
          beta[i, 2] <- z1_coef[2, 4]
        }

        z2 <- lm(as.numeric(expression[i, ]) ~ as.numeric(SV[i, ]) + as.numeric(copy[i, ]))
        z2_coef <- summary(z2)$coefficients
        if (nrow(z2_coef) >= 2 && !anyNA(z2_coef[2, 3:4])) {
          beta[i, 3] <- z2_coef[2, 3]
          beta[i, 4] <- z2_coef[2, 4]
        }
      }, error = function(e) {
        cat("Model failed for gene", rownames(expression)[i], "-", conditionMessage(e), "
")
      })
    }

    if (i %% 500 == 0) cat("Processed gene", i, "of", I, "
")
  }

  cat("Ending time:", date(), "
")
  return(beta)
}

# ------------------------
# Run Integration and Save Results
# ------------------------
beta <- SV_integration(expression, copy, SV)
rownames(beta) <- rownames(expression)
colnames(beta) <- c("no correction (t-stat)", "no correction (p-value)",
                    "CNA correction (t-stat)", "CNA correction (p-value)")

write.table(beta, file = "lm_test_results.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# ------------------------
# Filter selected Genes
# ------------------------
# Load full results
results <- read.table("lm_test_results.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Define gene list
genes_of_interest <- c("ESR1", "ARNT2", "TAB3", "DMD", "AP2B1", 
                       "SKAP1", "EDC3", "PREX1", "ARHGEF7", 
                       "ATP11A", "DENND5A", "SORL1")

# Subset the results
selected_genes <- results[rownames(results) %in% genes_of_interest, ]

# Save to file
write.table(selected_genes, file = "selected_genes_output.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Optionally, print it in the console
print(selected_genes)

# ------------------------
# Filter Significant Genes (Raw p-value < 0.05)
# ------------------------
results <- read.table("lm_test_results.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Filter genes in either model where either p-value is < 0.05
sig_genes <- subset(results,
                    `no correction (p-value)` < 0.05 |
                    `CNA correction (p-value)` < 0.05)

# Save significant genes 
write.table(sig_genes,
            file = "lm_test_results_significant_genes.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Filter genes in CNA correction model where CNA correction p-value is < 0.05
sig_genes <- subset(results,
                    `CNA correction (p-value)` < 0.05)
                    

# Save significant genes 
write.table(sig_genes,
            file = "lm_test_results_sig_genes_CNA_model.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)
			
# Filter genes in No correction model where no correction p-value is < 0.05
sig_genes <- subset(results,
                    `no correction (p-value)` < 0.05)
                    
# Save significant genes 
write.table(sig_genes,
            file = "lm_test_results_sig_genes_No_Correction_model.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

#####################################################################################################################
