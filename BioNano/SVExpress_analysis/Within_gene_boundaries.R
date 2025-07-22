```r
#####################################################################################################################
# This can generate the gene-to-sample matrix directly from OverlapGenes (within)
# May miss regulatory effects (e.g., enhancer hijacking).

library(rtracklayer)
library(tidyr)
library(dplyr)
library(readr)

# Load the SV CSV
sv <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/SV_df_unqiue.csv", stringsAsFactors = FALSE)

# Extract relevant columns
sv_clean <- sv[, c("sample", "OverlapGenes")]
sv_clean <- sv_clean[sv_clean$OverlapGenes != "", ]  # Remove rows without gene annotations

# Normalize sample names to match other matrices
sv_clean$sample <- gsub("_", "-", sv_clean$sample)

# Rename columns for consistency
colnames(sv_clean) <- c("Sample", "OverlapGenes")

# Split semicolon-separated genes
sv_clean <- separate_rows(sv_clean, OverlapGenes, sep = ";")

# Remove whitespace
sv_clean$OverlapGenes <- trimws(sv_clean$OverlapGenes)

# Create binary matrix: gene X sample
gene_sample_mat <- with(sv_clean, table(OverlapGenes, Sample))
gene_sample_mat[gene_sample_mat > 1] <- 1

# Export to breakpoints.txt
write.table(as.matrix(gene_sample_mat),
            file = "breakpoints.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

head(gene_sample_mat)

SV <- read.table("breakpoints.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Check header 
head(colnames(SV), 10)

#####################################################################################################################
# Generate CNA matrix from Merged cnv File
library(dplyr)
library(tidyr)
library(readr)

# Load the merged CNA CSV
cnv_data <- read_csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation/merged_cnv_data.csv")

# Filter out rows with missing or empty OverlapGenes
cnv_data <- cnv_data %>%
  filter(!is.na(OverlapGenes), OverlapGenes != "")

# Expand semicolon-separated genes into rows
cnv_expanded <- cnv_data %>%
  mutate(OverlapGenes = strsplit(OverlapGenes, ";")) %>%
  unnest(OverlapGenes)

# Aggregate: mean CopyNumber per gene per sample
cna_matrix <- cnv_expanded %>%
  group_by(OverlapGenes, sample) %>%
  summarise(CopyNumber = mean(CopyNumber, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = CopyNumber)

# Convert to matrix
cna_df <- as.data.frame(cna_matrix)
rownames(cna_df) <- cna_df$OverlapGenes
cna_df <- cna_df[, -1]  # remove gene column from data

# Replace NA with 2 (diploid)
cna_df[is.na(cna_df)] <- 2

# Save to CNA.txt
write.table(cna_df, file = "CNA.txt", sep = "\t", quote = FALSE, col.names = NA)

CNA <- read.table("CNA.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Check header 
head(colnames(CNA), 10)

#####################################################################################################################
# Generate Expression matrix

# ------------------------------------------------
# Load and Prepare Expression Matrix (log2CPM with HGNC symbols)
# ------------------------------------------------
log2cpm_full <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/RNASeq/log2cpm_matrix.csv",
                         header = TRUE, check.names = FALSE)

log2cpm_full$Gene_Symbol <- trimws(log2cpm_full$Gene_Symbol)
log2cpm_full <- log2cpm_full[log2cpm_full$Gene_Symbol != "" & !is.na(log2cpm_full$Gene_Symbol), ]

# Aggregate duplicate gene symbols
log2cpm_aggregated <- log2cpm_full %>%
  group_by(Gene_Symbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

log2cpm_aggregated <- as.data.frame(log2cpm_aggregated)
rownames(log2cpm_aggregated) <- log2cpm_aggregated$Gene_Symbol
log2cpm_aggregated$Gene_Symbol <- NULL

write.table(log2cpm_aggregated, file = "expression_matrix_cleaned.txt", sep = "	", quote = FALSE, row.names = TRUE)

log2cpm <- as.matrix(log2cpm_aggregated)

#####################################################################################################################
# Load SV, expression & CNA Matrices followed SV-expression association analysis

SV <- as.matrix(read.table("breakpoints.txt", header = TRUE, row.names = 1, check.names = FALSE))
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
# To Check the presence or validate the breakpoint matrix

# Load the breakpoint matrix
bp_mat <- read.table("breakpoints.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Check if ESR1 is present # Change the gene name
if ("ESR1" %in% rownames(bp_mat)) {
  # Count how many samples have a breakpoint near ESR1
  num_samples <- sum(bp_mat["ESR1", ] == 1)
  cat("Gene ESR1 has SVs in", num_samples, "samples.\n")

  # List sample names
  samples_with_sv <- colnames(bp_mat)[which(bp_mat["ESR1", ] == 1)]
  print(samples_with_sv)
} else {
  cat("Gene ESR1 is not present in the breakpoints matrix.\n")
}
#####################################################################################################################
