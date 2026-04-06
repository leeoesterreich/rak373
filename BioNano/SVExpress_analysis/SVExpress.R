
```r
#####################################################################################################################
# This can generate the gene-to-sample matrix directly from OverlapGenes and skips coordinate-based SV-to-gene proximity (±1 Mb).
# May miss regulatory effects (e.g., enhancer hijacking).

library(rtracklayer)
library(tidyr)

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

# Create binary matrix: gene × sample
gene_sample_mat <- with(sv_clean, table(OverlapGenes, Sample))
gene_sample_mat[gene_sample_mat > 1] <- 1

# Export to breakpoints.txt
write.table(as.matrix(gene_sample_mat),
            file = "breakpoints.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

SV <- read.table("breakpoints.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Check header 
head(colnames(SV), 10)

#####################################################################################################################

# This can generates the gene-to-sample breakpoint matrix (breakpoints.txt) directly from structural variant and gene annotation files
# SV breakpoints to nearby genes using coordinate overlaps
# sv_bedpe.txt: BEDPE-like file of SV breakpoints (chrom1, start1, end1, chrom2, start2, end2, sample_id), From SV calling
# gene_coords.bed: Gene annotations (chrom, start, end, gene_name), Convert GTF → BED
# window_size: ±1 Mb from TSS to define "proximity"

######################################################
# ------------------------
# Convert GTF → BED
# ------------------------
library(rtracklayer)

# Set path to your GTF file
gtf_file <- "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/gencode.v47.annotation.gtf"

# Import GTF and extract gene features
gtf <- import(gtf_file)
genes <- gtf[gtf$type == "gene"]

# Create BED-style data frame: chr, start, end, gene_name
bed <- data.frame(
  chr = as.character(seqnames(genes)),
  start = start(genes) - 1,  # BED is 0-based start
  end = end(genes),
  gene = mcols(genes)$gene_name
)

# Save as gene_coords.bed
write.table(bed, "gene_coords.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# --------------------------------------------------
# generates the gene-to-sample breakpoint matrix  
# --------------------------------------------------
library(GenomicRanges)

# Load SV CSV
sv_raw <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/SV_df_unqiue.csv", stringsAsFactors = FALSE)

# Fix sample names: replace underscores with dashes
sv_raw$sample <- gsub("_", "-", sv_raw$sample)

# Create BEDPE-like data frame
sv_bedpe <- data.frame(
  chr1 = sv_raw$RefcontigID1,
  start1 = sv_raw$RefStartPos,
  end1 = sv_raw$RefStartPos + 1,
  chr2 = sv_raw$RefcontigID2,
  start2 = sv_raw$RefEndPos,
  end2 = sv_raw$RefEndPos + 1,
  sample = sv_raw$sample
)

# Load gene coordinates (BED format: chr, start, end, gene name)
gene_coords <- read.table("gene_coords.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_coords) <- c("chr", "start", "end", "gene")

gene_gr <- GRanges(seqnames = gene_coords$chr,
                   ranges = IRanges(start = gene_coords$start, end = gene_coords$end),
                   gene = gene_coords$gene)

# Extend gene windows ±1Mb
window_size <- 1e6
flank_gene_gr <- resize(gene_gr, width = width(gene_gr) + 2 * window_size, fix = "center")

# Convert SV breakpoints to GRanges
sv_gr1 <- GRanges(seqnames = sv_bedpe$chr1,
                  ranges = IRanges(start = sv_bedpe$start1, end = sv_bedpe$end1),
                  sample = sv_bedpe$sample)

sv_gr2 <- GRanges(seqnames = sv_bedpe$chr2,
                  ranges = IRanges(start = sv_bedpe$start2, end = sv_bedpe$end2),
                  sample = sv_bedpe$sample)

# Find overlaps between breakpoints and extended gene regions
overlap1 <- findOverlaps(sv_gr1, flank_gene_gr)
overlap2 <- findOverlaps(sv_gr2, flank_gene_gr)

sv_hits <- rbind(
  data.frame(gene = mcols(flank_gene_gr)$gene[subjectHits(overlap1)],
             sample = mcols(sv_gr1)$sample[queryHits(overlap1)]),
  data.frame(gene = mcols(flank_gene_gr)$gene[subjectHits(overlap2)],
             sample = mcols(sv_gr2)$sample[queryHits(overlap2)])
)

# Create binary gene-by-sample matrix 
gene_sample_mat <- with(sv_hits, table(gene, sample))
gene_sample_mat[gene_sample_mat > 1] <- 1

# Export
write.table(as.matrix(gene_sample_mat), "breakpoints.txt", sep = "\t", quote = FALSE)


SV <- read.table("breakpoints.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Check header 
head(colnames(SV), 10)

#####################################################################################################################
# Generate expression.txt from log2CPM.csv
# ------------------------
# Load expression and SV
# ------------------------
# Load expression matrix (log2CPM) from CSV
log2cpm_raw <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/RNASeq/log2cpm_matrix.csv",
                        header = TRUE, row.names = 1, check.names = FALSE)

# Drop the empty second column if present (common issue)
if (colnames(log2cpm_raw)[1] == "") {
  log2cpm <- log2cpm_raw[, -1]
} else {
  log2cpm <- log2cpm_raw
}

# Remove extra 'external_gene_name' column if it was included
log2cpm <- log2cpm[, colnames(log2cpm) != "external_gene_name"]

# Convert to matrix
log2cpm <- as.matrix(log2cpm)

# Load breakpoints matrix
SV <- read.table("breakpoints.txt", header = TRUE, row.names = 1, sep = "\t")
SV <- as.matrix(SV)

# ------------------------------------------------------------------------
# Define filtering function
# Keep genes present in both datasets (expression and SV matrices)
# and with variation in SV (i.e., SD > 0.001 across samples)
# ------------------------------------------------------------------------

filter_genes_by_sv_presence <- function(expression, SV, min_sd = 0.001) {
  common_genes <- intersect(rownames(expression), rownames(SV))
  variable_genes <- common_genes[apply(SV[common_genes, , drop = FALSE], 1, sd) > min_sd]
  
  expression_filtered <- expression[variable_genes, , drop = FALSE]
  SV_filtered <- SV[variable_genes, , drop = FALSE]
  
  list(expression = expression_filtered, SV = SV_filtered)
}

# ------------------------
# Apply filtering
# ------------------------
filtered_data <- filter_genes_by_sv_presence(log2cpm, SV)

expression_filtered <- filtered_data$expression
SV_filtered <- filtered_data$SV

# Optional: sort by gene name
expression_filtered <- expression_filtered[order(rownames(expression_filtered)), ]

# ------------------------
# Save final expression matrix as TSV
# ------------------------
write.table(expression_filtered,
            file = "expression.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

dim(expression)        # Should be genes × samples
head(rownames(expression))  # Should show gene IDs
expression <- read.table("expression.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Check header 
head(colnames(expression), 10)

#####################################################################################################################
# Generate CNA.txt from Merged cnv File
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
# Ensuring that sample names align across all input matrices (expression, SV, CNA) is critical before running the SV_integration function
# Load all matrices:
expression <- as.matrix(read.table("expression.txt", header = TRUE, row.names = 1))
SV <- as.matrix(read.table("breakpoints.txt", header = TRUE, row.names = 1))
CNA <- as.matrix(read.table("CNA.txt", header = TRUE, row.names = 1))

# Check sample name consistency:
all_samples <- Reduce(intersect, list(colnames(expression), colnames(SV), colnames(CNA)))
length(all_samples)  # Should be equal to number of shared samples

# Subset all matrices to shared samples in the same order:
expression <- expression[, all_samples]
SV <- SV[, all_samples]
CNA <- CNA[, all_samples]

# Recheck dimensions
dim(expression)
dim(SV)
dim(CNA)

#####################################################################################################################
# ------------------------
# Load matrices
# ------------------------
expression <- as.matrix(read.table("expression.txt", header = TRUE, row.names = 1))
SV <- as.matrix(read.table("breakpoints.txt", header = TRUE, row.names = 1))
copy <- as.matrix(read.table("CNA.txt", header = TRUE, row.names = 1))

# -----------------------------------------------------
# Align by shared sample names across all 3 matrices 
# -----------------------------------------------------
shared_samples <- Reduce(intersect, list(colnames(expression), colnames(SV), colnames(copy)))
expression <- expression[, shared_samples]
SV <- SV[, shared_samples]
copy <- copy[, shared_samples]

# ------------------------
# SV Integration Function 
# ------------------------
SV_integration <- function(expression, copy, SV) {
  cat("Starting time:", date(), "\n")
  J <- ncol(expression)
  I <- nrow(expression)
  beta <- matrix(-1, nrow = I, ncol = 4)
  
  for (i in 1:I) {
    if (sd(SV[i, ]) > 0.001) {
      # SV Only Model, No correction, Tests if expression is different in samples with vs without an SV near the gene
      z1 <- lm(as.numeric(expression[i, ]) ~ as.numeric(SV[i, ]))
      beta[i, 1] <- summary(z1)$coefficients[2, 3]  # t-statistic
      beta[i, 2] <- summary(z1)$coefficients[2, 4]  # p-value
      
      # SV + CNA Model, CNA correction, Controls for the possibility that changes in gene expression 
	  # might be due to CNAs (gene dosage) rather than SVs
      z2 <- lm(as.numeric(expression[i, ]) ~ as.numeric(SV[i, ]) + as.numeric(copy[i, ]))
      beta[i, 3] <- summary(z2)$coefficients[2, 3]
      beta[i, 4] <- summary(z2)$coefficients[2, 4]
    }
    
    if (i %% 500 == 0) cat("Processed gene", i, "of", I, "\n")
  }
  
  cat("Ending time:", date(), "\n")
  return(beta)
}

# ------------------------
# Run and Save Results
# ------------------------
beta <- SV_integration(expression, copy, SV)
rownames(beta) <- rownames(expression)
colnames(beta) <- c("no correction (t-stat)", "no correction (p-value)", 
                    "CNA correction (t-stat)", "CNA correction (p-value)")

write.table(beta, file = "lm_test_results.txt", sep = "\t", quote = FALSE)

# ------------------------
# Annotation
# ------------------------

# Load required libraries
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load your results
results <- read.table("lm_test_results.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Extract Ensembl gene IDs
ensembl_ids <- rownames(results)

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add to results
results$Gene_Symbol <- gene_symbols

# Reorder columns (optional)
results <- results[, c("Gene_Symbol", colnames(results)[1:(ncol(results) - 1)])]

# Save annotated results
write.table(results, file = "lm_test_results_annotated.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# -----------------------------
# Apply FDR Correction (BH)
# -----------------------------

# Load the annotated results file
results <- read.table("lm_test_results_annotated.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Apply FDR correction to both p-value columns
results$FDR_no_correction <- p.adjust(results$`no correction (p-value)`, method = "BH")
results$FDR_CNA_correction <- p.adjust(results$`CNA correction (p-value)`, method = "BH")

# Save the updated results with both FDR columns
write.table(results, file = "lm_test_results_annotated_with_BH_fdr.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# ------------------------------------------------------------------
# Filter genes where the p-value is less than 0.05 in either model.
# ------------------------------------------------------------------
# Load the annotated results
results <- read.table("lm_test_results_annotated.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter genes with p-value < 0.05 in either model
sig_genes <- subset(results, 
                    `no.correction..p.value.` < 0.05 | `CNA.correction..p.value.` < 0.05)

# Sort by minimum p-value
sig_genes <- sig_genes[order(pmin(sig_genes$`no.correction..p.value.`, 
                                  sig_genes$`CNA.correction..p.value.`, 
                                  na.rm = TRUE)), ]

# Save the filtered results
write.table(sig_genes, file = "significant_genes_rawp.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------
# Filter significant genes based on FDR values
# ------------------------------------------------
# Load libraries
library(stats)

# Load the annotated results
results <- read.table("lm_test_results_annotated.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Apply FDR correction for both p-value columns
results$FDR_no_correction <- p.adjust(results$`no.correction..p.value.`, method = "fdr")
results$FDR_CNA_correction <- p.adjust(results$`CNA.correction..p.value.`, method = "fdr")

# Filter genes with FDR < 0.05 in either model
sig_genes <- subset(results, FDR_no_correction < 0.05 | FDR_CNA_correction < 0.05)

# Sort by the most significant FDR
sig_genes <- sig_genes[order(pmin(sig_genes$FDR_no_correction, sig_genes$FDR_CNA_correction, na.rm = TRUE)), ]

# Save significant genes to file
write.table(sig_genes, file = "significant_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################################################################################################

```

