```r
# Tried to generate geneXsample matrix by Excel macro file (provided by author) but not working
# one possible reason of not having ori information on bionano data
# So i used my own R script to generate breakpoint matrix
# SVexpress Preprocessing and Integration Pipelin
# Generate CNA and expression matrix as earlier used in SVExpress_TAD and SVExpress_gene_boundaries
#####################################################################################################################

# This can generates the gene-to-sample breakpoint matrix (breakpoints.txt) directly from structural variant and gene annotation files
# SV breakpoints to nearby genes using coordinate overlaps
# sv_bedpe.txt: BEDPE-like file of SV breakpoints (chrom1, start1, end1, chrom2, start2, end2, sample_id), From SV calling
# gene_coords.bed: Gene annotations (chrom, start, end, gene_name), Convert GTF to BED
# window_size: Â±1 Mb from TSS to define "proximity"

######################################################
# ------------------------
# Convert GTF to BED
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

# ------------------------------------------------
# generates the gene-to-sample breakpoint matrix 
# -------------------------------------------------
library(GenomicRanges)

# Load SV CSV
sv_raw <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/SV_df_unique.csv", stringsAsFactors = FALSE)

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

###################### Based on extended window size on both side ######################

# Extend gene windows ±10kb, ±100kb, ±1Mb # Adjust window size accordingly
window_size <- 1e5
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

# Get gene sample hit table
sv_hits <- rbind(
  data.frame(gene = mcols(flank_gene_gr)$gene[subjectHits(overlap1)],
             sample = mcols(sv_gr1)$sample[queryHits(overlap1)]),
  data.frame(gene = mcols(flank_gene_gr)$gene[subjectHits(overlap2)],
             sample = mcols(sv_gr2)$sample[queryHits(overlap2)])
)

# Create binary gene-by-sample matrix 
gene_sample_mat <- with(sv_hits, table(gene, sample))
gene_sample_mat[gene_sample_mat > 1] <- 1

# Export # Adjust file name acoording to window size
write.table(as.matrix(gene_sample_mat), "breakpoints_10kb.txt", sep = "\t", quote = FALSE)

head(gene_sample_mat)

SV <- read.table("breakpoints_10kb.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# Check header 
head(colnames(SV), 10)


###################### Based on extended window size on either side ######################
library(GenomicRanges)

# Extend upstream and downstream separately (100kb) # Adjust size accordingly
upstream_100kb_gr <- flank(gene_gr, width = 100000, start = TRUE, both = FALSE)
downstream_100kb_gr <- flank(gene_gr, width = 100000, start = FALSE, both = FALSE)

# Convert SV breakpoints to GRanges
sv_gr1 <- GRanges(seqnames = sv_bedpe$chr1,
                  ranges = IRanges(start = sv_bedpe$start1, end = sv_bedpe$end1),
                  sample = sv_bedpe$sample)

sv_gr2 <- GRanges(seqnames = sv_bedpe$chr2,
                  ranges = IRanges(start = sv_bedpe$start2, end = sv_bedpe$end2),
                  sample = sv_bedpe$sample)

# --------------------------
# For upstream overlaps
# --------------------------
overlap1_up <- findOverlaps(sv_gr1, upstream_100kb_gr)
overlap2_up <- findOverlaps(sv_gr2, upstream_100kb_gr)

sv_hits_up <- rbind(
  data.frame(gene = mcols(upstream_100kb_gr)$gene[subjectHits(overlap1_up)],
             sample = mcols(sv_gr1)$sample[queryHits(overlap1_up)]),
  data.frame(gene = mcols(upstream_100kb_gr)$gene[subjectHits(overlap2_up)],
             sample = mcols(sv_gr2)$sample[queryHits(overlap2_up)])
)

gene_sample_mat_up <- with(sv_hits_up, table(gene, sample))
gene_sample_mat_up[gene_sample_mat_up > 1] <- 1

write.table(as.matrix(gene_sample_mat_up), 
            file = "breakpoints_upstream100kb.txt", 
            sep = "\t", 
            quote = FALSE)

# --------------------------
# For downstream overlaps
# --------------------------
overlap1_down <- findOverlaps(sv_gr1, downstream_100kb_gr)
overlap2_down <- findOverlaps(sv_gr2, downstream_100kb_gr)

sv_hits_down <- rbind(
  data.frame(gene = mcols(downstream_100kb_gr)$gene[subjectHits(overlap1_down)],
             sample = mcols(sv_gr1)$sample[queryHits(overlap1_down)]),
  data.frame(gene = mcols(downstream_100kb_gr)$gene[subjectHits(overlap2_down)],
             sample = mcols(sv_gr2)$sample[queryHits(overlap2_down)])
)

gene_sample_mat_down <- with(sv_hits_down, table(gene, sample))
gene_sample_mat_down[gene_sample_mat_down > 1] <- 1

write.table(as.matrix(gene_sample_mat_down), 
            file = "breakpoints_downstream100kb.txt", 
            sep = "\t", 
            quote = FALSE)

#####################################################################################################################
# Load SV, expression & CNA Matrices followed SV-expression association analysis

SV <- as.matrix(read.table("breakpoints_downstream10kb.txt", header = TRUE, row.names = 1, check.names = FALSE))
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

# Save significant genes to file
write.table(sig_genes,
            file = "lm_test_results_significant_genes.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# Filter genes in CNA correction model where CNA correction p-value is < 0.05
sig_genes <- subset(results,
                    `CNA correction (p-value)` < 0.05)
                    

# Save significant genes to file
write.table(sig_genes,
            file = "lm_test_results_sig_genes_CNA_model.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)
			
# Filter genes in No correction model where no correction p-value is < 0.05
sig_genes <- subset(results,
                    `no correction (p-value)` < 0.05)
                    
# Save significant genes to file
write.table(sig_genes,
            file = "lm_test_results_sig_genes_No_Correction_model.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)
