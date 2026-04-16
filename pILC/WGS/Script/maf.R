# Load required libraries
library(maftools)
library(data.table)

# Define input/output directories
input_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/PCGR/pcgr_outputs_2.2.5"
output_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/PCGR/MAFtool"

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Find all MAF files recursively
maf_files <- list.files(input_dir, pattern = "\\.maf$", recursive = TRUE, full.names = TRUE)

# Confirm what files were found
cat("Found MAF files:\n")
print(maf_files)

# Read and combine all MAFs into a single data.table
maf_list <- lapply(maf_files, function(f) {
  tryCatch({
    fread(f, data.table = FALSE)
  }, error = function(e) {
    cat("Error reading file:", f, "\n", e$message, "\n")
    return(NULL)
  })
})

# Remove any NULL entries (failed reads)
maf_list <- Filter(Negate(is.null), maf_list)

# Combine into one big MAF data.frame
combined_maf_df <- do.call(rbind, maf_list)

# Now read it into a maftools MAF object
combined_maf <- read.maf(maf = combined_maf_df)

saveRDS(combined_maf, file = "combined_maf.rds")

combined_maf <- readRDS("combined_maf.rds")


# Plot summary dashboard for all 10 samples together
pdf(file.path(output_dir, "Cohort_MAF_Summary.pdf"), width = 10, height = 8)
plotmafSummary(maf = combined_maf,
               rmOutlier = TRUE,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE)
dev.off()

# Save as PNG
png(file.path(output_dir, "Cohort_MAF_Summary.png"), width = 1000, height = 800, res = 150)
plotmafSummary(maf = combined_maf,
               rmOutlier = TRUE,
               addStat = 'median',
               dashboard = TRUE,
               titvRaw = FALSE)
dev.off()

# Save summary tables
# write.mafSummary(maf = combined_maf, basename = file.path(output_dir, "Cohort"))

#-------------------------------------------
# Plot Transition and Transversions.
#-------------------------------------------
# Calculate TiTv ratios
combined_maf.titv <- titv(maf = combined_maf, plot = FALSE, useSyn = TRUE)

# Plot TiTv summary
plotTiTv(res = combined_maf.titv)

# PNG
png(file.path(output_dir, "Cohort_TiTv_Plot.png"), width = 1000, height = 800, res = 150)
plotTiTv(res = combined_maf.titv)
dev.off()

# PDF
pdf(file.path(output_dir, "Cohort_TiTv_Plot.pdf"), width = 10, height = 8)
plotTiTv(res = combined_maf.titv)
dev.off()

#-------------------------------------------
# Check Available Columns
#-------------------------------------------
getFields(combined_maf)

#-------------------------------------------
# check how many samples mutated CDH1
#-------------------------------------------
cdh1_mutations <- subset(getGeneSummary(combined_maf), Hugo_Symbol == "CDH1")
print(cdh1_mutations)

cdh1_variants <- subset(combined_maf@data, Hugo_Symbol == "CDH1")
unique(cdh1_variants$Tumor_Sample_Barcode)


# Extract all mutation records for CDH1
cdh1_variants <- subset(combined_maf@data, Hugo_Symbol == "CDH1")

# Check how many unique samples have CDH1 mutations
cdh1_samples <- unique(cdh1_variants$Tumor_Sample_Barcode)
cat("Number of samples with CDH1 mutations:", length(cdh1_samples), "\n")

# Optional: print sample names
print(cdh1_samples)

# Save to CSV for downstream reference

write.csv(
  cdh1_variants,
  "CDH1_mutations_all_samples.csv",
  row.names = FALSE
)

# Quick summary: mutation types in CDH1

table(cdh1_variants$Variant_Classification)
table(cdh1_variants$Variant_Type)

# Table with just the key features
library(dplyr)

cdh1_summary <- cdh1_variants %>%
  select(
    Tumor_Sample_Barcode,
    Chromosome,
    Start_Position,
    Variant_Classification,
    Variant_Type,
    HGVSc,
    HGVSp_Short,
    dbSNP_RS,
    Existing_variation,
    CLIN_SIG,
    IMPACT,
    t_depth, t_ref_count, t_alt_count,
    gnomAD_AF,
    DOMAINS
  )

write.csv(cdh1_summary, "CDH1_summary_for_presentation.csv", row.names = FALSE)

#-------------------------------------------
# lollipop plot
#-------------------------------------------
png(file.path(output_dir, "CDH1_LollipopPlot.png"), width = 1000, height = 600, res = 150)

lollipopPlot(
  maf = combined_maf,
  gene = "CDH1",
  AACol = "HGVSp_Short",
  showMutationRate = TRUE
)

dev.off()

# Subset MAF data for CDH1
cdh1_mutations <- subset(combined_maf@data, Hugo_Symbol == "CDH1")

# Extract unique protein change positions from HGVSp_Short
# Example format: "p.E200K" → extract "200"
library(stringr)

# Use regex to extract numbers from the middle of strings like "p.E200K"
cdh1_positions <- unique(as.numeric(str_extract(cdh1_mutations$HGVSp_Short, "\\d+")))

# Remove NAs if any malformed entries exist
cdh1_positions <- cdh1_positions[!is.na(cdh1_positions)]

# Save PDF
pdf(file.path(output_dir, "CDH1_LollipopPlot_Labeled.pdf"), width = 10, height = 6)

lollipopPlot(
  maf = combined_maf,
  gene = "CDH1",
  AACol = "HGVSp_Short",
  showMutationRate = TRUE,
  labelPos = cdh1_positions
)

dev.off()

cat("Labeled positions:", sort(cdh1_positions), "\n")

#-------------------------------------------
# Compare mutation load against TCGA cohorts
#-------------------------------------------
laml.mutload = tcgaCompare(
  maf = laml,              # Your MAF object (should be 'combined_maf' in your case)
  cohortName = 'pILC',     # Your cohort name for the legend/title
  logscale = TRUE,         # Log10 scale for Y-axis (good for wide TMB range)
  capture_size = 50        # Size of the coding region (Mb), affects TMB calculation
)

png(file.path(output_dir, "MutationLoad_Comparison.png"), width = 1200, height = 800, res = 150)
tcgaCompare(
  maf = combined_maf,
  cohortName = 'pILC',
  logscale = TRUE,
  capture_size = 50
)
dev.off()

pdf(file.path(output_dir, "MutationLoad_Comparison.pdf"), width = 12, height = 8)
tcgaCompare(
  maf = combined_maf,
  cohortName = 'pILC',
  logscale = TRUE,
  capture_size = 50
)
dev.off()

#-------------------------------------------
# Compare mutation load against TCGA cohorts
#-------------------------------------------
# Create a new VAF column
combined_maf@data$CustomVAF <- combined_maf@data$t_alt_count / combined_maf@data$t_depth
plotVaf(maf = combined_maf, vafCol = "CustomVAF")

png(file.path(output_dir, "VAF_Distribution.png"), width = 1000, height = 800, res = 150)
plotVaf(maf = combined_maf, vafCol = "CustomVAF")
dev.off()

pdf(file.path(output_dir, "VAF_Distribution.pdf"), width = 10, height = 8)
plotVaf(maf = combined_maf, vafCol = "CustomVAF")
dev.off()

# Peaks near 0.5 → heterozygous clonal mutations
# Peaks near 1.0 → LOH or copy-neutral homozygous
# Broader spread or low VAFs (<0.3) → subclonal populations or noise

#-------------------------------------------
# Plot Somatic Interactions
#-------------------------------------------
somaticInteractions(
  maf = combined_maf,
  top = 50,
  pvalue = c(0.05, 0.1)  # Two thresholds shown on the plot
)
# PNG
png(file.path(output_dir, "Gene_Interaction_CoOccurrence.png"), width = 1200, height = 1000, res = 150)
somaticInteractions(maf = combined_maf, top = 25, pvalue = c(0.05, 0.1))
dev.off()

# PDF
pdf(file.path(output_dir, "Gene_Interaction_CoOccurrence.pdf"), width = 12, height = 10)
somaticInteractions(maf = combined_maf, top = 25, pvalue = c(0.05, 0.1))
dev.off()

#-------------------------------------------
# Drug-Gene Interactions
#-------------------------------------------
dgi <- drugInteractions(maf = combined_maf)
head(dgi)

png(
  file = file.path(output_dir, "DrugGeneInteractions.png"),
  width = 12,
  height = 8,
  units = "in",
  res = 150
)
drugInteractions(maf = combined_maf, fontSize = 1.0)
dev.off()

pdf(file.path(output_dir, "DrugGeneInteractions.pdf"), width = 12, height = 8)
drugInteractions(maf = combined_maf, fontSize = 1.0)
dev.off()

#-------------------------------------------
# Oncogenic Signaling Pathways
#-------------------------------------------
pws <- pathways(maf = combined_maf, plotType = "treemap")

png(file.path(output_dir, "Pathway_Treemap.png"), width = 1200, height = 1000, res = 150)
pathways(maf = combined_maf, plotType = "treemap")
dev.off()

pdf(file.path(output_dir, "Pathway_Treemap.pdf"), width = 12, height = 10)
pathways(maf = combined_maf, plotType = "treemap")
dev.off()

#-------------------------------------------
# Tumor heterogeneity and MATH scores
#-------------------------------------------

library(mclust)

sample_ids <- getSampleSummary(combined_maf)$Tumor_Sample_Barcode

for (s in sample_ids) {
  cat("Processing", s, "\n")
  result <- tryCatch({
    inferHeterogeneity(maf = combined_maf, tsb = s)
  }, error = function(e) {
    cat("Failed for", s, ":", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    png(
      filename = file.path(output_dir, paste0("Heterogeneity_", s, ".png")),
      width = 10,
      height = 8,
      units = "in",
      res = 150
    )
    
    plotClusters(clusters = result)
    
    dev.off()
  }
  
  if (!is.null(result)) {
    pdf(file.path(output_dir, paste0("Heterogeneity_", s, ".pdf")), width = 10, height = 8)
    plotClusters(clusters = result)
    dev.off()
  }
}

## Compute MATH Score per Sample
library(data.table)

# Extract the MAF data into a data.table
maf_dt <- as.data.table(combined_maf@data)

# Create VAF column
maf_dt[, t_vaf := t_alt_count / t_depth]

# Filter out NA or invalid values
maf_dt <- maf_dt[!is.na(t_vaf) & t_vaf > 0]

# Calculate MATH Score per sample
math_scores <- maf_dt[, .(
  MedianVAF = median(t_vaf, na.rm = TRUE),
  MAD_VAF = mad(t_vaf, na.rm = TRUE),
  MATH_Score = round((mad(t_vaf, na.rm = TRUE) / median(t_vaf, na.rm = TRUE)) * 100, 2)
), by = Tumor_Sample_Barcode]

print(math_scores)

fwrite(math_scores, file = file.path(output_dir, "MATH_scores.tsv"), sep = "\t")

#-------------------------------------------
# Plot Signature analysis
#-------------------------------------------
# Compute trinucleotide matrix
# Modify combined_maf@data in-place (Chromosome names in MAF must match chromosome names in reference genome)
combined_maf@data$Chromosome <- paste0("chr", combined_maf@data$Chromosome)

library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
combined_tnm <- trinucleotideMatrix(maf = combined_maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", useSyn = TRUE)

plotApobecDiff(tnm = combined_tnm, maf = combined_maf, pVal = 0.2)

library('NMF')
combined.sign = estimateSignatures(mat = combined_tnm, nTry = 6)

plotCophenetic(res = combined.sign)
png(
  filename = file.path(output_dir, "Signature_Elbow_Plot.png"),
  width = 1000,
  height = 800,
  res = 150
)
plotCophenetic(res = combined.sign)
dev.off()

combined.sig = extractSignatures(mat = combined_tnm, n = 2)

#Compate against original 30 signatures 
combined.og30.cosm = compareSignatures(nmfRes = combined.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures 
combined.og60.cosm = compareSignatures(nmfRes = combined.sig, sig_db = "SBS")

library(pheatmap)

# Save heatmap of cosine similarities
heatmap_file <- file.path(output_dir, "COSMIC_v3_SBS_Heatmap.png")
png(
  filename = heatmap_file,
  width = 1200,
  height = 1000,
  units = "px",
  res = 150
)
pheatmap::pheatmap(
  mat = combined.og60.cosm$cosine_similarities,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Cosine Similarity with COSMIC v3 SBS Signatures",
  fontsize_row = 10,
  fontsize_col = 8
)
dev.off()  # 

signature_file <- file.path(output_dir, "Combined_Mutational_Signatures_COSMICv3.png")

png(
  filename = signature_file,
  width = 1200,
  height = 900,
  res = 150
)

maftools::plotSignatures(
  nmfRes = combined.sig,
  title_size = 1.4,
  sig_db = "SBS"
)

dev.off() 

#-------------------------------------------
# Prepare MAF file for MutSigCV analysis
#-------------------------------------------
mutsig_ready_maf <- prepareMutSig(maf = combined_maf)

write.table(
  mutsig_ready_maf,
  file = file.path(output_dir, "MutSig_Input.maf"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==========================================================
# Oncoplot for Top 20 Mutated Genes (TCGA + ICGC Unique)
# ==========================================================
library(data.table)
library(maftools)
library(dplyr)

# ---- 1. File paths ----
tcga_file <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Reference/TCGA.BRCA.onco.genes.OncoVar.tsv"
icgc_file <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Reference/ICGC.BRCA.onco.genes.OncoVar.tsv"
output_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/PCGR/MAFtool"

# ---- 2. Read and combine gene lists ----
tcga_genes <- unique(fread(tcga_file)$Gene_symbol)
icgc_genes <- unique(fread(icgc_file)$Gene_symbol)
unique_genes <- union(tcga_genes, icgc_genes)

# ---- 3. Summarize mutations from MAF ----
maf_gene_summary <- getGeneSummary(combined_maf)
col_to_use <- ifelse("total" %in% colnames(maf_gene_summary), "total", "Total")

maf_filtered <- maf_gene_summary %>%
  filter(Hugo_Symbol %in% unique_genes) %>%
  arrange(desc(.data[[col_to_use]]))

# ---- 4. Select top 20 genes ----
top20_genes <- maf_filtered %>%
  slice_head(n = 20) %>%
  pull(Hugo_Symbol)

cat("Top 20 mutated unique TCGA+ICGC OncoVar genes:\n")
print(top20_genes)

# ---- 5. Plot high-quality oncoplot ----
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

png(file.path(output_dir, "Oncoplot_Top20_TCGA_ICGC_Unique.png"),
    width = 2600, height = 1800, res = 200)

oncoplot(
  maf = combined_maf,
  genes = top20_genes,
  removeNonMutated = TRUE,
  showTumorSampleBarcodes = FALSE,  # hide sample labels for cleaner look
  fontSize = 1.6,                   # enlarge gene name text
  annotationFontSize = 1.4,         # enlarge annotation labels
  gene_mar = 12,                    # increase spacing for gene labels
  draw_titv = FALSE,                # disable Ti/Tv summary panel
  top = 20
)
dev.off()

pdf(file.path(output_dir, "Oncoplot_Top20_TCGA_ICGC_Unique.pdf"),
    width = 18, height = 14)

oncoplot(
  maf = combined_maf,
  genes = top20_genes,
  removeNonMutated = TRUE,
  showTumorSampleBarcodes = FALSE,
  fontSize = 1.6,
  annotationFontSize = 1.4,
  gene_mar = 12,
  draw_titv = FALSE,
  top = 20
)
dev.off()

#-------------------------------------------
# Plot common Breast Cancer Genes
#-------------------------------------------

# Define breast cancer driver genes of interest
brca_genes <- c(
  "BRCA1", "BRCA2", "PALB2", "TP53", "PTEN", "CDH1", "ATM", "CHEK2",
  "BARD1", "BRIP1", "NF1", "RAD51C", "RAD51D", "STK11", "EPCAM", "MSH2",
  "MLH1", "MSH6", "PMS2", "ESR1", "CYP19A1", "ERBB2", "ERBB3"
)

# Draw oncoplot
oncoplot(
  maf = combined_maf,
  genes = brca_genes,
  removeNonMutated = TRUE,
  draw_titv = FALSE,
  sortByAnnotation = FALSE,
  titleText = "Breast Cancer-Associated Mutations"
)
png(file.path(output_dir, "Breast_Cancer_Genes_Oncoplot.png"), width = 1400, height = 1000, res = 150)

oncoplot(
  maf = combined_maf,
  genes = brca_genes,
  removeNonMutated = TRUE,
  draw_titv = FALSE,
  sortByAnnotation = FALSE,
  titleText = "Breast Cancer-Associated Mutations"
)

dev.off()
