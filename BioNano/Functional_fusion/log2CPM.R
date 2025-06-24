##### This script was used to convert HTSeq to log2CPM for further downstreaming analysis

library(edgeR)
library(biomaRt)  # For gene annotation

# Paths
input_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/4_Count_file"
output_csv <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/4_Count_file/pILC_log2cpm_matrix.csv"
output_rdata <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/4_Count_file/pILC_log2CPM_matrix.Rdata"

# Step 1: Read all HTSeq files
count_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
count_list <- lapply(count_files, function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("Gene", basename(file))
  df
})

# Step 2: Merge all count files
merged <- Reduce(function(x, y) merge(x, y, by = "Gene"), count_list)
rownames(merged) <- merged$Gene
merged <- merged[, -1]

# Step 3: Remove HTSeq special rows
exclude <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
merged <- merged[!rownames(merged) %in% exclude, ]

# Step 4: Convert to integer matrix
gene_ids <- rownames(merged)
merged <- apply(merged, 2, as.integer)
rownames(merged) <- gene_ids

# Step 5: Clean column names
colnames(merged) <- gsub("\\.txt$", "", colnames(merged))

# Step 6: Calculate log2CPM
dge <- DGEList(counts = merged)
log2_cpm <- cpm(dge, log = TRUE, prior.count = 1)
dim(log2_cpm)
head(log2_cpm)

# Step 7: Annotate Ensembl IDs with gene names
ensembl_ids <- gsub("\\..*", "", rownames(log2_cpm))

# Use biomaRt to fetch gene symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  
gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "ensembl_gene_id",
                  values = ensembl_ids,
                  mart = mart)
dim(gene_map)
head(gene_map)

# Merge annotation
log2_cpm_df <- as.data.frame(log2_cpm)
log2_cpm_df$Ensembl_ID <- gsub("\\..*", "", rownames(log2_cpm_df))
annotated_log2cpm <- merge(gene_map, log2_cpm_df, by.x = "ensembl_gene_id", by.y = "Ensembl_ID")

# Reorder: Ensembl ID, Gene Name, then samples
final_log2cpm <- annotated_log2cpm[, c("ensembl_gene_id", "external_gene_name", setdiff(names(annotated_log2cpm), c("ensembl_gene_id", "external_gene_name")))]

# Step 8: Save output
write.csv(final_log2cpm, file = output_csv, row.names = FALSE)
save(final_log2cpm, file = output_rdata)

cat("Saved log2CPM with Ensembl IDs and gene names.\n")

#################################################################################################
###### Convert .csv into .Rdata
# Set the input and output file paths
input_path <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/4_Count_file/pILC_log2cpm_matrix.csv"
output_path <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/4_Count_file/pILC_log2CPM_matrix.Rdata"

# Read the CSV file
SV_df <- read.csv(input_path, stringsAsFactors = FALSE)

# Save it as .Rdata
save(SV_df, file = output_path)

# Optional: Confirm it's saved
cat("Saved as:", output_path, "\n")
