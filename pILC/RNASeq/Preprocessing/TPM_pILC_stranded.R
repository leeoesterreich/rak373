# ============================================================
# TPM matrix for pILC RNA-seq from the corrected (htseq-count -s reverse) counts
# Method (same as the earlier TPM): gene length = union of exons (GENCODE v47 basic),
# TPM, Ensembl -> HGNC via biomaRt, one Ensembl ID per HGNC (highest max TPM).
# Output: RNASeq/4_Count_file_HTSeq_stranded/Max_TPM/
# ============================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(txdbmaker)    # makeTxDbFromGFF() moved here from GenomicFeatures (Bioconductor >= 3.19)
  library(biomaRt)
  library(dplyr)
})

# ---------------- 1. Paths ----------------
input_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded"  # per-sample HTSeq counts
gtf_file  <- "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/gencode.v47.basic.annotation.gtf"

tpm_dir              <- file.path(input_dir, "Max_TPM")
output_tpm_file      <- file.path(tpm_dir, "tpm_matrix_pILC.txt")
output_log2_tpm_file <- file.path(tpm_dir, "tpm_matrix_pILC_log2TPM.txt")
output_map_file      <- file.path(tpm_dir, "HGNC_to_Ensembl_mapping.txt")
dir.create(tpm_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- 2. Read all HTSeq files (top level of the folder only) ----------------
count_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
message("Count files found: ", length(count_files))
print(basename(count_files))   # CHECK: should be the 31 per-sample files only
if (length(count_files) != 31)
  stop("Expected 31 per-sample count files, found ", length(count_files),
       ". Check the list above (extra or missing .txt files?).")

count_list <- lapply(count_files, function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("Gene", basename(file))
  df
})

# ---------------- 3. Merge, drop HTSeq summary rows, integer matrix ----------------
merged <- Reduce(function(x, y) merge(x, y, by = "Gene"), count_list)
rownames(merged) <- merged$Gene
merged <- merged[, -1]

exclude <- c("__no_feature", "__ambiguous", "__too_low_aQual",
             "__not_aligned", "__alignment_not_unique")
merged <- merged[!rownames(merged) %in% exclude, ]

gene_ids <- rownames(merged)
merged <- apply(merged, 2, as.integer)
rownames(merged) <- gene_ids
colnames(merged) <- gsub("\\.txt$", "", colnames(merged))
rownames(merged) <- gsub("\\.\\d+$", "", rownames(merged))   # strip Ensembl version

# ---------------- 4. Gene lengths from the GTF (non-overlapping exon length per gene) ----------------
txdb  <- makeTxDbFromGFF(gtf_file, format = "gtf")
exons <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(exons)))
length_df <- data.frame(ensembl_id     = gsub("\\.\\d+$", "", names(gene_lengths)),
                        gene_length_kb = as.numeric(gene_lengths) / 1000)

# ---------------- 5. TPM = (counts / length), scaled so each sample sums to 1e6 ----------------
merged_df <- data.frame(ensembl_id = rownames(merged), merged)
tpm_input <- merge(merged_df, length_df, by = "ensembl_id")

raw_counts <- tpm_input[, !(names(tpm_input) %in% c("ensembl_id", "gene_length_kb"))]
rpk <- sweep(raw_counts, 1, tpm_input$gene_length_kb, "/")
tpm <- as.data.frame(sweep(rpk, 2, colSums(rpk), "/") * 1e6)
tpm$ensembl_id <- tpm_input$ensembl_id

# ---------------- 6. Ensembl -> HGNC (same biomaRt version as the earlier TPM) ----------------
USE_BIOMART <- FALSE   # FALSE = symbols from the GTF (no internet needed); TRUE = biomaRt 109 as before

if (USE_BIOMART) {
  message("Querying biomaRt for HGNC symbols...")
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 109)
  id_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id", values = tpm$ensembl_id, mart = mart)
} else {
  # Gene symbols (gene_name) come straight from the same GENCODE v47 GTF used for the counts
  message("Reading gene symbols from the GTF...")
  gtf_genes <- rtracklayer::import(gtf_file, feature.type = "gene")
  id_map <- data.frame(ensembl_gene_id = gsub("\\.\\d+$", "", gtf_genes$gene_id),
                       hgnc_symbol     = gtf_genes$gene_name,
                       stringsAsFactors = FALSE)
  id_map <- id_map[!is.na(id_map$hgnc_symbol), ]
}
id_map <- id_map[!duplicated(id_map$ensembl_gene_id), ]
id_map <- id_map[id_map$hgnc_symbol != "", ]
message("Genes with a symbol: ", nrow(id_map))

# ---------------- 7. One Ensembl ID per HGNC symbol: highest max TPM ----------------
tpm <- merge(tpm, id_map, by.x = "ensembl_id", by.y = "ensembl_gene_id")
sample_cols <- setdiff(names(tpm), c("ensembl_id", "hgnc_symbol"))
tpm[sample_cols] <- lapply(tpm[sample_cols], function(x) as.numeric(as.character(x)))
tpm$score_max <- apply(tpm[, sample_cols, drop = FALSE], 1, max, na.rm = TRUE)

# explicit dplyr:: because GenomicFeatures/biomaRt mask select(), rename(), slice()
tpm_best <- tpm %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::arrange(dplyr::desc(score_max), ensembl_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

write.table(dplyr::select(tpm_best, hgnc_symbol, ensembl_id), output_map_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

final_matrix <- tpm_best %>%
  dplyr::select(hgnc_symbol, dplyr::all_of(sample_cols)) %>%
  dplyr::rename(Gene = hgnc_symbol)
colnames(final_matrix)[-1] <- gsub("\\.", "-", colnames(final_matrix)[-1])

# ---------------- 8. Save TPM, log2(TPM + 1) and per-sample TPM files ----------------
write.table(final_matrix, output_tpm_file, sep = "\t", quote = FALSE, row.names = FALSE)

log2_tpm <- final_matrix
log2_tpm[, -1] <- log2(log2_tpm[, -1] + 1)
write.table(log2_tpm, output_log2_tpm_file, sep = "\t", quote = FALSE, row.names = FALSE)

for (s in colnames(final_matrix)[-1]) {
  write.table(final_matrix[, c("Gene", s)],
              file.path(tpm_dir, paste0(s, "_TPM.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("TPM saved: ", output_tpm_file, "  (", nrow(final_matrix), " genes x ",
        ncol(final_matrix) - 1, " samples)")
