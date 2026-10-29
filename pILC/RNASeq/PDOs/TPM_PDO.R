# ============================================================
# TPM matrix for the 4 PDO RNA-seq samples (same method as TPM_pILC_stranded.R)
# Gene length = union of exons (GENCODE v47 basic) -> TPM -> HGNC symbol from the GTF,
# one Ensembl ID per HGNC symbol (highest max TPM).
# Input : RNASeq/PDO/PDO_raw_counts.csv  (column "gene" + one column per sample, RAW counts)
# Output: RNASeq/PDO/Max_TPM/
# ============================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(txdbmaker)    # makeTxDbFromGFF() lives here (Bioconductor >= 3.19)
  library(rtracklayer)
  library(dplyr)
})

# ---------------- 1. Paths ----------------
pdo_dir    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/PDO"
count_file <- file.path(pdo_dir, "PDO_raw_counts.csv")
gtf_file   <- "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/gencode.v47.basic.annotation.gtf"

tpm_dir              <- file.path(pdo_dir, "Max_TPM")
output_tpm_file      <- file.path(tpm_dir, "tpm_matrix_PDO.txt")
output_log2_tpm_file <- file.path(tpm_dir, "tpm_matrix_PDO_log2TPM.txt")
output_map_file      <- file.path(tpm_dir, "HGNC_to_Ensembl_mapping.txt")
dir.create(tpm_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- 2. Read the single count matrix ----------------
if (!file.exists(count_file)) stop("Count file not found: ", count_file)
raw <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- as.character(raw[[1]])                       # first column = Ensembl ID (with version)
merged   <- as.matrix(raw[, -1, drop = FALSE])
mode(merged) <- "numeric"
rownames(merged) <- gene_ids
message("Samples found: ", paste(colnames(merged), collapse = ", "))
if (ncol(merged) != 4) stop("Expected 4 PDO samples, found ", ncol(merged))

# Drop HTSeq-style summary rows if present, check values, make integers
merged <- merged[!grepl("^__", rownames(merged)), , drop = FALSE]
if (anyNA(merged)) stop("Missing values in the count matrix")
merged <- round(merged)
rownames(merged) <- gsub("\\.\\d+$", "", rownames(merged))   # strip Ensembl version
if (anyDuplicated(rownames(merged))) stop("Duplicated Ensembl IDs after stripping versions")
message("Count matrix: ", nrow(merged), " genes x ", ncol(merged), " samples")

# ---------------- 3. Gene lengths from the GTF (non-overlapping exon length per gene) ----------------
txdb  <- makeTxDbFromGFF(gtf_file, format = "gtf")
exons <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(exons)))
length_df <- data.frame(ensembl_id     = gsub("\\.\\d+$", "", names(gene_lengths)),
                        gene_length_kb = as.numeric(gene_lengths) / 1000)

# ---------------- 4. TPM = (counts / length), scaled so each sample sums to 1e6 ----------------
merged_df <- data.frame(ensembl_id = rownames(merged), merged, check.names = FALSE)
tpm_input <- merge(merged_df, length_df, by = "ensembl_id")
message("Genes with a GTF length: ", nrow(tpm_input), " of ", nrow(merged_df))

raw_counts <- tpm_input[, !(names(tpm_input) %in% c("ensembl_id", "gene_length_kb")), drop = FALSE]
rpk <- sweep(raw_counts, 1, tpm_input$gene_length_kb, "/")
tpm <- as.data.frame(sweep(rpk, 2, colSums(rpk), "/") * 1e6, check.names = FALSE)
tpm$ensembl_id <- tpm_input$ensembl_id

# ---------------- 5. Ensembl -> HGNC symbol (from the same GENCODE v47 GTF) ----------------
message("Reading gene symbols from the GTF...")
gtf_genes <- rtracklayer::import(gtf_file, feature.type = "gene")
id_map <- data.frame(ensembl_gene_id = gsub("\\.\\d+$", "", gtf_genes$gene_id),
                     hgnc_symbol     = gtf_genes$gene_name,
                     stringsAsFactors = FALSE)
id_map <- id_map[!is.na(id_map$hgnc_symbol), ]
id_map <- id_map[!duplicated(id_map$ensembl_gene_id), ]
id_map <- id_map[id_map$hgnc_symbol != "", ]
message("Genes with a symbol: ", nrow(id_map))

# ---------------- 6. One Ensembl ID per HGNC symbol: highest max TPM ----------------
tpm <- merge(tpm, id_map, by.x = "ensembl_id", by.y = "ensembl_gene_id")
sample_cols <- setdiff(names(tpm), c("ensembl_id", "hgnc_symbol"))
tpm[sample_cols] <- lapply(tpm[sample_cols], function(x) as.numeric(as.character(x)))
tpm$score_max <- apply(tpm[, sample_cols, drop = FALSE], 1, max, na.rm = TRUE)

# explicit dplyr:: because GenomicFeatures/rtracklayer mask select(), rename(), slice()
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

# ---------------- 7. Save TPM, log2(TPM + 1) and per-sample TPM files ----------------
write.table(final_matrix, output_tpm_file, sep = "\t", quote = FALSE, row.names = FALSE)

log2_tpm <- final_matrix
log2_tpm[, -1] <- log2(log2_tpm[, -1] + 1)
write.table(log2_tpm, output_log2_tpm_file, sep = "\t", quote = FALSE, row.names = FALSE)

for (s in colnames(final_matrix)[-1]) {
  write.table(final_matrix[, c("Gene", s)],
              file.path(tpm_dir, paste0(s, "_TPM.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Sanity check: every sample's TPM should sum to ~1e6 before the one-ID-per-symbol step
message("TPM column sums before symbol collapsing (expect 1e6): ",
        paste(round(colSums(tpm[sample_cols])), collapse = ", "))
message("TPM saved: ", output_tpm_file, "  (", nrow(final_matrix), " genes x ",
        ncol(final_matrix) - 1, " samples)")
