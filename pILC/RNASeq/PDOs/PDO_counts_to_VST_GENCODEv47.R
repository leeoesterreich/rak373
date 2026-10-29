############################################################
# PDO raw counts -> gene-level raw counts and VST matrices
# 4 PDO samples; gene labels come from the GTF used for counting
# (GENCODE v47 basic, GRCh38) - no biomaRt, no version stripping
#
# Input : RNASeq/PDO/PDO_raw_counts.csv  (column "gene" + one column per sample, RAW counts)
# Output: RNASeq/PDO/VST_data_GENCODEv47/
############################################################

set.seed(42)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

############################################################
## 1. Paths
############################################################

pdo_dir    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/PDO"
count_file <- file.path(pdo_dir, "PDO_raw_counts.csv")

genome_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7"
gtf_file   <- file.path(genome_dir, "gencode.v47.basic.annotation.gtf")

out_dir <- file.path(pdo_dir, "VST_data_GENCODEv47")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(gtf_file))   stop("GTF not found: ", gtf_file)
if (!file.exists(count_file)) stop("Count file not found: ", count_file)

############################################################
## 2. Read gene annotation from the GTF (gene lines only)
############################################################

gtf <- fread(
  cmd = paste("grep -v '^#'", shQuote(gtf_file)),
  sep = "\t", header = FALSE, quote = "",
  select = c(3, 9), col.names = c("feature", "attributes")
)
gtf <- gtf[feature == "gene"]

get_attr <- function(x, key) {
  pattern <- paste0('.*', key, ' "([^"]+)".*')
  ifelse(grepl(paste0(key, ' "'), x), sub(pattern, "\\1", x), NA_character_)
}

gene_annot <- data.table(
  ensembl_id   = get_attr(gtf$attributes, "gene_id"),
  gene_symbol  = get_attr(gtf$attributes, "gene_name"),
  gene_biotype = get_attr(gtf$attributes, "gene_type")  # GENCODE uses gene_type
)
gene_annot <- unique(gene_annot, by = "ensembl_id")

cat("Genes in GTF:", nrow(gene_annot), "\n\n")

############################################################
## 3. Read the PDO count matrix (strict checks, no silent NA -> 0)
############################################################

raw <- fread(count_file, header = TRUE, check.names = FALSE)
ids <- as.character(raw[[1]])                    # first column "gene": Ensembl ID with version
counts_mat <- as.matrix(raw[, -1])
rownames(counts_mat) <- ids

cat("Samples:", paste(colnames(counts_mat), collapse = ", "), "\n")
if (ncol(counts_mat) != 4) stop("Expected 4 PDO samples, found ", ncol(counts_mat))

# Remove HTSeq summary rows if present (__no_feature, __ambiguous, ...)
counts_mat <- counts_mat[!grepl("^__", rownames(counts_mat)), , drop = FALSE]

if (anyNA(counts_mat))                   stop("Missing count values detected.")
if (anyDuplicated(rownames(counts_mat))) stop("Duplicated gene IDs in the count file.")
if (any(counts_mat != round(counts_mat))) stop("Counts are not integers.")
cat("Count matrix: ", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n\n")

############################################################
## 4. Build unique gene labels from the GTF
##    Exact ID match first (versioned IDs); if a gene is not
##    found, fall back to the ID without the version suffix.
############################################################

strip_ver <- function(x) sub("\\.\\d+$", "", x)

idx <- match(rownames(counts_mat), gene_annot$ensembl_id)
n_exact_miss <- sum(is.na(idx))
if (n_exact_miss > 0) {
  fallback <- match(strip_ver(rownames(counts_mat)), strip_ver(gene_annot$ensembl_id))
  idx[is.na(idx)] <- fallback[is.na(idx)]
}
cat("Count IDs not matched exactly:", n_exact_miss, "\n")
cat("Count IDs not found in GTF   :", sum(is.na(idx)), "\n")

symbol <- gene_annot$gene_symbol[idx]
symbol <- ifelse(is.na(symbol) | symbol == "", rownames(counts_mat), symbol)

# Where a symbol is shared by several ENSG IDs, use SYMBOL|ENSG; otherwise keep SYMBOL
dup_symbol <- duplicated(symbol) | duplicated(symbol, fromLast = TRUE)
gene_label <- ifelse(dup_symbol, paste(symbol, rownames(counts_mat), sep = "|"), symbol)
stopifnot(!anyDuplicated(gene_label))

annot_out <- data.table(
  Gene         = gene_label,
  ensembl_id   = rownames(counts_mat),
  gene_symbol  = symbol,
  gene_biotype = gene_annot$gene_biotype[idx]
)
fwrite(annot_out, file.path(out_dir, "gene_annotation_from_GTF.tsv"), sep = "\t")

rownames(counts_mat) <- gene_label

cat("Genes still labelled by ENSG ID only:", sum(grepl("^ENSG", gene_label) & !grepl("\\|", gene_label)), "\n\n")

############################################################
## 5. Save raw counts (input for the PCA / DESeq2 scripts)
############################################################

raw_out <- data.table(Gene = rownames(counts_mat), counts_mat)
fwrite(raw_out, file.path(out_dir, "htseq_raw_counts_with_gene_names.csv"))
saveRDS(counts_mat, file.path(out_dir, "htseq_raw_counts_with_gene_names.rds"))

############################################################
## 6. Low-expression filter before VST
##    Same rule as the tumors: count >= 5 in at least 10% of samples
##    (with 4 samples: ceiling(0.10 * 4) = 1 sample)
############################################################

min_count   <- 5
min_samples <- ceiling(ncol(counts_mat) * 0.10)

keep_genes      <- rowSums(counts_mat >= min_count) >= min_samples
counts_filtered <- counts_mat[keep_genes, , drop = FALSE]

cat("Genes after low-expression filter:", nrow(counts_filtered), "\n\n")

fwrite(
  data.table(Gene = rownames(counts_filtered), counts_filtered),
  file.path(out_dir, "htseq_expression_filtered_counts_gene_names.csv")
)
saveRDS(counts_filtered, file.path(out_dir, "htseq_expression_filtered_counts_gene_names.rds"))

############################################################
## 7. VST (blind = TRUE, design ~ 1: unsupervised use)
############################################################

storage.mode(counts_filtered) <- "integer"

coldata <- data.frame(sample = colnames(counts_filtered),
                      row.names = colnames(counts_filtered))

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = coldata,
  design    = ~ 1
)

# vst() needs >= 1000 genes (fine here); fall back as in the tumor script
vst_obj <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("vst() failed: ", conditionMessage(e),
            "\nUsing varianceStabilizingTransformation() instead.")
    varianceStabilizingTransformation(dds, blind = TRUE)
  }
)
vst_mat <- assay(vst_obj)

fwrite(
  data.table(Gene = rownames(vst_mat), vst_mat),
  file.path(out_dir, "htseq_VST_gene_names_all_expression_filtered_genes.csv")
)
saveRDS(vst_obj, file.path(out_dir, "DESeq2_VST_object.rds"))
saveRDS(vst_mat, file.path(out_dir, "htseq_VST_matrix_all_expression_filtered_genes.rds"))

############################################################
## 8. Top 5000 variable genes after VST
############################################################

top_n_genes   <- 5000
vst_gene_vars <- apply(vst_mat, 1, var, na.rm = TRUE)
top_genes     <- names(sort(vst_gene_vars, decreasing = TRUE))[1:min(top_n_genes, length(vst_gene_vars))]
vst_mat_topvar <- vst_mat[top_genes, , drop = FALSE]

fwrite(
  data.table(Gene = rownames(vst_mat_topvar), vst_mat_topvar),
  file.path(out_dir, "htseq_VST_gene_names_top5000_variable_genes.csv")
)
saveRDS(vst_mat_topvar, file.path(out_dir, "htseq_VST_matrix_top5000_variable_genes.rds"))

gene_variance_df <- data.table(Gene = names(vst_gene_vars),
                               VST_variance = as.numeric(vst_gene_vars))
setorder(gene_variance_df, -VST_variance)
fwrite(gene_variance_df, file.path(out_dir, "htseq_VST_gene_variance_table.csv"))

fwrite(data.table(Sample = colnames(vst_mat)),
       file.path(out_dir, "VST_sample_names.csv"))

############################################################
## 9. GCT file (same format as the tumor file, all expression-filtered genes)
############################################################

gct_file <- file.path(out_dir, "vst_transformed_data.gct")

desc <- annot_out$ensembl_id[match(rownames(vst_mat), annot_out$Gene)]  # Description = Ensembl ID

con <- file(gct_file, "w")
writeLines("#1.2", con)
writeLines(paste(nrow(vst_mat), ncol(vst_mat), sep = "\t"), con)
close(con)

gct_body <- data.table(Name = rownames(vst_mat), Description = desc, vst_mat)
fwrite(gct_body, gct_file, sep = "\t", append = TRUE, col.names = TRUE)

############################################################
## 10. Quick check: assigned reads per sample
############################################################

lib_check <- data.table(Sample = colnames(counts_mat),
                        assigned_reads = colSums(counts_mat))
print(lib_check)
fwrite(lib_check, file.path(out_dir, "PDO_assigned_reads_per_sample.csv"))
stopifnot(!anyNA(vst_mat), all(is.finite(vst_mat)), min(vst_mat) >= 0)

############################################################
## 11. Final summary
############################################################

cat("\n============================================================\n")
cat("DONE: PDO VST RNA-seq matrix generated\n")
cat("============================================================\n")
cat("Output folder:", out_dir, "\n\n")
cat("Raw counts:", file.path(out_dir, "htseq_raw_counts_with_gene_names.csv"), "\n\n")
cat("Genes after removing summary rows :", nrow(counts_mat), "\n")
cat("Genes after low-expression filter :", nrow(counts_filtered), "\n")
cat("Genes in top-variable VST matrix  :", nrow(vst_mat_topvar), "\n")
cat("Number of samples                 :", ncol(vst_mat), "\n")
