############################################################
## ILC_HTSeq_counts_to_VST_GENCODEv47.R
##
## PURPOSE (no batch correction)
##   HTSeq count files of the two ILC sets, processed SEPARATELY:
##     Chelsea  21 samples, FFPE          -> own filter + VST + top 5000 + GCT
##     Bionano   3 samples, fresh frozen  -> raw counts only (a VST on 3 samples
##                                           alone is not meaningful; use them
##                                           together with pILC, see DESeq2 script)
##   GENCODE v47 gene labels, filter rule as in the pILC script.
##
## OUTPUT  2_Exploratory_strandedness/<setting>/
##   (root)            sample_manifest.csv, gene_annotation_from_GTF.tsv,
##                     htseq_assignment_QC_per_sample.csv,
##                     htseq_raw_counts_with_gene_names.csv / .rds  (all 24 samples,
##                     input for DESeq2_ILC_vs_pILC.R and the diagnostics script)
##   Chelsea_FFPE_21/  raw, filtered, VST (all genes), top 5000, GCT, variance table
##   Bionano_FF_3/     raw counts (csv + rds)
##
## RUN
##   Rscript ILC_HTSeq_counts_to_VST_GENCODEv47.R Reverse     (or No)
##   In RStudio the default below is used.
############################################################

SEED <- 42
set.seed(SEED)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

############################################################
## 0. Setting to run
############################################################

args <- commandArgs(trailingOnly = TRUE)
strand_setting <- if (length(args) >= 1) args[1] else "Reverse"   # "Reverse", "No" or "Yes"
stopifnot(strand_setting %in% c("Reverse", "No", "Yes"))

min_count <- 5      # same rule as the pILC script:
min_frac  <- 0.10   # count >= 5 in at least 10% of the samples

############################################################
## 1. Paths and explicit sample list
############################################################

base_dir    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_ILC"
bionano_dir <- file.path(base_dir, "1_Counts", "Bionano", paste0("Stranded_", strand_setting))
chelsea_dir <- file.path(base_dir, "1_Counts", "Chelsea", paste0("Stranded_", strand_setting))

genome_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7"
gtf_file   <- file.path(genome_dir, "gencode.v47.basic.annotation.gtf")
if (!file.exists(gtf_file)) stop("GTF not found: ", gtf_file)

out_dir      <- file.path(base_dir, "2_Exploratory_strandedness", strand_setting)
out_chelsea  <- file.path(out_dir, "Chelsea_FFPE_21")
out_bionano  <- file.path(out_dir, "Bionano_FF_3")
for (d in c(out_dir, out_chelsea, out_bionano)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("Strandedness setting:", strand_setting, "\n")
cat("Bionano folder:", bionano_dir, "\n")
cat("Chelsea folder:", chelsea_dir, "\n")
cat("Output folder :", out_dir, "\n\n")

chelsea_ids <- c("TP17_M378", "TP17_M383", "TP17_M384", "TP17_M510",
                 "TP17_M559", "TP17_M605", "TP17_M622", "TP17_M802",
                 "TP17_M882", "TP17_M899", "TP17_M901",
                 "TP18_M95",  "TP18_M121", "TP18_M202", "TP18_M251",
                 "TP18_M278", "TP18_M302", "TP18_M329", "TP18_M355",
                 "TP18_M372", "TP18_M403")
bionano_ids <- c("TP10-1655", "TP12-813", "TP14-79")

manifest <- rbind(
  data.table(file = file.path(chelsea_dir, paste0(chelsea_ids, ".txt")),
             raw_id = chelsea_ids, count_folder = "Chelsea"),
  data.table(file = file.path(bionano_dir, paste0(bionano_ids, ".txt")),
             raw_id = bionano_ids, count_folder = "Bionano")
)
manifest[, strand_setting := strand_setting]
manifest[, preservation := fifelse(count_folder == "Chelsea", "FFPE", "Fresh frozen")]
manifest[, sample := make.names(raw_id)]     # TP10-1655 -> TP10.1655 (same as the pILC script)

missing_files <- manifest$file[!file.exists(manifest$file)]
if (length(missing_files) > 0) stop("Missing count files:\n", paste(missing_files, collapse = "\n"))
if (anyDuplicated(manifest$sample)) stop("Duplicated sample names in the sample list.")
if (nrow(manifest) != 24) stop("Expected 24 samples, found ", nrow(manifest))
fwrite(manifest, file.path(out_dir, "sample_manifest.csv"))

############################################################
## 2. Gene annotation from the GTF (gene lines only)
############################################################

gtf <- fread(cmd = paste("grep -v '^#'", shQuote(gtf_file)),
             sep = "\t", header = FALSE, quote = "",
             select = c(3, 9), col.names = c("feature", "attributes"))
gtf <- gtf[feature == "gene"]

get_attr <- function(x, key) {
  pattern <- paste0('.*', key, ' "([^"]+)".*')
  ifelse(grepl(paste0(key, ' "'), x), sub(pattern, "\\1", x), NA_character_)
}
gene_annot <- data.table(ensembl_id   = get_attr(gtf$attributes, "gene_id"),
                         gene_symbol  = get_attr(gtf$attributes, "gene_name"),
                         gene_biotype = get_attr(gtf$attributes, "gene_type"))
gene_annot <- unique(gene_annot, by = "ensembl_id")
cat("Genes in GTF:", nrow(gene_annot), "\n\n")

############################################################
## 3. Read HTSeq files (strict checks, no silent NA -> 0)
############################################################

read_htseq <- function(file) fread(file, header = FALSE, sep = "\t", col.names = c("ensembl_id", "count"))
count_list <- lapply(manifest$file, read_htseq)
names(count_list) <- manifest$sample

ref_ids  <- count_list[[1]]$ensembl_id
same_ids <- vapply(count_list, function(x) identical(x$ensembl_id, ref_ids), logical(1))
if (!all(same_ids)) stop("Gene IDs differ between files: ", paste(names(same_ids)[!same_ids], collapse = ", "))

counts_all <- do.call(cbind, lapply(count_list, function(x) x$count))
rownames(counts_all) <- ref_ids
colnames(counts_all) <- manifest$sample

# HTSeq summary rows (__no_feature, __ambiguous, ...) are kept for QC, then removed
is_summary   <- grepl("^__", rownames(counts_all))
summary_rows <- counts_all[is_summary, , drop = FALSE]
counts_mat   <- counts_all[!is_summary, , drop = FALSE]

assign_qc <- data.table(Sample = colnames(counts_all), cohort = manifest$count_folder,
                        assigned_reads = colSums(counts_mat), total_reads = colSums(counts_all),
                        assigned_fraction = round(colSums(counts_mat) / colSums(counts_all), 4))
assign_qc <- cbind(assign_qc, t(summary_rows)[assign_qc$Sample, , drop = FALSE])
print(assign_qc)
fwrite(assign_qc, file.path(out_dir, "htseq_assignment_QC_per_sample.csv"))
if (anyNA(counts_mat)) stop("Missing count values detected.")

############################################################
## 4. Unique gene labels from the GTF
############################################################

idx <- match(rownames(counts_mat), gene_annot$ensembl_id)
cat("Count IDs not found in GTF:", sum(is.na(idx)), "\n")

symbol <- gene_annot$gene_symbol[idx]
symbol <- ifelse(is.na(symbol) | symbol == "", rownames(counts_mat), symbol)
dup_symbol <- duplicated(symbol) | duplicated(symbol, fromLast = TRUE)
gene_label <- ifelse(dup_symbol, paste(symbol, rownames(counts_mat), sep = "|"), symbol)
stopifnot(!anyDuplicated(gene_label))

annot_out <- data.table(Gene = gene_label, ensembl_id = rownames(counts_mat),
                        gene_symbol = symbol, gene_biotype = gene_annot$gene_biotype[idx])
fwrite(annot_out, file.path(out_dir, "gene_annotation_from_GTF.tsv"), sep = "\t")
rownames(counts_mat) <- gene_label

############################################################
## 5. Raw counts: all 24 samples (no filter, no correction)
############################################################

fwrite(data.table(Gene = rownames(counts_mat), counts_mat),
       file.path(out_dir, "htseq_raw_counts_with_gene_names.csv"))
saveRDS(counts_mat, file.path(out_dir, "htseq_raw_counts_with_gene_names.rds"))

# Bionano (3): raw counts only
bio_cols <- manifest$sample[manifest$count_folder == "Bionano"]
fwrite(data.table(Gene = rownames(counts_mat), counts_mat[, bio_cols, drop = FALSE]),
       file.path(out_bionano, "htseq_raw_counts_with_gene_names.csv"))
saveRDS(counts_mat[, bio_cols, drop = FALSE], file.path(out_bionano, "htseq_raw_counts_with_gene_names.rds"))

############################################################
## 6. Chelsea (21, FFPE): filter, VST, top 5000, GCT
############################################################

che_cols   <- manifest$sample[manifest$count_folder == "Chelsea"]
counts_che <- counts_mat[, che_cols, drop = FALSE]

fwrite(data.table(Gene = rownames(counts_che), counts_che),
       file.path(out_chelsea, "htseq_raw_counts_with_gene_names.csv"))
saveRDS(counts_che, file.path(out_chelsea, "htseq_raw_counts_with_gene_names.rds"))

min_samples <- ceiling(ncol(counts_che) * min_frac)          # 3 of 21
keep_genes  <- rowSums(counts_che >= min_count) >= min_samples
counts_filtered <- counts_che[keep_genes, , drop = FALSE]
storage.mode(counts_filtered) <- "integer"
cat("\nChelsea: genes after filter (count >=", min_count, "in >=", min_samples, "samples):",
    nrow(counts_filtered), "\n")

fwrite(data.table(Gene = rownames(counts_filtered), counts_filtered),
       file.path(out_chelsea, "htseq_expression_filtered_counts_gene_names.csv"))
saveRDS(counts_filtered, file.path(out_chelsea, "htseq_expression_filtered_counts_gene_names.rds"))

# VST (blind = TRUE, design ~ 1)
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = data.frame(sample = colnames(counts_filtered),
                                                   row.names = colnames(counts_filtered)),
                              design = ~ 1)
vst_obj <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("vst() failed: ", conditionMessage(e), "\nUsing varianceStabilizingTransformation() instead.")
    varianceStabilizingTransformation(dds, blind = TRUE)
  }
)
vst_mat <- assay(vst_obj)

fwrite(data.table(Gene = rownames(vst_mat), vst_mat),
       file.path(out_chelsea, "htseq_VST_gene_names_all_expression_filtered_genes.csv"))
saveRDS(vst_obj, file.path(out_chelsea, "DESeq2_VST_object.rds"))
saveRDS(vst_mat, file.path(out_chelsea, "htseq_VST_matrix_all_expression_filtered_genes.rds"))

# top 5000 variable genes
vst_gene_vars  <- apply(vst_mat, 1, var, na.rm = TRUE)
top_genes      <- names(sort(vst_gene_vars, decreasing = TRUE))[1:min(5000, length(vst_gene_vars))]
vst_mat_topvar <- vst_mat[top_genes, , drop = FALSE]
fwrite(data.table(Gene = rownames(vst_mat_topvar), vst_mat_topvar),
       file.path(out_chelsea, "htseq_VST_gene_names_top5000_variable_genes.csv"))
saveRDS(vst_mat_topvar, file.path(out_chelsea, "htseq_VST_matrix_top5000_variable_genes.rds"))

gene_variance_df <- data.table(Gene = names(vst_gene_vars), VST_variance = as.numeric(vst_gene_vars))
setorder(gene_variance_df, -VST_variance)
fwrite(gene_variance_df, file.path(out_chelsea, "htseq_VST_gene_variance_table.csv"))
fwrite(data.table(Sample = colnames(vst_mat)), file.path(out_chelsea, "VST_sample_names.csv"))

# GCT file for GenePattern NMFConsensus
gct_file <- file.path(out_chelsea, "vst_transformed_data.gct")
desc <- annot_out$ensembl_id[match(rownames(vst_mat), annot_out$Gene)]
con <- file(gct_file, "w")
writeLines("#1.2", con)
writeLines(paste(nrow(vst_mat), ncol(vst_mat), sep = "\t"), con)
close(con)
fwrite(data.table(Name = rownames(vst_mat), Description = desc, vst_mat),
       gct_file, sep = "\t", append = TRUE, col.names = TRUE)

stopifnot(!anyNA(vst_mat), all(is.finite(vst_mat)), min(vst_mat) >= 0)

############################################################
## 7. Summary
############################################################

cat("\n============================================================\n")
cat("DONE (", strand_setting, " counts, no batch correction)\n", sep = "")
cat("============================================================\n")
cat("Root (24 samples)    :", out_dir, "\n")
cat("Chelsea FFPE (21)    :", out_chelsea, "\n")
cat("Bionano fresh frozen :", out_bionano, "(raw counts only)\n")
cat("Genes after removing HTSeq summary rows:", nrow(counts_mat), "\n")
cat("Chelsea genes after filter             :", nrow(counts_filtered), "\n")
cat("Chelsea genes in top-variable matrix   :", nrow(vst_mat_topvar), "\n")
