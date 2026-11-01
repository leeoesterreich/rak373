# =============================================================================
# RNA_top5000_gene_class_check_GENCODEv47.R
#
# PURPOSE
#   Count how many of the top 5000 variable RNA features (the matrix already used
#   for MOFA2) are:
#     - unannotated  (only an Ensembl ID, no gene symbol)
#     - antisense    (antisense lncRNA)
#     - annotated    (all other genes, split by GENCODE gene type)
#   and show how many genes would remain under different restriction options
#   (for example "protein-coding + IG/TR genes only").
#
#   RESTRICTED SET (used for the restricted MOFA2 model in Biased_analysis):
#     KEEP    protein-coding, IG / TR genes, lncRNA that is not antisense
#     EXCLUDE antisense lncRNA, unannotated (Ensembl ID only), pseudogenes, rRNA,
#             snRNA, snoRNA and every other class not listed above, and all
#             mitochondrial genes (chrM or names starting MT-)
#   The allowed genes are written to RNA_genes_restricted_set.txt, which
#   01_Prepare_MOFA2_input_GENCODEv47.Rmd reads.
#
# HOW IT WORKS (summary)
#   1. Read the feature names from the top-5000 VST matrix (.rds).
#   2. Read the GENCODE v47 GTF (gene records only): gene ID, gene name, gene type,
#      position and strand.
#   3. Match every feature to the GTF (by gene symbol, or by Ensembl ID if the
#      feature has no symbol).
#   4. Flag antisense lncRNAs (by name and by opposite-strand overlap).
#   5. Assign each gene ONE class, count genes per class, save tables and a plot.
#   6. Optionally check which classes the strongest MOFA2 weights belong to.
#
# NOTES
#   - No random steps are used. SEED is kept only for consistency with other scripts.
#   - Run it section by section the first time and check the printed tables.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------
USER_LIB <- "/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(USER_LIB, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)   # fast reading of the GTF and table handling
  library(ggplot2)      # bar plot
})

SEED <- 42
set.seed(SEED)

# -----------------------------------------------------------------------------
# 1. PATHS AND SETTINGS  (the only part you normally need to edit)
# -----------------------------------------------------------------------------

# Output folder (created if missing)
OUT_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/MOFA2/Biased_analysis/RNA_gene_filter"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# RNA matrix to classify. It MUST be the same file as RNA_VST_FILE in script 01.
# This is the top-5000 matrix already used for the MOFA2 model. The restricted set is
# the part of these 5000 genes that passes the filter, so no new genes are added.
RNA_VST_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47/htseq_VST_matrix_top5000_variable_genes.rds"

# GENCODE v47 annotation = the GTF that was used for the HTSeq counts (.gtf or .gtf.gz).
# Confirm this is the GTF used for counting (full annotation, not the basic one).
GTF_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/gencode.v47.annotation.gtf"

# Optional: MOFA2 weight table written by the downstream script (section 9 uses it)
WEIGHTS_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/MOFA2/MOFA2_Downstream/MOFA2_feature_weights_long.tsv"
RNA_VIEW <- "RNA_VST_top5000var"   # name of the RNA view in the MOFA2 model
TOP_N_WEIGHTS_CHECK <- 100         # top weights per factor to classify in section 9

# TRUE = stop if GenomicRanges is missing, so antisense is never defined by name only
REQUIRE_OVERLAP_ANTISENSE <- TRUE

# --- Find the GTF if the path was not set ---
if (is.null(GTF_FILE) || !file.exists(GTF_FILE)) {
  candidates <- suppressWarnings(system(
    "find /ix1/alee/LO_LAB/Personal/Rahul -maxdepth 6 -iname '*gencode*v47*gtf*' 2>/dev/null | head -20",
    intern = TRUE))
  message("GTF_FILE is not set or the file was not found. Candidate files:\n",
          paste(candidates, collapse = "\n"))
  stop("Set GTF_FILE (section 1) to the GENCODE v47 GTF used for counting.")
}

# -----------------------------------------------------------------------------
# 2. READ THE FEATURE NAMES OF THE RNA MATRIX
# -----------------------------------------------------------------------------
# In your matrix, rows are gene symbols; genes WITHOUT a symbol keep their Ensembl ID
# (for example ENSG00000289983). This is how "unannotated" genes are recognised.

rna_obj <- readRDS(RNA_VST_FILE)

if (is.data.frame(rna_obj)) {
  # Data frame: the gene names may be in a column or in the row names
  rna_dt <- as.data.table(rna_obj)
  gene_col <- intersect(c("Gene", "gene", "gene_name", "Hugo_Symbol",
                          "hgnc_symbol", "symbol", "ensembl_gene_id"),
                        colnames(rna_dt))[1]
  features <- if (!is.na(gene_col)) as.character(rna_dt[[gene_col]]) else rownames(rna_obj)
} else {
  # Matrix: gene names are the row names
  features <- rownames(rna_obj)
}
stopifnot(!is.null(features))
features <- trimws(features)

message("Features in the RNA matrix: ", length(features),
        " | unique names: ", uniqueN(features))
print(head(features))

# -----------------------------------------------------------------------------
# 3. READ GENE RECORDS FROM THE GENCODE GTF
# -----------------------------------------------------------------------------
# 'zcat -f' reads both .gtf and .gtf.gz; 'grep -v ^#' removes the header lines.
gtf <- fread(cmd = paste("zcat -f", shQuote(GTF_FILE), "| grep -v '^#'"),
             sep = "\t", header = FALSE, quote = "")
gtf <- gtf[V3 == "gene"]          # keep gene records only (not transcripts/exons)

# Helper: pull one field (for example gene_type) out of the GTF attribute column
get_attr <- function(x, key) {
  has <- grepl(paste0(key, ' "'), x, fixed = TRUE)
  out <- rep(NA_character_, length(x))
  out[has] <- sub(paste0('.*', key, ' "([^"]+)".*'), "\\1", x[has])
  out
}

anno <- data.table(
  gene_id   = get_attr(gtf$V9, "gene_id"),
  gene_name = get_attr(gtf$V9, "gene_name"),
  gene_type = get_attr(gtf$V9, "gene_type"),
  chr = gtf$V1, start = gtf$V4, end = gtf$V5, strand = gtf$V7
)

# Ensembl ID without version: ENSG00000289983.1 -> ENSG00000289983
# (this also removes the "_PAR_Y" suffix of pseudo-autosomal genes)
anno[, gene_id_nov := sub("\\..*$", "", gene_id)]

message("Genes in the GTF: ", nrow(anno))

# Check which gene types exist in THIS GENCODE release (look for "antisense")
print(anno[, .N, by = gene_type][order(-N)])

# -----------------------------------------------------------------------------
# 4. FLAG ANTISENSE lncRNAs
# -----------------------------------------------------------------------------
# GENCODE v47 has no separate "antisense" gene type; antisense genes are lncRNAs.
# Two definitions are combined:
#   (1) NAME:    lncRNA whose name ends in -AS, -AS1, -AS2 ...   (e.g. THBS1-AS1)
#   (2) OVERLAP: lncRNA that overlaps a protein-coding gene on the OPPOSITE strand
#                (needs the Bioconductor package GenomicRanges)
anno[, antisense_name    := gene_type == "lncRNA" & grepl("-AS[0-9]*$", gene_name)]
anno[, antisense_overlap := FALSE]

if (requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("IRanges", quietly = TRUE) &&
    requireNamespace("S4Vectors", quietly = TRUE)) {

  gr  <- GenomicRanges::GRanges(anno$chr,
                                IRanges::IRanges(anno$start, anno$end),
                                strand = anno$strand)
  pc  <- which(anno$gene_type == "protein_coding")
  lnc <- which(anno$gene_type == "lncRNA")

  # all lncRNA - protein-coding pairs that overlap (strand ignored for now)
  ov <- GenomicRanges::findOverlaps(gr[lnc], gr[pc], ignore.strand = TRUE)
  q  <- lnc[S4Vectors::queryHits(ov)]      # lncRNA row numbers
  s  <- pc[S4Vectors::subjectHits(ov)]     # overlapping protein-coding row numbers

  # keep only the pairs on opposite strands = antisense by overlap
  opposite <- anno$strand[q] != anno$strand[s]
  anno[unique(q[opposite]), antisense_overlap := TRUE]
} else {
  if (REQUIRE_OVERLAP_ANTISENSE) {
    stop("GenomicRanges/IRanges/S4Vectors not available. Install them, or set ",
         "REQUIRE_OVERLAP_ANTISENSE <- FALSE to define antisense by gene name only.")
  }
  message("GenomicRanges not available: antisense is defined by gene name only.")
}

anno[, is_antisense := antisense_name | antisense_overlap]

# How many lncRNAs are antisense by name, by overlap, or both (whole genome)
print(anno[gene_type == "lncRNA", .N, by = .(antisense_name, antisense_overlap)])

# -----------------------------------------------------------------------------
# 5. MATCH THE FEATURES TO THE ANNOTATION
# -----------------------------------------------------------------------------
# Three kinds of feature names occur in the matrix:
#   1. plain Ensembl ID                  ENSG00000289983      -> matched by gene_id
#   2. symbol with an Ensembl ID suffix  MATR3|ENSG00000015479.20
#      (used when one symbol belongs to several Ensembl IDs)  -> matched by gene_id
#   3. plain gene symbol                 TSPAN6               -> matched by gene_name
# If a gene name occurs twice in the GTF, one record is kept (non-PAR_Y preferred).
anno_name <- anno[order(grepl("_PAR_Y", gene_id))][!duplicated(gene_name)]
anno_id   <- anno[!duplicated(gene_id_nov)]

feat_dt <- data.table(feature = features)
feat_dt[, is_ensembl_id := grepl("^ENSG[0-9]+(\\.[0-9]+)?$", feature)]        # type 1 only
feat_dt[, has_id_suffix := grepl("\\|ENSG[0-9]+(\\.[0-9]+)?$", feature)]      # type 2
feat_dt[, match_by_id   := is_ensembl_id | has_id_suffix]
# Ensembl ID without the symbol prefix and without the version
feat_dt[, key_id := sub("\\.[0-9]+$", "", sub("^.*\\|", "", feature))]

message("Features written as SYMBOL|ENSEMBL_ID (matched by Ensembl ID): ", feat_dt[, sum(has_id_suffix)])

# Row number of the matching GTF record (NA = not found)
feat_dt[, idx := fifelse(match_by_id,
                         match(key_id, anno_id$gene_id_nov),
                         match(feature, anno_name$gene_name))]

# Copy a column from the right annotation table (by ID or by name)
pick <- function(col) {
  fifelse(feat_dt$match_by_id,
          anno_id[[col]][feat_dt$idx],
          anno_name[[col]][feat_dt$idx])
}
feat_dt[, `:=`(gene_id           = pick("gene_id"),
               gtf_gene_name     = pick("gene_name"),
               chr               = pick("chr"),
               gene_type         = pick("gene_type"),
               is_antisense      = pick("is_antisense"),
               antisense_name    = pick("antisense_name"),
               antisense_overlap = pick("antisense_overlap"))]
feat_dt[is.na(is_antisense), is_antisense := FALSE]
feat_dt[, idx := NULL]

message("Features not found in the GTF: ", feat_dt[is.na(gene_type), .N])

# -----------------------------------------------------------------------------
# 6. ASSIGN ONE CLASS PER GENE
# -----------------------------------------------------------------------------
# Order of the rules matters (first match wins):
#   1. no GTF match                -> "Not found in GTF"
#   2. Ensembl ID only             -> "Unannotated (Ensembl ID only)"
#   3. lncRNA and antisense        -> "Antisense lncRNA"
#   4. other lncRNA                -> "Other lncRNA"
#   5. any pseudogene type         -> "Pseudogene"
#   6. IG_* / TR_* gene types      -> "IG / TR gene"   (immunoglobulin / T-cell receptor)
#   7. protein_coding              -> "Protein-coding"
#   8. everything else (miRNA, snoRNA, TEC ...) -> "Other ncRNA / other"
# NOTE: IG/TR genes (IGHG1, IGKC ...) are NOT labelled protein_coding in GENCODE,
#       so they need their own class or they would be lost by a coding-only filter.
feat_dt[, primary_class :=
  fifelse(is.na(gene_type), "Not found in GTF",
  fifelse(is_ensembl_id, "Unannotated (Ensembl ID only)",
  fifelse(gene_type == "lncRNA" & is_antisense, "Antisense lncRNA",
  fifelse(gene_type == "lncRNA", "Other lncRNA",
  fifelse(grepl("pseudogene", gene_type), "Pseudogene",
  fifelse(grepl("^(IG|TR)_", gene_type), "IG / TR gene",
  fifelse(gene_type == "protein_coding", "Protein-coding",
          "Other ncRNA / other")))))))]

# Simple three-way split (as asked): unannotated / antisense / other annotated
feat_dt[, three_way :=
  fifelse(primary_class == "Unannotated (Ensembl ID only)", "Unannotated",
  fifelse(primary_class == "Antisense lncRNA", "Antisense", "Annotated (other)"))]

# Flags for the possible restriction options
feat_dt[, keep_annotated := !primary_class %in% c("Unannotated (Ensembl ID only)", "Not found in GTF")]
feat_dt[, keep_PC_IG_TR  := primary_class %in% c("Protein-coding", "IG / TR gene")]
feat_dt[, keep_PC_only   := primary_class == "Protein-coding"]

# Restricted set for the restricted MOFA2 model:
#   keep protein-coding, IG/TR and non-antisense lncRNA; remove mitochondrial genes.
# Mitochondrial genes: chromosome chrM or a name starting with MT-
# (mitochondrial protein-coding genes are labelled protein_coding in GENCODE).
feat_dt[, is_mito := (!is.na(chr) & chr == "chrM") |
                     (!is.na(gtf_gene_name) & grepl("^MT-", gtf_gene_name))]
feat_dt[, keep_restricted := primary_class %in% c("Protein-coding", "IG / TR gene", "Other lncRNA") &
                             !is_mito]
# Reason for every excluded gene (mitochondrial is reported first)
feat_dt[, restricted_status := fifelse(keep_restricted, "Kept",
                               fifelse(is_mito, "Excluded: mitochondrial",
                                       paste0("Excluded: ", primary_class)))]

# Per-gene table (use this to filter genes later)
fwrite(feat_dt, file.path(OUT_DIR, "RNA_gene_class_per_gene.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 7. SUMMARY TABLES
# -----------------------------------------------------------------------------
n_total <- nrow(feat_dt)

# 7a. Genes per class
class_summary <- feat_dt[, .(n_genes = .N), by = primary_class][order(-n_genes)]
class_summary[, percent := round(100 * n_genes / n_total, 1)]
print(class_summary)
fwrite(class_summary, file.path(OUT_DIR, "RNA_gene_class_summary.tsv"), sep = "\t")

# 7b. The three-way split you asked for
three_way_summary <- feat_dt[, .(n_genes = .N), by = three_way][order(-n_genes)]
three_way_summary[, percent := round(100 * n_genes / n_total, 1)]
print(three_way_summary)
fwrite(three_way_summary,
       file.path(OUT_DIR, "RNA_unannotated_antisense_annotated_summary.tsv"), sep = "\t")

# 7c. Exact GENCODE gene types behind each class
biotype_table <- feat_dt[, .N, by = .(primary_class, gene_type)][order(primary_class, -N)]
print(biotype_table)
fwrite(biotype_table, file.path(OUT_DIR, "RNA_gene_type_breakdown.tsv"), sep = "\t")

# 7d. How many genes remain under each restriction option
restriction_options <- data.table(
  option = c("All genes (current model)",
             "Annotated genes only (drop Ensembl-ID-only genes)",
             "Protein-coding + IG/TR genes",
             "Protein-coding only",
             "RESTRICTED SET (protein-coding + IG/TR + non-antisense lncRNA, no mitochondrial)"),
  n_genes_kept = c(n_total,
                   feat_dt[, sum(keep_annotated)],
                   feat_dt[, sum(keep_PC_IG_TR)],
                   feat_dt[, sum(keep_PC_only)],
                   feat_dt[, sum(keep_restricted)])
)
restriction_options[, percent_kept := round(100 * n_genes_kept / n_total, 1)]
print(restriction_options)
fwrite(restriction_options, file.path(OUT_DIR, "RNA_restriction_options.tsv"), sep = "\t")

# 7e. Gene lists for each option (one gene per line), ready to use as a filter
writeLines(feat_dt[keep_annotated == TRUE, feature],
           file.path(OUT_DIR, "RNA_genes_annotated_only.txt"))
writeLines(feat_dt[keep_PC_IG_TR == TRUE, feature],
           file.path(OUT_DIR, "RNA_genes_protein_coding_plus_IG_TR.txt"))
writeLines(feat_dt[keep_PC_only == TRUE, feature],
           file.path(OUT_DIR, "RNA_genes_protein_coding_only.txt"))

# Restricted set: this file is read by script 01 (one feature name per line, exactly as
# the names appear in the RNA matrix)
writeLines(feat_dt[keep_restricted == TRUE, feature],
           file.path(OUT_DIR, "RNA_genes_restricted_set.txt"))

# Genes kept / excluded, with the reason (paired table for the restricted set)
restricted_summary <- feat_dt[, .(n_genes = .N), by = restricted_status][order(-n_genes)]
restricted_summary[, percent := round(100 * n_genes / n_total, 1)]
print(restricted_summary)
fwrite(restricted_summary, file.path(OUT_DIR, "RNA_restricted_set_summary.tsv"), sep = "\t")

# Mitochondrial genes found (check)
print(feat_dt[is_mito == TRUE, .(feature, gene_type, chr)])

# 7f. Immunoglobulin / T-cell receptor genes in the matrix (lost by a coding-only filter)
print(feat_dt[primary_class == "IG / TR gene", .(feature, gene_type)][order(gene_type, feature)])

# -----------------------------------------------------------------------------
# 8. BAR PLOT OF GENES PER CLASS  (saved as PDF and PNG)
# -----------------------------------------------------------------------------
plot_dt <- copy(class_summary)
plot_dt[, primary_class := factor(primary_class, levels = rev(primary_class))]

p_class <- ggplot(plot_dt, aes(x = n_genes, y = primary_class)) +
  geom_col(fill = "#3C78A8", width = 0.7) +
  geom_text(aes(label = paste0(n_genes, " (", percent, "%)")), hjust = -0.1, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +   # room for the labels
  labs(title = paste0("Gene classes among the ", n_total, " RNA features"),
       x = "Number of genes", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot")

print(p_class)

pdf(file.path(OUT_DIR, "RNA_gene_class_counts.pdf"), width = 8, height = 4.5)
print(p_class); dev.off()
png(file.path(OUT_DIR, "RNA_gene_class_counts.png"),
    width = 8, height = 4.5, units = "in", res = 300)
print(p_class); dev.off()

# -----------------------------------------------------------------------------
# 9. WHICH CLASSES DO THE STRONGEST MOFA2 WEIGHTS BELONG TO?  (optional)
# -----------------------------------------------------------------------------
# Shows whether the factors depend on non-coding genes. Skipped if the weight table
# from the downstream script is not available.
if (file.exists(WEIGHTS_FILE)) {

  w <- fread(WEIGHTS_FILE)
  w <- w[view == RNA_VIEW]

  # use the cleaned feature name if the downstream script created it
  feat_col <- if ("feature_clean" %in% colnames(w)) "feature_clean" else "feature"
  w[, feature_use := as.character(get(feat_col))]

  # attach the gene class
  w <- merge(w, feat_dt[, .(feature_use = feature, primary_class)],
             by = "feature_use", all.x = TRUE)

  # top N genes per factor by absolute weight, then count genes per class
  top_w <- w[order(factor, -abs_weight)][, .SD[seq_len(min(.N, TOP_N_WEIGHTS_CHECK))], by = factor]
  weight_class <- top_w[, .N, by = .(factor, primary_class)]
  weight_class[, percent := round(100 * N / sum(N), 1), by = factor]
  weight_class <- weight_class[order(factor, -N)]

  print(weight_class)
  fwrite(weight_class,
         file.path(OUT_DIR, paste0("Top", TOP_N_WEIGHTS_CHECK, "_RNA_weights_gene_class_by_factor.tsv")),
         sep = "\t")
} else {
  message("Weights file not found, section 9 skipped: ", WEIGHTS_FILE)
}

# -----------------------------------------------------------------------------
# 10. SESSION INFORMATION
# -----------------------------------------------------------------------------
sessionInfo()