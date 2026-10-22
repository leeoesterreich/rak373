############################################################
## PURPOSE
##   Oncoplots (WES, non-silent, tumor VAF >= 3%) of the leading-edge genes
##   from the estrogen-vs-chromatin correlation analysis, by RNA cluster.
##   One oncoplot for the chromatin gene set and one for the estrogen gene set.
##
## LAYOUT (same style as the response-to-estradiol oncoplot)
##   - Genes: leading-edge genes (top N by |rho| with the other family's score);
##     only genes mutated in at least one tumor are drawn, sorted by frequency.
##   - Columns: all WES tumors, split by RNA cluster (1, 2, 3, No RNA).
##   - Cell colour = variant class; Multi_Hit = more than one variant in the
##     gene in that tumor.
##   - Stars = PCGR Oncogenic, triangles = PCGR Likely_Oncogenic (exact tumor
##     + variant match to the PCGR workbooks).
##   - Top: mutations per sample (stacked by class) and the RNA cluster strip.
##     Right: altered samples by class and % altered. Bottom: clinical annotation.
##
## INPUT
##   - chromatin_leading_edge_genes.tsv and estrogen_leading_edge_genes.tsv
##     written by pILC_chromatin_estrogen_correlation.R (run that script first)
##
## NOTE
##   The per-gene cluster test is exploratory (few mutated tumors). Clusters and
##   leading-edge genes were derived from the same expression data.
############################################################

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
graphics.off()

# ----------------------------
# Paths (fixed names, no date or time stamps in any output)
# ----------------------------
root          <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"
wes_audit     <- file.path(root, paste0("WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/",
                                        "Combined_MAF/Full_MAF_nonsilent_variant_PCGR_VAF_match_audit_2026-09-23.tsv"))
pcgr_dir      <- file.path(root, "WES/Downstream/PCGR/With_PON/pcgr_outputs")
cluster_csv_file <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
cluster_gct_file <- file.path(root, "RNASeq/NMF/GenePattern_NMFConsensus/vst_transformed_data.consensus.k.3.gct")  # fallback
clinical_file <- file.path(root, "Clinical_correlation/pILC_ClinicoPathological_details.csv")

corr_dir   <- file.path(root, "RNASeq/Downstream/Chromatin_Estrogen_correlation")
gene_files <- c(Chromatin = file.path(corr_dir, "tables", "chromatin_leading_edge_genes.tsv"),
                Estrogen  = file.path(corr_dir, "tables", "estrogen_leading_edge_genes.tsv"))

out_dir   <- file.path(corr_dir, "oncoplots")
table_dir <- file.path(out_dir, "tables")
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(wes_audit, clinical_file, gene_files)) {
  if (!file.exists(f)) stop("Input file not found: ", f,
                            "\n(The leading-edge gene tables come from pILC_chromatin_estrogen_correlation.R.)")
}
if (!dir.exists(pcgr_dir)) stop("PCGR folder not found: ", pcgr_dir)
if (!file.exists(cluster_csv_file) && !file.exists(cluster_gct_file)) stop("No cluster file found.")

# ----------------------------
# Settings
# ----------------------------
N_GENES <- 60                   # leading-edge genes per set: top N by |rho| with the other family's score (NULL = all)
MIN_VAF <- 0.03
VAF_TAG <- paste0("VAFge", round(100 * MIN_VAF), "pct")
SHOW_UNMUTATED_GENES <- FALSE   # TRUE = also draw genes with no mutation
CELL_W_MM <- 4.6
CELL_H_MM <- 5.8
SAVE_PNG  <- TRUE

nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Splice_Region", "Translation_Start_Site"
)
cluster_cols <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3", "No RNA" = "#BDBDBD")
class_cols <- c(
  "Missense_Mutation" = "#1F77B4", "Nonsense_Mutation" = "#D62728", "Frame_Shift_Del" = "#9467BD",
  "Frame_Shift_Ins" = "#8C564B", "In_Frame_Del" = "#2CA02C", "In_Frame_Ins" = "#17BECF",
  "Splice_Site" = "#FF7F0E", "Splice_Region" = "#FFBB78", "Translation_Start_Site" = "#E377C2",
  "Nonstop_Mutation" = "#BCBD22", "Multi_Hit" = "#000000"
)

# ----------------------------
# Helpers
# ----------------------------
std_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^X(?=TP)", "", x, perl = TRUE)
  x <- sub("_S[0-9]+$", "", x)
  x <- gsub("\\.", "-", x)
  x <- sub("^JC287-", "", x)
  x <- sub("(_T|-T)$", "", x)
  x <- sub("(?<=[0-9])T$", "", x, perl = TRUE)
  sub("-0+(?=[0-9])", "-", x, perl = TRUE)      # TP09-09 = TP09-9
}
natural_sample_order <- function(s) {
  s <- unique(as.character(s))
  tp <- suppressWarnings(as.integer(sub("^TP([0-9]+).*", "\\1", s))); tp[is.na(tp)] <- 9999L
  rest <- suppressWarnings(as.integer(gsub("[^0-9]", "", sub("^TP[0-9]+[-_]*", "", s)))); rest[is.na(rest)] <- 9999L
  s[order(tp, rest, s)]
}
parse_percent_value <- function(x) {
  x <- trimws(as.character(x)); x[x %in% c("", "NA", "N/A", "NaN", "nan")] <- NA_character_
  lt <- grepl("^<", x)
  v <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  v[lt & !is.na(v)] <- v[lt & !is.na(v)] / 2
  v
}
clean_receptor_status <- function(x) {
  x <- gsub("\\s+", " ", trimws(as.character(x)))
  out <- rep("Unknown", length(x))
  out[grepl("^Positive", x, ignore.case = TRUE)] <- "Positive"
  out[grepl("^Negative", x, ignore.case = TRUE)] <- "Negative"
  out
}
clean_grade_status <- function(x) {
  x <- gsub("\\s+", " ", trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NaN", "nan") | is.na(x)] <- "Unknown"
  x
}
normalize_oncogenicity <- function(x) {
  x <- toupper(gsub("[ -]+", "_", trimws(as.character(x))))
  out <- rep(NA_character_, length(x))
  out[x == "ONCOGENIC"] <- "Oncogenic"
  out[x == "LIKELY_ONCOGENIC"] <- "Likely_Oncogenic"
  out
}

# Variant key chr:pos:ref:alt; identical PREFIX removed first (shifts pos), then suffix.
edit_key <- function(chr, pos, ref, alt, maf_format = FALSE) {
  chr <- sub("^CHR", "", toupper(trimws(as.character(chr))))
  pos <- suppressWarnings(as.integer(pos))
  ref <- toupper(trimws(as.character(ref))); alt <- toupper(trimws(as.character(alt)))
  if (anyNA(c(chr, pos, ref, alt)) || !nzchar(chr) || pos < 1L ||
      !grepl("^([ACGTN]+|-)$", ref) || !grepl("^([ACGTN]+|-)$", alt)) return(NA_character_)
  if (maf_format && ref == "-") pos <- pos + 1L
  if (ref == "-") ref <- ""
  if (alt == "-") alt <- ""
  while (nzchar(ref) && nzchar(alt) && substr(ref, 1, 1) == substr(alt, 1, 1)) {
    ref <- substring(ref, 2); alt <- substring(alt, 2); pos <- pos + 1L
  }
  while (nzchar(ref) && nzchar(alt) && substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))) {
    ref <- substr(ref, 1, nchar(ref) - 1L); alt <- substr(alt, 1, nchar(alt) - 1L)
  }
  if (!nzchar(ref) && !nzchar(alt)) return(NA_character_)
  paste(chr, pos, ifelse(nzchar(ref), ref, "-"), ifelse(nzchar(alt), alt, "-"), sep = ":")
}
pcgr_genomic_key <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  hit <- regmatches(x, regexec("^(?:chr)?([^:]+):g\\.([0-9]+)([ACGTN]+|-)>([ACGTN]+|-)$", x, perl = TRUE, ignore.case = TRUE))[[1]]
  if (length(hit) == 5L) return(edit_key(hit[2], hit[3], hit[4], hit[5]))
  NA_character_
}
save_plot <- function(file_base, width, height, draw_fn, res = 300) {
  pdf(paste0(file_base, ".pdf"), width = width, height = height); draw_fn(); dev.off()
  if (isTRUE(SAVE_PNG)) {
    png(paste0(file_base, ".png"), width = width, height = height, units = "in", res = res); draw_fn(); dev.off()
  }
  cat("Saved:", paste0(file_base, c(".pdf", ".png")), sep = "\n  "); cat("\n")
}

# ----------------------------
# 1. Leading-edge gene lists (both families)
# ----------------------------
read_genes <- function(file) {
  g <- fread(file)
  if (!("Gene" %in% names(g))) stop("No 'Gene' column in ", file)
  g[, Gene := toupper(trimws(Gene))]
  if ("Rho_with_other_score" %in% names(g)) g <- g[order(-abs(Rho_with_other_score))]
  out <- unique(g$Gene)
  if (!is.null(N_GENES)) out <- head(out, N_GENES)
  out
}
gene_sets <- lapply(gene_files, read_genes)
cat("Leading-edge genes:", paste(names(gene_sets), lengths(gene_sets), sep = " = ", collapse = " | "), "\n")
all_genes <- unique(unlist(gene_sets))

# ----------------------------
# 2. WES calls: non-silent, VAF >= cutoff, in any leading-edge gene (read once)
# ----------------------------
wes <- fread(wes_audit, sep = "\t", na.strings = c("", "NA", "."),
             select = c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "VAF_TUMOR_fraction",
                        "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
wes[, Tumor := std_id(Tumor_Sample_Barcode)]
wes_tumors <- natural_sample_order(unique(wes$Tumor))
cat("Tumors with WES:", length(wes_tumors), "\n")

calls_all <- wes[toupper(Hugo_Symbol) %in% all_genes & Variant_Classification %in% nonsilent_classes &
                   !is.na(VAF_TUMOR_fraction) & VAF_TUMOR_fraction >= MIN_VAF]
rm(wes)
calls_all[, `:=`(Gene = toupper(Hugo_Symbol), Class = Variant_Classification, VAF = VAF_TUMOR_fraction)]
calls_all[, Variant_key := mapply(edit_key, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
                                  MoreArgs = list(maf_format = TRUE), USE.NAMES = FALSE)]
calls_all[, Match_key := paste(Tumor, Variant_key, sep = "|")]
cat("Non-silent variants in the leading-edge genes at VAF >=", MIN_VAF * 100, "%:", nrow(calls_all), "\n")
if (nrow(calls_all) == 0) stop("No WES mutations found in the gene sets at this cutoff.")

# ----------------------------
# 3. PCGR oncogenicity (exact tumor + variant match; read once)
# ----------------------------
pcgr_files <- list.files(pcgr_dir, pattern = "\\.pcgr\\.grch38\\.xlsx$", recursive = TRUE, full.names = TRUE)
cat("PCGR workbooks found:", length(pcgr_files), "\n")
onc <- rbindlist(lapply(pcgr_files, function(f) {
  if (!"SOMATIC_SNV_INDEL" %in% tryCatch(excel_sheets(f), error = function(e) character())) return(NULL)
  dt <- tryCatch(as.data.table(read_excel(f, sheet = "SOMATIC_SNV_INDEL", guess_max = 100000)), error = function(e) NULL)
  if (is.null(dt) || !all(c("GENOMIC_CHANGE", "ONCOGENICITY") %in% names(dt)) || nrow(dt) == 0) return(NULL)
  data.table(Tumor = std_id(basename(dirname(f))),
             Variant_key = vapply(dt$GENOMIC_CHANGE, pcgr_genomic_key, character(1)),
             ONC = normalize_oncogenicity(dt$ONCOGENICITY))
}))
onc <- onc[!is.na(Variant_key) & !is.na(ONC)]
onc[, rank := fifelse(ONC == "Oncogenic", 2L, 1L)]
onc_best <- onc[order(-rank), .(ONC = ONC[1]), by = .(Match_key = paste(Tumor, Variant_key, sep = "|"))]
calls_all <- merge(calls_all, onc_best, by = "Match_key", all.x = TRUE, sort = FALSE)
cat("Variants with a PCGR Oncogenic / Likely_Oncogenic call:\n"); print(calls_all[!is.na(ONC), .N, by = ONC])

# ----------------------------
# 4. RNA clusters (same file as the correlation script)
# ----------------------------
if (file.exists(cluster_csv_file)) {
  cc <- fread(cluster_csv_file)
  samp_col <- intersect(c("sample", "Sample", "Name", "name"), names(cc))
  samp_col <- if (length(samp_col)) samp_col[1] else names(cc)[1]
  clus_col <- setdiff(grep("cluster|membership", names(cc), ignore.case = TRUE, value = TRUE), samp_col)[1]
  if (is.na(clus_col)) stop("No cluster column found in ", cluster_csv_file)
  message("Cluster CSV: sample column '", samp_col, "', cluster column '", clus_col, "'")
  cl <- data.table(Tumor = std_id(cc[[samp_col]]),
                   Cluster = paste("Cluster", as.integer(gsub("[^0-9]", "", as.character(cc[[clus_col]])))))
} else {
  message("Cluster CSV not found; using the GCT file: ", cluster_gct_file)
  cg <- fread(cluster_gct_file, skip = 2)
  cl_col <- grep("membership", names(cg), value = TRUE)[1]
  cl <- cg[, .(Tumor = std_id(Name), Cluster = paste("Cluster", get(cl_col)))]
}
grp <- cl$Cluster[match(wes_tumors, cl$Tumor)]
grp[is.na(grp)] <- "No RNA"
cat("WES tumors by RNA cluster:\n"); print(table(grp))

tumor_dt <- data.table(Tumor = wes_tumors, Group = factor(grp, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "No RNA")))
tumor_dt <- tumor_dt[order(Group, match(Tumor, wes_tumors))]
tumors <- tumor_dt$Tumor
group  <- droplevels(tumor_dt$Group)

# ----------------------------
# 5. Clinical annotation
# ----------------------------
clin <- fread(clinical_file, sep = ",", header = TRUE, data.table = TRUE, fill = TRUE,
              na.strings = c("", "NA", "NaN", "N/A"))
clin[, Tumor := std_id(Sample)]
clin <- unique(clin[, .(Tumor,
  tumor_cellularity_pct = parse_percent_value(`tumor_cellularity (%)`),
  pleomorphic_pct = parse_percent_value(`pleomorphic_invasive_compoent_cellularity (%)`),
  TILs_pct = parse_percent_value(`TILs (%)`),
  Grade = clean_grade_status(`Grade/Differentiation-Desc`),
  HER2 = clean_receptor_status(`BREAST: Her-2/Neu-Desc`),
  ER = clean_receptor_status(`Tumor Marker #1 (ER)`),
  PR = clean_receptor_status(`Tumor Marker #2 (PR)`))], by = "Tumor")
clin <- clin[match(tumors, Tumor)]
cat("Tumors matched to the clinical file:", sum(tumors %in% clin$Tumor), "of", length(tumors), "\n")

# Clinical annotation colours (shared by both oncoplots)
grade_cols <- c("Grade I" = "#A6CEE3", "Grade II" = "#33A02C", "Grade III" = "#1F78B4",
                "Cell type not determined" = "#E31A1C", "Unknown" = "#BDBDBD")
rec_cols <- c("Negative" = "#1F78B4", "Positive" = "#E31A1C", "Unknown" = "#BDBDBD")
extra_g <- setdiff(unique(na.omit(clin$Grade)), names(grade_cols))
if (length(extra_g)) { eg <- hcl.colors(length(extra_g), "Dark 3"); names(eg) <- extra_g; grade_cols <- c(grade_cols, eg) }
anno_df <- data.frame(
  `Tumor cell. (%)` = clin$tumor_cellularity_pct, `Pleomorphic (%)` = clin$pleomorphic_pct,
  `TILs (%)` = clin$TILs_pct, Grade = clin$Grade, HER2 = clin$HER2, ER = clin$ER, PR = clin$PR,
  check.names = FALSE, row.names = tumors
)

# ----------------------------
# 6. One oncoplot per gene set
# ----------------------------
make_oncoplot <- function(family, leading_edge_genes) {
  GENESET_LABEL <- paste0(family, "_leading_edge")
  calls <- calls_all[Gene %in% leading_edge_genes]
  cat("\n=====", family, ": non-silent variants at VAF >=", MIN_VAF * 100, "%:", nrow(calls), "=====\n")
  if (nrow(calls) == 0) { message("No WES mutations in the ", family, " gene set; no oncoplot drawn."); return(invisible(NULL)) }

  cell <- calls[, .(n = .N, Class = if (.N > 1) "Multi_Hit" else Class[1],
                    ONC = if ("Oncogenic" %in% ONC) "Oncogenic" else if ("Likely_Oncogenic" %in% ONC) "Likely_Oncogenic" else NA_character_,
                    Variants = paste(unique(Variant_key), collapse = ";"),
                    VAFs = paste(round(100 * VAF, 1), collapse = ";")), by = .(Tumor, Gene)]

  # genes: only mutated ones (unless SHOW_UNMUTATED_GENES), most frequent first
  gfreq <- merge(data.table(Gene = leading_edge_genes), cell[, .(n_tumors = uniqueN(Tumor), n_var = sum(n)), by = Gene],
                 by = "Gene", all.x = TRUE)
  gfreq[is.na(n_tumors), `:=`(n_tumors = 0L, n_var = 0L)]
  gfreq <- gfreq[order(-n_tumors, -n_var, match(Gene, leading_edge_genes))]
  cat("Genes without a mutation:", sum(gfreq$n_tumors == 0), "of", nrow(gfreq), "\n")
  genes <- if (SHOW_UNMUTATED_GENES) gfreq$Gene else gfreq[n_tumors > 0, Gene]

  mat <- matrix("None", nrow = length(genes), ncol = length(tumors), dimnames = list(genes, tumors))
  onc_mat <- matrix("", nrow = length(genes), ncol = length(tumors), dimnames = list(genes, tumors))
  d <- cell[Gene %in% genes & Tumor %in% tumors]
  mat[cbind(d$Gene, d$Tumor)] <- d$Class
  onc_mat[cbind(d$Gene, d$Tumor)] <- fifelse(is.na(d$ONC), "", d$ONC)
  mutated <- mat != "None"

  # per-gene counts by cluster and an exploratory Fisher test across clusters 1-3 (tumors with RNA only)
  rna <- group != "No RNA"
  st <- rbindlist(lapply(genes, function(g) {
    m <- mutated[g, ]
    row <- data.table(Gene = g, Mutated_all = sum(m), Tumors_all = length(m))
    for (k in levels(group)) {
      row[[paste0(sub(" ", "", k), "_mutated")]] <- sum(m[group == k]); row[[paste0(sub(" ", "", k), "_n")]] <- sum(group == k)
    }
    row[["Fisher_p_across_clusters_1_3"]] <- tryCatch(
      fisher.test(table(factor(m[rna], c(FALSE, TRUE)), droplevels(group[rna])))$p.value, error = function(e) NA_real_)
    row
  }), fill = TRUE)
  st[, BH_adjusted_p := p.adjust(Fisher_p_across_clusters_1_3, method = "BH")]
  fwrite(st,   file.path(table_dir, paste0(GENESET_LABEL, "_WES_", VAF_TAG, "_mutations_by_cluster.tsv")), sep = "\t")
  fwrite(cell, file.path(table_dir, paste0(GENESET_LABEL, "_WES_", VAF_TAG, "_cell_level_table.tsv")), sep = "\t", na = "")
  print(st[, .(Gene, Mutated_all, p = signif(Fisher_p_across_clusters_1_3, 3), BH = signif(BH_adjusted_p, 3))])

  used <- intersect(names(class_cols), unique(as.vector(mat)))

  # ---- top: mutations per sample (stacked by class) + RNA cluster strip
  top_mat <- do.call(cbind, lapply(used, function(cn) colSums(mat == cn)))
  colnames(top_mat) <- used; rownames(top_mat) <- tumors
  top_ha <- HeatmapAnnotation(
    `Mutations per sample` = anno_barplot(top_mat, gp = gpar(fill = class_cols[used], col = NA), border = FALSE,
                                          bar_width = 0.97, height = unit(2.2, "cm"), axis_param = list(gp = gpar(fontsize = 8))),
    `RNA cluster` = as.character(group),
    col = list(`RNA cluster` = cluster_cols),
    annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 9),
    simple_anno_size = unit(0.35, "cm"), gap = unit(1.5, "mm"),
    show_legend = TRUE
  )

  # ---- right: altered samples by class and % altered
  bar_mat <- do.call(cbind, lapply(used, function(cn) rowSums(mat == cn)))
  colnames(bar_mat) <- used; rownames(bar_mat) <- genes
  pct <- paste0(round(100 * rowSums(mutated) / ncol(mutated)), "%")
  right_ha <- rowAnnotation(
    `Altered samples by class` = anno_barplot(bar_mat, gp = gpar(fill = class_cols[used], col = NA), border = FALSE,
                                              bar_width = 0.90, axis = TRUE, axis_param = list(gp = gpar(fontsize = 7), side = "bottom")),
    `% altered` = anno_text(pct, gp = gpar(fontsize = 10), just = "left", location = 0),
    width = unit(4.2, "cm"), show_annotation_name = FALSE
  )

  # ---- bottom: clinical annotation
  bottom_ha <- HeatmapAnnotation(
    df = anno_df,
    col = list(`Tumor cell. (%)` = colorRamp2(c(0, 50, 100), c("#F7FBFF", "#6BAED6", "#08306B")),
               `Pleomorphic (%)` = colorRamp2(c(0, 50, 100), c("#F7FCF5", "#74C476", "#00441B")),
               `TILs (%)` = colorRamp2(c(0, 30, 100), c("#FFF5F0", "#FB6A4A", "#67000D")),
               Grade = grade_cols, HER2 = rec_cols, ER = rec_cols, PR = rec_cols),
    na_col = "grey90", annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(0.30, "cm"), gap = unit(1.1, "mm")
  )

  # ---- legends
  lgd_mut <- Legend(title = "Mutation type", at = used, legend_gp = gpar(fill = class_cols[used]),
                    labels_gp = gpar(fontsize = 8), title_gp = gpar(fontsize = 9, fontface = "bold"))
  lgd_onc <- Legend(title = "PCGR oncogenicity", type = "points", labels = c("Oncogenic", "Likely_Oncogenic"),
                    pch = c(8, 24), size = unit(2, "mm"), legend_gp = gpar(col = "black", fill = "black"),
                    labels_gp = gpar(fontsize = 8), title_gp = gpar(fontsize = 9, fontface = "bold"))

  n_samples <- ncol(mat)
  ht <- Heatmap(
    mat, name = "Mutation type", col = c(class_cols, "None" = "grey80"), na_col = "grey80",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = genes, column_split = group, column_gap = unit(1.5, "mm"), column_title = NULL,
    show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 10),
    show_column_names = FALSE,
    width = unit(CELL_W_MM * n_samples, "mm"), height = unit(CELL_H_MM * length(genes), "mm"),
    rect_gp = gpar(type = "none"),
    top_annotation = top_ha, bottom_annotation = bottom_ha, right_annotation = right_ha,
    row_gap = unit(0.05, "mm"),
    show_heatmap_legend = FALSE,
    layer_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x, y, w * 0.995, h * 0.995, gp = gpar(fill = "grey80", col = "white", lwd = 0.1))
      m <- pindex(mat, i, j)
      ok <- m != "None"
      if (any(ok)) grid.rect(x[ok], y[ok], w[ok] * 0.96, h[ok] * 0.96, gp = gpar(fill = unname(class_cols[m[ok]]), col = NA))
      o <- pindex(onc_mat, i, j)
      sc <- ifelse(m == "Multi_Hit", "white", "black")
      star <- o == "Oncogenic"; tri <- o == "Likely_Oncogenic"
      if (any(star)) grid.points(x[star], y[star], pch = 8, size = unit(1.8, "mm"), gp = gpar(col = sc[star], lwd = 1))
      if (any(tri)) grid.points(x[tri], y[tri], pch = 24, size = unit(1.8, "mm"), gp = gpar(col = sc[tri], fill = sc[tri], lwd = 1))
    }
  )

  draw_fig <- function() {
    draw(ht, heatmap_legend_list = list(lgd_mut, lgd_onc),
         heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = FALSE,
         column_title = paste0(GENESET_LABEL, " | genes: ", length(genes), " | samples: ", n_samples,
                               " | WES VAF = ", MIN_VAF * 100, "%"),
         column_title_gp = gpar(fontsize = 11, fontface = "bold"),
         padding = unit(c(4, 4, 6, 4), "mm"))
  }
  fig_w <- (CELL_W_MM * n_samples + 3 * 1.5 + 35 + 42 + 22 + 85) / 25.4
  fig_h <- (CELL_H_MM * length(genes) + 30 + 4 + 7 * 4.5 + 25) / 25.4
  save_plot(file.path(out_dir, paste0(GENESET_LABEL, "_oncoplot_WES_", VAF_TAG)),
            width = max(fig_w, 11), height = max(fig_h, 6), draw_fn = draw_fig)
  graphics.off()
  invisible(NULL)
}

for (fam in names(gene_sets)) make_oncoplot(fam, gene_sets[[fam]])

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("All outputs are in:", out_dir, "\n")