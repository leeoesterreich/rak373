############################################################
## PURPOSE
##   Oncoplots of the top mutated genes in the 31 RNA-clustered pILC tumors,
##   split by the k = 3 NMF RNA clusters. Mutations: nonsilent, tumor VAF >= 3%.
##
## WHAT IS PLOTTED (for every set in GENE_SETS and every value in TOP_N_VALUES)
##   - Combined : a gene is altered in a tumor if WES, WGS or MSK-IMPACT calls it
##   - WES / WGS / IMPACT : one modality only (no mixing of assays)
##   - PER_CLUSTER_TOP_N : one extra oncoplot per cluster, genes ranked
##     WITHIN that cluster only
##
## HOW GENES ARE CHOSEN
##   Rank = number of altered tumors (ties: number of variant rows, then A-Z).
##   This ranking uses the tumors shown, so it is descriptive, not a test.
##
## HOW CELLS ARE DRAWN
##   Colour = variant class; Multi_Hit = more than one variant in that gene and
##   tumor. If several modalities call the same gene in a tumor, the modality
##   with the most variant rows is shown (no double counting).
##   Star = PCGR Oncogenic, triangle = PCGR Likely_Oncogenic (WES and WGS only;
##   MSK-IMPACT has no PCGR call).
##
## STATISTICS
##   Fisher exact test of "mutated vs not" across clusters, per gene, BH-adjusted
##   over the plotted genes. Exploratory: about 10 tumors per cluster.
##
## OUTPUTS (fixed file names, no date stamps) go to
##   .../RNASeq/Downstream/Oncoplot_by_RNA_cluster
############################################################


# ======================================================================
# 0. Libraries and seed
# ======================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

SEED <- 42
set.seed(SEED)
graphics.off()
gc()

# ======================================================================
# 1. Paths and parameters
# ======================================================================
# ============================================================
# Paths
# ============================================================
root <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC"

cluster_file <- file.path(root, "RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv")
# WES: combined MAF + PCGR VAF audit written by TopGenes + PCGR xlsx folders
wes_maf_file       <- file.path(root, "WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Tier_1_common_all_samples.combined.maf")
wes_vaf_audit_file <- file.path(root, "WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Full_MAF_nonsilent_variant_PCGR_VAF_match_audit_2026-09-23.tsv")
wes_pcgr_dir       <- file.path(root, "WES/Downstream/PCGR/With_PON/pcgr_outputs")

# WGS: PCGR VAF audit + PCGR xlsx folders (one folder per sample)
wgs_vaf_audit_file <- file.path(root, "WGS/PCGR/Downstream/Variant_classification/Top_Mutated_Genes_Cohort/Full_MAF_nonsilent_variant_PCGR_VAF_match_audit_2026-09-29.tsv")
wgs_pcgr_dir       <- file.path(root, "WGS/PCGR/pcgr_outputs")

# MSK-IMPACT (filtered table)
impact_file <- file.path(root, "MSK_IMPACT/pILC_IMPACT_filtered.tsv")

clinical_file <- file.path(root, "Clinical_correlation/pILC_ClinicoPathological_details.csv")

out_dir <- file.path(root, "RNASeq/Downstream/Oncoplot_by_RNA_cluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Tumors that were sequenced but have no row in a filtered mutation table
# (e.g. an IMPACT tumor with no retained variant). Use standardized IDs, e.g. "TP05-123".
extra_universe <- list(WES = character(), WGS = character(), IMPACT = character())

# ============================================================
# Parameters
# ============================================================
MIN_VAF <- 0.03
MIN_VAF_PERCENT <- 100 * MIN_VAF
VAF_TAG   <- paste0("VAFge", MIN_VAF_PERCENT, "pct")
VAF_LABEL <- paste0("VAF \u2265", MIN_VAF_PERCENT, "%")

TOP_N_VALUES <- c(50, 100)

# Separate oncoplot per RNA cluster, with genes ranked WITHIN that cluster only
PER_CLUSTER_TOP_N <- c(50)

# name -> modalities pooled in that plot
GENE_SETS <- list(
  Combined = c("WES", "WGS", "IMPACT"),
  WES      = "WES",
  WGS      = "WGS",
  IMPACT   = "IMPACT"
)

# When one tumor-gene is called by several modalities, the modality with the most
# variant rows is shown; ties go to the first modality in this order.
modality_priority <- c("WES", "WGS", "IMPACT")

nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Splice_Region", "Translation_Start_Site"
)

allow_sample_gene_fallback <- TRUE   # oncogenicity match on Sample+Gene if HGVS keys fail
SHOW_SAMPLE_NAMES <- TRUE
SAVE_PNG <- TRUE

input_files <- c(cluster_file = cluster_file, wes_maf_file = wes_maf_file,
                 wes_vaf_audit_file = wes_vaf_audit_file,
                 wgs_vaf_audit_file = wgs_vaf_audit_file,
                 impact_file = impact_file, clinical_file = clinical_file)
missing_inputs <- names(input_files)[!file.exists(input_files)]
if (length(missing_inputs) > 0) {
  stop("Input file(s) not found:\n",
       paste0(missing_inputs, ": ", input_files[missing_inputs], collapse = "\n"))
}

# ======================================================================
# 2. Colours
# ======================================================================
# ============================================================
# Colors (same as the earlier oncoplot script)
# ============================================================
grade_cols <- c(
  "Grade I" = "#A6CEE3", "Grade II" = "#33A02C", "Grade III" = "#1F78B4",
  "Cell type not determined" = "#E31A1C", "Unknown" = "#BDBDBD"
)

binary_receptor_cols <- c("Negative" = "#1F78B4", "Positive" = "#E31A1C", "Unknown" = "#BDBDBD")

base_mutation_colors <- c(
  "Missense_Mutation"       = "#1F77B4",
  "Nonsense_Mutation"       = "#D62728",
  "Frame_Shift_Del"         = "#9467BD",
  "Frame_Shift_Ins"         = "#8C564B",
  "In_Frame_Del"            = "#2CA02C",
  "In_Frame_Ins"            = "#17BECF",
  "Splice_Site"             = "#FF7F0E",
  "Splice_Region"           = "#FFBB78",
  "Translation_Start_Site"  = "#E377C2",
  "Nonstop_Mutation"        = "#BCBD22",
  "Multi_Hit"               = "#000000"
)

get_mutation_colors <- function(types) {
  types <- unique(types)
  types <- types[!is.na(types) & types != ""]
  known <- intersect(types, names(base_mutation_colors))
  missing <- setdiff(types, names(base_mutation_colors))
  cols <- base_mutation_colors[known]
  if (length(missing) > 0) {
    extra_cols <- hcl.colors(length(missing), palette = "Dark 3")
    names(extra_cols) <- missing
    cols <- c(cols, extra_cols)
  }
  cols[types]
}

modality_cols <- c("Sequenced" = "black", "Not sequenced" = "grey85")
cluster_palette <- c("#E7298A", "#66A61E", "#E6AB02", "#666666", "#A6761D")

# ======================================================================
# 3. Helper functions
# ======================================================================
# ============================================================
# Helpers
# ============================================================
clean_key <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", ".", "NA", "N/A", "nan", "NaN")] <- NA_character_
  toupper(x)
}

# One common sample ID for every modality, the RNA IDs and the clinical file:
#   TP13.1150_S40, TP13_1150, TP13-1150 -> TP13-1150
#   TP09-09T (IMPACT), TP09-09_T (WES), JC287-TP09-9 (WGS) -> TP09-9
std_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^JC287-", "", x)
  x <- sub("_S[0-9]+$", "", x)                  # RNA library suffix (_S40)
  x <- sub("(_T|-T|\\.T)$", "", x)              # tumor suffix
  x <- sub("^(TP[0-9]+)[._-]", "\\1-", x)       # TP13.1150 / TP13_1150 -> TP13-1150
  x <- sub("(?<=[0-9])T$", "", x, perl = TRUE)  # TP03-117T
  sub("-0+(?=[0-9])", "-", x, perl = TRUE)      # TP09-09 -> TP09-9
}

natural_sample_order <- function(samples) {
  samples <- unique(as.character(samples))
  tp <- suppressWarnings(as.integer(sub("^TP([0-9]+).*", "\\1", samples)))
  tp[is.na(tp)] <- 9999L
  rest <- suppressWarnings(as.integer(gsub("[^0-9]", "", sub("^TP[0-9]+[-_]*", "", samples))))
  rest[is.na(rest)] <- 9999L
  samples[order(tp, rest, samples)]
}

clean_hgv <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", ".", "NA", "N/A", "nan", "NaN")] <- NA_character_
  x <- sub(".*:", "", x)
  toupper(x)
}

normalize_oncogenicity <- function(x) {
  x <- toupper(gsub("[ -]+", "_", trimws(as.character(x))))
  out <- rep(NA_character_, length(x))
  out[x == "ONCOGENIC"] <- "Oncogenic"
  out[x == "LIKELY_ONCOGENIC"] <- "Likely_Oncogenic"
  out
}

normalize_colname <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

find_col <- function(cols, candidates) {
  hit <- match(normalize_colname(candidates), normalize_colname(cols))
  hit <- hit[!is.na(hit)]
  if (length(hit) > 0) cols[hit[1]] else NA_character_
}

get_col_or_na <- function(dt, candidates) {
  col <- find_col(colnames(dt), candidates)
  if (is.na(col)) rep(NA_character_, nrow(dt)) else as.character(dt[[col]])
}

safe_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

parse_percent_value <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NaN", "nan")] <- NA_character_
  is_less_than <- grepl("^<", x)
  val <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  val[is_less_than & !is.na(val)] <- val[is_less_than & !is.na(val)] / 2
  val
}

clean_receptor_status <- function(x) {
  x <- gsub("\\s+", " ", trimws(as.character(x)))
  out <- rep("Unknown", length(x))
  out[grepl("^Positive", x, ignore.case = TRUE)] <- "Positive"
  out[grepl("^Negative", x, ignore.case = TRUE)] <- "Negative"
  out[x %in% c("", "NA", "N/A", "NaN", "nan") | is.na(x)] <- "Unknown"
  out
}

clean_grade_status <- function(x) {
  x <- gsub("\\s+", " ", trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NaN", "nan") | is.na(x)] <- "Unknown"
  x
}

# ------------------------------------------------------------
# Variant-key helpers copied UNCHANGED from the TopGenes script.
# They must stay identical so Match_key agrees with the WES VAF audit table.
# ------------------------------------------------------------
vaf_clean_sample <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("(_T|-T|\\.T)$", "", x)
  x[x %in% c("", ".", "NA", "N/A")] <- NA_character_
  x
}

vaf_clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", ".", "NA", "N/A")] <- NA_character_
  x
}

edit_key <- function(chr, pos, ref, alt, maf_format = FALSE) {
  chr <- sub("^CHR", "", toupper(trimws(as.character(chr))))
  pos <- suppressWarnings(as.integer(pos))
  ref <- toupper(trimws(as.character(ref)))
  alt <- toupper(trimws(as.character(alt)))

  if (anyNA(c(chr, pos, ref, alt)) ||
      !nzchar(chr) || pos < 1L ||
      !grepl("^([ACGTN]+|-)$", ref) ||
      !grepl("^([ACGTN]+|-)$", alt)) {
    return(NA_character_)
  }

  if (maf_format && ref == "-") pos <- pos + 1L
  if (ref == "-") ref <- ""
  if (alt == "-") alt <- ""

  while (nzchar(ref) && nzchar(alt) &&
         substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))) {
    ref <- substr(ref, 1, nchar(ref) - 1L)
    alt <- substr(alt, 1, nchar(alt) - 1L)
  }

  while (nzchar(ref) && nzchar(alt) && substr(ref, 1, 1) == substr(alt, 1, 1)) {
    ref <- substring(ref, 2)
    alt <- substring(alt, 2)
    pos <- pos + 1L
  }

  if (!nzchar(ref) && !nzchar(alt)) return(NA_character_)

  paste(chr, pos, ifelse(nzchar(ref), ref, "-"), ifelse(nzchar(alt), alt, "-"), sep = ":")
}

# ======================================================================
# 4. RNA clusters
# ======================================================================
# ============================================================
# RNA clusters: the .gct is a membership table (Name, Description, membership.ordered)
# ============================================================
if (grepl("\\.csv$", cluster_file, ignore.case = TRUE)) {
  # CSV (Cluster_membership_k3_GPlabels.csv): one sample column + one cluster column
  cc <- fread(cluster_file)
  samp_col <- intersect(c("sample", "Sample", "Name", "name"), names(cc))
  samp_col <- if (length(samp_col)) samp_col[1] else names(cc)[1]
  clus_col <- setdiff(grep("cluster|membership", names(cc), ignore.case = TRUE, value = TRUE), samp_col)[1]
  if (is.na(clus_col)) stop("No cluster column found in ", cluster_file)
  message("Cluster CSV: sample column '", samp_col, "', cluster column '", clus_col, "'")
  cluster_dt <- data.table(
    RNA_ID = as.character(cc[[samp_col]]),
    Cluster_num = as.integer(gsub("[^0-9]", "", as.character(cc[[clus_col]])))
  )
} else {
  # GenePattern .gct membership table (Name, Description, membership.ordered)
  gct_raw <- fread(cluster_file, sep = "\t", skip = 2, header = TRUE, data.table = TRUE)
  memb_col <- setdiff(colnames(gct_raw), c("Name", "Description"))
  stopifnot(length(memb_col) == 1, "Name" %in% colnames(gct_raw))
  cluster_dt <- data.table(
    RNA_ID = as.character(gct_raw$Name),
    Cluster_num = as.integer(gct_raw[[memb_col]])
  )
}
cluster_dt[, Sample := std_id(RNA_ID)]
cluster_dt[, Cluster := paste0("Cluster ", Cluster_num)]

stopifnot(!anyNA(cluster_dt$Cluster_num), !anyDuplicated(cluster_dt$Sample))

cat("RNA-clustered samples:", nrow(cluster_dt), "\n")
print(cluster_dt[, .N, by = Cluster][order(Cluster)])

cluster_levels_all <- paste0("Cluster ", sort(unique(cluster_dt$Cluster_num)))
cluster_cols <- setNames(cluster_palette[seq_along(cluster_levels_all)], cluster_levels_all)

# ======================================================================
# 5. Clinical annotation
# ======================================================================
# ============================================================
# Clinical file (same parsing as the earlier oncoplot script)
# ============================================================
clinical_dt <- fread(clinical_file, sep = ",", header = TRUE, data.table = TRUE,
                     fill = TRUE, na.strings = c("", "NA", "NaN", "N/A"))

clinical_dt[, Sample_key := std_id(Sample)]
clinical_dt[, tumor_cellularity_pct := parse_percent_value(`tumor_cellularity (%)`)]
clinical_dt[, pleomorphic_component_cellularity_pct := parse_percent_value(
  `pleomorphic_invasive_compoent_cellularity (%)`)]
clinical_dt[, TILs_pct := parse_percent_value(`TILs (%)`)]
clinical_dt[, Grade_clean := clean_grade_status(`Grade/Differentiation-Desc`)]
clinical_dt[, HER2_status := clean_receptor_status(`BREAST: Her-2/Neu-Desc`)]
clinical_dt[, ER_status := clean_receptor_status(`Tumor Marker #1 (ER)`)]
clinical_dt[, PR_status := clean_receptor_status(`Tumor Marker #2 (PR)`)]

clinical_plot_dt <- unique(
  clinical_dt[, .(Sample_key, tumor_cellularity_pct, pleomorphic_component_cellularity_pct,
                  TILs_pct, Grade_clean, HER2_status, ER_status, PR_status)],
  by = "Sample_key"
)
cat("Clinical samples:", nrow(clinical_plot_dt), "\n")

# ======================================================================
# 6. Variants from WES, WGS and MSK-IMPACT (nonsilent, VAF >= cutoff)
# ======================================================================
# ============================================================
# Read the three modalities into ONE variant table
# (nonsilent, tumor VAF >= MIN_VAF). universe[[m]] = tumors sequenced by modality m
# BEFORE filtering, so tumors with no retained variant are still "sequenced".
# ============================================================
make_variant_table <- function(modality, sample, gene, cls, vaf_fraction,
                               hgvsp, hgvsp_short, hgvsc) {
  dt <- data.table(
    Modality = modality,
    Sample = std_id(sample),
    Gene = as.character(gene),
    Gene_key = clean_key(gene),
    Variant_Classification = as.character(cls),
    VAF_percent = round(100 * as.numeric(vaf_fraction), 2),
    HGVSp = hgvsp, HGVSp_Short = hgvsp_short, HGVSc = hgvsc
  )
  dt[, Sample_key := Sample]
  dt[, HGVSp_key := clean_hgv(HGVSp)]
  dt[, HGVSp_Short_key := clean_hgv(HGVSp_Short)]
  dt[, HGVSc_key := clean_hgv(HGVSc)]
  dt
}

# ---- WES: MAF + variant key + PCGR VAF audit (same logic as the earlier oncoplot script)
maf_full <- fread(wes_maf_file, sep = "\t", header = TRUE, data.table = TRUE,
                  na.strings = c("", ".", "NA", "N/A"))
req <- c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification",
         "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
if (length(setdiff(req, colnames(maf_full))) > 0) {
  stop("WES MAF lacks: ", paste(setdiff(req, colnames(maf_full)), collapse = ", "))
}

maf_full[, Variant_key := mapply(
  edit_key, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
  MoreArgs = list(maf_format = TRUE), USE.NAMES = FALSE
)]
maf_full[, Match_key := paste(vaf_clean_sample(Tumor_Sample_Barcode),
                              vaf_clean_gene(Hugo_Symbol), Variant_key, sep = "|")]

wes_audit <- fread(wes_vaf_audit_file, sep = "\t", header = TRUE, data.table = TRUE)
if (!all(c("Match_key", "VAF_TUMOR_fraction") %in% colnames(wes_audit))) {
  stop("WES VAF audit must contain Match_key and VAF_TUMOR_fraction.")
}
wes_vaf_by_key <- unique(
  wes_audit[, .(Match_key, VAF_TUMOR_fraction = suppressWarnings(as.numeric(VAF_TUMOR_fraction)))],
  by = "Match_key"
)
maf_full <- merge(maf_full, wes_vaf_by_key, by = "Match_key", all.x = TRUE, sort = FALSE)

if (!"HGVSp" %in% colnames(maf_full)) maf_full[, HGVSp := NA_character_]
if (!"HGVSp_Short" %in% colnames(maf_full)) maf_full[, HGVSp_Short := NA_character_]
if (!"HGVSc" %in% colnames(maf_full)) maf_full[, HGVSc := NA_character_]

wes_sub <- maf_full[
  Variant_Classification %in% nonsilent_classes &
    !is.na(vaf_clean_gene(Hugo_Symbol)) &
    !is.na(VAF_TUMOR_fraction) & VAF_TUMOR_fraction >= MIN_VAF
]
wes_var <- make_variant_table("WES", wes_sub$Tumor_Sample_Barcode, wes_sub$Hugo_Symbol,
                              wes_sub$Variant_Classification, wes_sub$VAF_TUMOR_fraction,
                              wes_sub$HGVSp, wes_sub$HGVSp_Short, wes_sub$HGVSc)

# ---- IMPACT (VAF is the TUMOR_MAF column, as in the cross-modality script)
imp_raw <- fread(impact_file, sep = "\t", na.strings = c("", "NA"))
need <- c("TUMOR_SAMPLE", "SYMBOL", "Variant_Classification", "TUMOR_MAF")
if (length(setdiff(need, colnames(imp_raw))) > 0) {
  stop("IMPACT file lacks: ", paste(setdiff(need, colnames(imp_raw)), collapse = ", "))
}
imp_raw[, TUMOR_MAF := suppressWarnings(as.numeric(TUMOR_MAF))]
imp_sub <- imp_raw[
  Variant_Classification %in% nonsilent_classes & !is.na(clean_key(SYMBOL)) &
    !is.na(TUMOR_MAF) & TUMOR_MAF >= MIN_VAF
]
impact_var <- make_variant_table(
  "IMPACT", imp_sub$TUMOR_SAMPLE, imp_sub$SYMBOL, imp_sub$Variant_Classification, imp_sub$TUMOR_MAF,
  get_col_or_na(imp_sub, c("HGVSp", "HGVSP")),
  get_col_or_na(imp_sub, c("HGVSp_Short", "HGVSP_SHORT")),
  get_col_or_na(imp_sub, c("HGVSc", "HGVSC"))
)

# ---- WGS (PCGR VAF audit: nonsilent + matched rows only, as in the cross-modality script)
wgs_raw <- fread(wgs_vaf_audit_file, sep = "\t", na.strings = c("", "NA"))
need <- c("Sample_key", "Hugo_Symbol", "Variant_Classification",
          "VAF_TUMOR_fraction", "matched", "is_nonsilent")
if (length(setdiff(need, colnames(wgs_raw))) > 0) {
  stop("WGS audit lacks: ", paste(setdiff(need, colnames(wgs_raw)), collapse = ", "))
}
wgs_raw[, VAF_TUMOR_fraction := suppressWarnings(as.numeric(VAF_TUMOR_fraction))]
wgs_sub <- wgs_raw[
  as.logical(is_nonsilent) & as.logical(matched) &
    !is.na(clean_key(Hugo_Symbol)) &
    !is.na(VAF_TUMOR_fraction) & VAF_TUMOR_fraction >= MIN_VAF
]
wgs_var <- make_variant_table(
  "WGS", wgs_sub$Sample_key, wgs_sub$Hugo_Symbol, wgs_sub$Variant_Classification,
  wgs_sub$VAF_TUMOR_fraction,
  get_col_or_na(wgs_sub, c("HGVSp", "HGVSP", "Protein_Change")),
  get_col_or_na(wgs_sub, c("HGVSp_Short", "HGVSP_SHORT")),
  get_col_or_na(wgs_sub, c("HGVSc", "HGVSC", "Coding_Change"))
)

# ---- Sequenced tumors per modality (before filtering)
wgs_folders <- std_id(list.dirs(wgs_pcgr_dir, recursive = FALSE, full.names = FALSE))
universe <- list(
  WES    = union(unique(std_id(maf_full$Tumor_Sample_Barcode)), extra_universe$WES),
  WGS    = union(union(unique(std_id(wgs_raw$Sample_key)), wgs_folders), extra_universe$WGS),
  IMPACT = union(unique(std_id(imp_raw$TUMOR_SAMPLE)), extra_universe$IMPACT)
)

variants_all <- rbindlist(list(wes_var, wgs_var, impact_var), use.names = TRUE)
variants_all[, variant_id := .I]

cat("\nTumors sequenced per modality:\n")
print(sapply(universe, length))
cat("\nVariants retained (nonsilent,", VAF_LABEL, "):\n")
print(variants_all[, .N, by = Modality])

# Display symbol per gene (WES spelling preferred, then WGS, then IMPACT)
gene_display <- variants_all[order(match(Modality, modality_priority)),
                             .(Gene = Gene[1]), by = Gene_key]
gene_display_vec <- setNames(gene_display$Gene, gene_display$Gene_key)

# ======================================================================
# 7. Sample-ID check across modalities
# ======================================================================
# ============================================================
# Check how the RNA-clustered tumors map onto each modality and the clinical file
# ============================================================
id_map <- copy(cluster_dt)
for (m in names(universe)) id_map[, (paste0("In_", m)) := Sample %in% universe[[m]]]
id_map[, In_clinical := Sample %in% clinical_plot_dt$Sample_key]
id_map[, N_modalities := In_WES + In_WGS + In_IMPACT]
setorder(id_map, Cluster_num, Sample)

cat("RNA-clustered tumors found per modality:\n")
print(colSums(id_map[, .(In_WES, In_WGS, In_IMPACT, In_clinical)]))
cat("\nRNA-clustered tumors with NO modality (check ID spelling):\n")
print(id_map[N_modalities == 0, .(RNA_ID, Sample, Cluster)])

# Sequenced tumors that are not in the RNA cluster file: near-duplicate IDs would show up here
for (m in names(universe)) {
  cat("\n", m, "tumors not in the RNA cluster file:\n", sep = "")
  print(natural_sample_order(setdiff(universe[[m]], cluster_dt$Sample)))
}

stopifnot(any(id_map$N_modalities > 0))

fwrite(id_map[, .(RNA_ID, Sample, Cluster, In_WES, In_WGS, In_IMPACT, In_clinical)],
       file.path(out_dir, "RNA_cluster_sample_ID_mapping_and_modalities.tsv"),
       sep = "\t", quote = FALSE, na = "")

# ======================================================================
# 8. PCGR oncogenicity (WES + WGS)
# ======================================================================
# ============================================================
# PCGR oncogenicity (WES + WGS only; IMPACT has none) and matching to variants
# ============================================================
read_pcgr_onc <- function(pcgr_dir, modality) {
  files <- list.files(pcgr_dir, pattern = "\\.pcgr\\.grch38\\.xlsx$",
                      recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    files <- list.files(pcgr_dir, pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)
  }
  cat(modality, "PCGR xlsx files found:", length(files), "\n")

  out <- list()
  for (f in files) {
    sheets <- tryCatch(excel_sheets(f), error = function(e) character())
    if (!"SOMATIC_SNV_INDEL" %in% sheets) {
      message("[SKIP] SOMATIC_SNV_INDEL sheet not found: ", f)
      next
    }
    dt <- tryCatch(as.data.table(read_excel(f, sheet = "SOMATIC_SNV_INDEL")),
                   error = function(e) { message("[SKIP] Could not read: ", f); NULL })
    if (is.null(dt) || nrow(dt) == 0) next

    gene_col <- find_col(colnames(dt), c("Hugo_Symbol", "HUGO_SYMBOL", "Gene", "Gene_Symbol",
                                         "SYMBOL", "GENE_NAME"))
    onc_col <- find_col(colnames(dt), c("ONCOGENICITY", "Oncogenicity"))
    if (is.na(gene_col) || is.na(onc_col)) {
      message("[SKIP] Could not find Gene or ONCOGENICITY column in: ", f)
      next
    }

    tmp <- data.table(
      Modality = modality,
      Sample_key = std_id(basename(dirname(f))),
      Gene_key = clean_key(dt[[gene_col]]),
      HGVSp_key = clean_hgv(get_col_or_na(dt, c("HGVSp", "HGVSP", "Protein_Change", "ProteinChange"))),
      HGVSp_Short_key = clean_hgv(get_col_or_na(dt, c("HGVSp_Short", "HGVSP_SHORT"))),
      HGVSc_key = clean_hgv(get_col_or_na(dt, c("HGVSc", "HGVSC", "Coding_Change", "CodingChange"))),
      ONCOGENICITY = normalize_oncogenicity(dt[[onc_col]])
    )
    tmp <- tmp[!is.na(ONCOGENICITY)]
    if (nrow(tmp) > 0) out[[length(out) + 1]] <- tmp
  }

  if (length(out) == 0) {
    warning("No Oncogenic/Likely_Oncogenic PCGR records found for ", modality)
    return(data.table())
  }
  d <- rbindlist(out, fill = TRUE)
  d[, onc_priority := fifelse(ONCOGENICITY == "Oncogenic", 2L, 1L)]
  d
}

pcgr_dt <- rbindlist(list(read_pcgr_onc(wes_pcgr_dir, "WES"),
                          read_pcgr_onc(wgs_pcgr_dir, "WGS")), fill = TRUE)

add_pcgr_match <- function(annot, pcgr, keys, method_name) {
  if (nrow(pcgr) == 0) return(annot)

  pcgr_sub <- copy(pcgr)
  for (k in keys) pcgr_sub <- pcgr_sub[!is.na(get(k)) & get(k) != ""]
  if (nrow(pcgr_sub) == 0) return(annot)

  pcgr_best <- pcgr_sub[
    order(-onc_priority),
    .(ONCOGENICITY_MATCH = first(ONCOGENICITY), PCGR_match_priority = first(onc_priority)),
    by = keys
  ]

  unmatched <- annot[is.na(ONCOGENICITY)]
  for (k in keys) unmatched <- unmatched[!is.na(get(k)) & get(k) != ""]
  if (nrow(unmatched) == 0) return(annot)

  tmp <- merge(unmatched[, c("variant_id", keys), with = FALSE], pcgr_best,
               by = keys, all.x = FALSE, all.y = FALSE, allow.cartesian = TRUE)
  if (nrow(tmp) == 0) return(annot)

  tmp <- tmp[order(-PCGR_match_priority), .SD[1], by = variant_id]

  annot[tmp, on = "variant_id",
        `:=`(ONCOGENICITY = i.ONCOGENICITY_MATCH, onc_match_method = method_name)]
  annot
}

onc_annot <- variants_all[, .(variant_id, Modality, Sample_key, Gene_key,
                              HGVSp_key, HGVSp_Short_key, HGVSc_key)]
onc_annot[, ONCOGENICITY := NA_character_]
onc_annot[, onc_match_method := NA_character_]

if (nrow(pcgr_dt) > 0) {
  onc_annot <- add_pcgr_match(onc_annot, pcgr_dt,
                              c("Modality", "Sample_key", "Gene_key", "HGVSp_key"),
                              "Sample+Gene+HGVSp")
  onc_annot <- add_pcgr_match(onc_annot, pcgr_dt,
                              c("Modality", "Sample_key", "Gene_key", "HGVSp_Short_key"),
                              "Sample+Gene+HGVSp_Short")
  onc_annot <- add_pcgr_match(onc_annot, pcgr_dt,
                              c("Modality", "Sample_key", "Gene_key", "HGVSc_key"),
                              "Sample+Gene+HGVSc")
  if (allow_sample_gene_fallback) {
    onc_annot <- add_pcgr_match(onc_annot, pcgr_dt,
                                c("Modality", "Sample_key", "Gene_key"),
                                "Sample+Gene_fallback")
  }
}

idx <- match(variants_all$variant_id, onc_annot$variant_id)
variants_all[, ONCOGENICITY := onc_annot$ONCOGENICITY[idx]]
variants_all[, onc_match_method := onc_annot$onc_match_method[idx]]

cat("\nVariants matched to PCGR Oncogenic / Likely_Oncogenic (by modality and match method):\n")
print(variants_all[!is.na(ONCOGENICITY), .N,
                   by = .(Modality, ONCOGENICITY, onc_match_method)][order(Modality, ONCOGENICITY)])

# ======================================================================
# 9. Oncoplot function
# ======================================================================
# ============================================================
# Plot function: one oncoplot split by RNA cluster
# ============================================================
make_cluster_oncoplot <- function(set_name, modalities, top_n, cluster_name = NULL) {

  # cluster_name = NULL : all clusters side by side (genes ranked across all shown tumors)
  # cluster_name = "Cluster 2" : only that cluster (genes ranked within that cluster)
  scope_label <- if (is.null(cluster_name)) "by_RNA_cluster" else safe_file_name(cluster_name)

  cat("============================================================\n")
  cat("Set:", set_name, "| modalities:", paste(modalities, collapse = "+"),
      "| top", top_n, "genes |", if (is.null(cluster_name)) "all clusters" else cluster_name, "\n")

  # ---- Samples: RNA-clustered tumors sequenced by >= 1 modality of this set,
  #      ordered by cluster, then by TP number within each cluster
  seq_in_set <- unique(unlist(universe[modalities]))
  keep <- cluster_dt[Sample %in% seq_in_set]
  if (!is.null(cluster_name)) keep <- keep[Cluster == cluster_name]
  if (nrow(keep) == 0) {
    cat("No RNA-clustered tumors for this set. Skipping.\n\n")
    return(NULL)
  }

  cluster_levels_used <- paste0("Cluster ", sort(unique(keep$Cluster_num)))
  plot_samples <- unlist(lapply(sort(unique(keep$Cluster_num)), function(cn) {
    natural_sample_order(keep[Cluster_num == cn, Sample])
  }))
  sample_cluster <- cluster_dt$Cluster[match(plot_samples, cluster_dt$Sample)]
  cat("Samples shown:", length(plot_samples), "\n")

  # ---- Variants of this set in the shown samples; rank genes by altered samples
  var_set <- variants_all[Modality %in% modalities & Sample %in% plot_samples]
  if (nrow(var_set) == 0) {
    cat("No variants for this set. Skipping.\n\n")
    return(NULL)
  }

  gene_summary <- var_set[, .(n_altered_samples = uniqueN(Sample), n_variant_rows = .N),
                          by = Gene_key]
  setorder(gene_summary, -n_altered_samples, -n_variant_rows, Gene_key)

  plot_gene_keys <- head(gene_summary$Gene_key, top_n)
  cutoff_n <- gene_summary$n_altered_samples[length(plot_gene_keys)]
  n_tied_out <- sum(gene_summary$n_altered_samples == cutoff_n) -
    sum(gene_summary$n_altered_samples[seq_along(plot_gene_keys)] == cutoff_n)
  cat("Genes with >= 1 nonsilent call:", nrow(gene_summary), "| plotted:", length(plot_gene_keys),
      "| lowest plotted frequency:", cutoff_n, "sample(s);",
      n_tied_out, "gene(s) with the same frequency are not plotted (cut by variant rows, then A-Z)\n")

  plot_genes <- unname(gene_display_vec[plot_gene_keys])

  # ---- One cell per sample x gene. If several modalities call the same gene in a
  #      tumor, the modality with the most variant rows is shown (no double counting).
  v_top <- var_set[Gene_key %in% plot_gene_keys]
  v_top[, mod_priority := match(Modality, modality_priority)]

  by_mod <- v_top[, .(n_variants = .N, class_first = Variant_Classification[1]),
                  by = .(Sample, Gene_key, Modality, mod_priority)]
  setorder(by_mod, Sample, Gene_key, -n_variants, mod_priority)
  rep_dt <- by_mod[, .SD[1], by = .(Sample, Gene_key)]

  v_rep <- merge(v_top, rep_dt[, .(Sample, Gene_key, Modality)],
                 by = c("Sample", "Gene_key", "Modality"))

  cell_info <- v_top[, {
    onc_vals <- ONCOGENICITY[!is.na(ONCOGENICITY)]
    cell_onc <- ""
    if ("Oncogenic" %in% onc_vals) {
      cell_onc <- "Oncogenic"
    } else if ("Likely_Oncogenic" %in% onc_vals) {
      cell_onc <- "Likely_Oncogenic"
    }
    list(
      Modalities_detected = paste(unique(Modality[order(mod_priority)]), collapse = "+"),
      Variant_Classification_List = paste(unique(Variant_Classification), collapse = ";"),
      HGVSp_List = paste(unique(na.omit(HGVSp)), collapse = ";"),
      VAF_percent_by_modality = paste(paste0(Modality, ":", VAF_percent), collapse = ";"),
      Cell_ONCOGENICITY = cell_onc
    )
  }, by = .(Sample, Gene_key)]

  cell_dt <- merge(
    rep_dt[, .(Sample, Gene_key, n_variants,
               Display = fifelse(n_variants > 1L, "Multi_Hit", class_first))],
    cell_info, by = c("Sample", "Gene_key")
  )
  cell_dt[, Gene := unname(gene_display_vec[Gene_key])]
  cell_dt[, Cluster := cluster_dt$Cluster[match(Sample, cluster_dt$Sample)]]

  # ---- Matrices
  mut_mat <- matrix("", nrow = length(plot_gene_keys), ncol = length(plot_samples),
                    dimnames = list(plot_genes, plot_samples))
  onc_mat <- mut_mat
  cell_idx <- cbind(match(cell_dt$Gene_key, plot_gene_keys), match(cell_dt$Sample, plot_samples))
  mut_mat[cell_idx] <- cell_dt$Display
  onc_mat[cell_idx] <- cell_dt$Cell_ONCOGENICITY

  plot_mat <- mut_mat
  plot_mat[plot_mat == ""] <- NA_character_

  # ---- Top stacked bar: variant rows per sample by class (shown modality only)
  cls_levels <- sort(unique(v_rep$Variant_Classification))
  sample_class_mat <- matrix(0, nrow = length(plot_samples), ncol = length(cls_levels),
                             dimnames = list(plot_samples, cls_levels))
  tmp <- v_rep[, .N, by = .(Sample, Variant_Classification)]
  sample_class_mat[cbind(match(tmp$Sample, plot_samples),
                         match(tmp$Variant_Classification, cls_levels))] <- tmp$N

  # ---- Right stacked bar: altered samples per gene by class
  disp_levels <- sort(unique(cell_dt$Display))
  gene_class_mat <- matrix(0, nrow = length(plot_gene_keys), ncol = length(disp_levels),
                           dimnames = list(plot_genes, disp_levels))
  tmp <- cell_dt[, .N, by = .(Gene_key, Display)]
  gene_class_mat[cbind(match(tmp$Gene_key, plot_gene_keys), match(tmp$Display, disp_levels))] <- tmp$N

  gene_mutated_counts <- rowSums(mut_mat != "")
  gene_mutated_pct <- round(100 * gene_mutated_counts / ncol(mut_mat))

  mutation_colors <- get_mutation_colors(c(na.omit(as.vector(plot_mat)), cls_levels, disp_levels))
  sample_bar_colors <- mutation_colors[colnames(sample_class_mat)]
  sample_bar_colors[is.na(sample_bar_colors)] <- "grey50"
  row_bar_colors <- mutation_colors[colnames(gene_class_mat)]
  row_bar_colors[is.na(row_bar_colors)] <- "grey50"

  top_bar_ha <- HeatmapAnnotation(
    `Mutations per sample` = anno_barplot(
      sample_class_mat,
      gp = gpar(fill = sample_bar_colors, col = NA),
      border = FALSE, bar_width = 0.97, height = unit(2.2, "cm"),
      axis_param = list(gp = gpar(fontsize = 8))
    ),
    annotation_name_gp = gpar(fontsize = 10)
  )

  right_bar_ha <- rowAnnotation(
    `Altered samples by class` = anno_barplot(
      gene_class_mat,
      gp = gpar(fill = row_bar_colors, col = NA),
      border = FALSE, bar_width = 0.90, axis = TRUE,
      axis_param = list(gp = gpar(fontsize = 7), side = "bottom")
    ),
    `% altered` = anno_text(
      paste0(gene_mutated_pct, "%"),
      gp = gpar(fontsize = 10), just = "left", location = 0
    ),
    width = unit(4.2, "cm"),
    show_annotation_name = FALSE
  )

  # ---- Bottom annotations: RNA cluster, modality sequenced, clinical
  clin <- merge(data.table(Sample = plot_samples), clinical_plot_dt,
                by.x = "Sample", by.y = "Sample_key", all.x = TRUE, sort = FALSE)
  clin <- clin[match(plot_samples, Sample)]

  anno_df <- data.frame(
    Cluster = sample_cluster,
    `WES seq.`    = ifelse(plot_samples %in% universe$WES,    "Sequenced", "Not sequenced"),
    `WGS seq.`    = ifelse(plot_samples %in% universe$WGS,    "Sequenced", "Not sequenced"),
    `IMPACT seq.` = ifelse(plot_samples %in% universe$IMPACT, "Sequenced", "Not sequenced"),
    `Tumor cell. (%)` = clin$tumor_cellularity_pct,
    `Pleomorphic (%)` = clin$pleomorphic_component_cellularity_pct,
    `TILs (%)` = clin$TILs_pct,
    Grade = clin$Grade_clean,
    HER2 = clin$HER2_status,
    ER = clin$ER_status,
    PR = clin$PR_status,
    check.names = FALSE
  )
  rownames(anno_df) <- plot_samples

  grade_cols_local <- grade_cols
  extra_grades <- setdiff(unique(na.omit(anno_df$Grade)), names(grade_cols))
  if (length(extra_grades) > 0) {
    eg <- hcl.colors(length(extra_grades), palette = "Dark 3")
    names(eg) <- extra_grades
    grade_cols_local <- c(grade_cols_local, eg)
  }

  anno_cols <- list(
    Cluster = cluster_cols[cluster_levels_used],
    `WES seq.` = modality_cols,
    `WGS seq.` = modality_cols,
    `IMPACT seq.` = modality_cols,
    `Tumor cell. (%)` = colorRamp2(c(0, 50, 100), c("#F7FBFF", "#6BAED6", "#08306B")),
    `Pleomorphic (%)` = colorRamp2(c(0, 50, 100), c("#F7FCF5", "#74C476", "#00441B")),
    `TILs (%)` = colorRamp2(c(0, 30, 100), c("#FFF5F0", "#FB6A4A", "#67000D")),
    Grade = grade_cols_local,
    HER2 = binary_receptor_cols,
    ER = binary_receptor_cols,
    PR = binary_receptor_cols
  )

  # one shared "Sequenced" legend for the three modality tracks
  anno_names <- colnames(anno_df)
  legend_flags <- c(FALSE, !(anno_names %in% c("WGS seq.", "IMPACT seq.")))

  bottom_ha <- HeatmapAnnotation(
    ` ` = anno_empty(height = unit(1, "mm"), border = FALSE),
    df = anno_df,
    col = anno_cols,
    na_col = "grey90",
    which = "column",
    annotation_name_side = "left",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(0.30, "cm"),
    gap = unit(c(3.5, rep(1.1, ncol(anno_df) - 1)), "mm"),
    show_annotation_name = setNames(c(FALSE, rep(TRUE, ncol(anno_df))), c(" ", anno_names)),
    show_legend = legend_flags,
    annotation_legend_param = list(
      Cluster = list(title = "RNA cluster"),
      `WES seq.` = list(title = "Sequenced")
    )
  )

  # ---- Layout
  col_split <- factor(sample_cluster, levels = cluster_levels_used)
  slice_titles <- paste0(cluster_levels_used, " (n = ",
                         as.integer(table(col_split)[cluster_levels_used]), ")")

  cell_size_mm <- 4.0
  body_w_mm <- max(55, cell_size_mm * ncol(mut_mat))
  body_h_mm <- max(38, cell_size_mm * nrow(mut_mat))
  plot_width  <- max(13, min(26, body_w_mm / 25.4 + 10.5))
  plot_height <- max(8, min(26, body_h_mm / 25.4 + 7))

  onc_lgd <- Legend(
    title = "PCGR oncogenicity", type = "points",
    labels = c("Oncogenic", "Likely_Oncogenic"), pch = c(8, 24),
    size = unit(c(2.0, 2.0), "mm"), legend_gp = gpar(col = "black", fill = "black")
  )
  extra_legends <- if (any(onc_mat != "")) list(onc_lgd) else list()

  ht <- Heatmap(
    plot_mat,
    name = "Mutation type",
    col = mutation_colors,
    na_col = "grey80",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = seq_len(nrow(plot_mat)),
    column_order = seq_len(ncol(plot_mat)),
    column_split = col_split,
    show_row_names = TRUE,
    show_column_names = SHOW_SAMPLE_NAMES,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 7),
    column_title = slice_titles,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    width = unit(body_w_mm, "mm"),
    height = unit(body_h_mm, "mm"),
    rect_gp = gpar(type = "none"),
    top_annotation = top_bar_ha,
    bottom_annotation = bottom_ha,
    right_annotation = right_bar_ha,
    column_gap = unit(2, "mm"),
    heatmap_legend_param = list(
      title = "Mutation type",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ),
    layer_fun = function(j, i, x, y, width, height, fill) {

      grid.rect(x = x, y = y, width = width * 0.995, height = height * 0.995,
                gp = gpar(fill = "grey80", col = "white", lwd = 0.10))

      for (kk in seq_along(x)) {
        this_mut <- mut_mat[i[kk], j[kk]]
        if (!is.na(this_mut) && this_mut != "") {
          grid.rect(x = x[kk], y = y[kk], width = width[kk] * 0.96, height = height[kk] * 0.96,
                    gp = gpar(fill = mutation_colors[this_mut], col = NA))
        }
      }

      for (kk in seq_along(x)) {
        this_mut <- mut_mat[i[kk], j[kk]]
        this_onc <- onc_mat[i[kk], j[kk]]
        symbol_col <- if (!is.na(this_mut) && this_mut == "Multi_Hit") "white" else "black"

        if (!is.na(this_onc) && this_onc == "Oncogenic") {
          grid.points(x = x[kk], y = y[kk], pch = 8, size = unit(1.8, "mm"),
                      gp = gpar(col = symbol_col, fill = symbol_col, lwd = 1.0))
        } else if (!is.na(this_onc) && this_onc == "Likely_Oncogenic") {
          grid.points(x = x[kk], y = y[kk], pch = 24, size = unit(1.8, "mm"),
                      gp = gpar(col = symbol_col, fill = symbol_col, lwd = 1.0))
        }
      }
    }
  )

  overall_title <- paste0("Top", top_n, " genes | ", set_name,
                          " (", paste(modalities, collapse = " + "), ")",
                          if (is.null(cluster_name)) "" else paste0(" | ", cluster_name, " only"),
                          " | genes: ", nrow(mut_mat),
                          " | samples: ", ncol(mut_mat), " | ", VAF_LABEL)

  draw_fn <- function() {
    draw(ht,
         column_title = overall_title,
         column_title_gp = gpar(fontsize = 11, fontface = "bold"),
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         heatmap_legend_list = extra_legends,
         merge_legends = FALSE,
         padding = unit(c(4, 28, 6, 4), "mm"))
  }

  base_name <- safe_file_name(paste0("Top", top_n, "_Genes_", set_name, "_", scope_label, "_", VAF_TAG))
  pdf_file <- file.path(out_dir, paste0(base_name, "_oncoplot_PCGR_oncogenicity.pdf"))
  png_file <- file.path(out_dir, paste0(base_name, "_oncoplot_PCGR_oncogenicity.png"))

  pdf(pdf_file, width = plot_width, height = plot_height)
  draw_fn()
  dev.off()

  if (isTRUE(SAVE_PNG)) {
    png(png_file, width = plot_width, height = plot_height, units = "in", res = 300)
    draw_fn()
    dev.off()
  } else {
    png_file <- NA_character_
  }

  # ---- Stats table: mutated tumors per gene and cluster, Fisher exact across clusters
  #      (exploratory: few tumors per cluster), BH-adjusted over the plotted genes
  set.seed(SEED)
  cl_f <- factor(sample_cluster, levels = cluster_levels_used)
  stats_dt <- rbindlist(lapply(seq_len(nrow(mut_mat)), function(r) {
    mut <- mut_mat[r, ] != ""
    tab <- table(factor(mut, levels = c(FALSE, TRUE)), cl_f)
    p <- if (nlevels(cl_f) < 2) NA_real_ else tryCatch(
      fisher.test(tab)$p.value,
      error = function(e) tryCatch(
        fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value,
        error = function(e2) NA_real_
      )
    )
    per_cluster <- as.list(tapply(mut, cl_f, sum))
    names(per_cluster) <- paste0("N_mutated_", gsub(" ", "_", names(per_cluster)))
    data.table(Gene = rownames(mut_mat)[r], N_mutated = sum(mut),
               Pct_mutated = round(100 * mean(mut), 1),
               as.data.table(per_cluster), Fisher_P = p)
  }))
  stats_dt[, FDR_BH := p.adjust(Fisher_P, method = "BH")]
  if (nlevels(cl_f) < 2) stats_dt[, c("Fisher_P", "FDR_BH") := NULL]   # no across-cluster test for one cluster

  size_row <- as.list(table(cl_f))
  names(size_row) <- paste0("N_mutated_", gsub(" ", "_", names(size_row)))
  cat("Tumors per cluster in this plot:", paste(names(size_row), unlist(size_row), sep = "=", collapse = ", "), "\n")

  fwrite(stats_dt, file.path(out_dir, paste0(base_name, "_gene_by_cluster_stats.tsv")),
         sep = "\t", quote = FALSE, na = "")

  setcolorder(cell_dt, c("Sample", "Cluster", "Gene", "Display", "n_variants",
                         "Modalities_detected", "Cell_ONCOGENICITY"))
  fwrite(cell_dt[order(match(Sample, plot_samples), match(Gene_key, plot_gene_keys))],
         file.path(out_dir, paste0(base_name, "_cell_level_table.tsv")),
         sep = "\t", quote = FALSE, na = "")

  cat("Saved:\n  ", pdf_file, "\n")
  if (!is.na(png_file)) cat("  ", png_file, "\n")
  cat("\n")

  graphics.off()
  gc()

  data.table(Set = set_name, Top_N = top_n, Genes_plotted = nrow(mut_mat),
             Samples_shown = ncol(mut_mat), Cells_mutated = sum(mut_mat != ""),
             Cells_Oncogenic = sum(onc_mat == "Oncogenic"),
             Cells_Likely_Oncogenic = sum(onc_mat == "Likely_Oncogenic"))
}

# ======================================================================
# 10. Run every gene set x top-N combination
# ======================================================================
# ============================================================
# Run every gene set x top-N combination
# ============================================================
graphics.off()
gc()

summaries <- list()

for (set_name in names(GENE_SETS)) {
  for (top_n in TOP_N_VALUES) {
    res <- tryCatch(
      make_cluster_oncoplot(set_name, GENE_SETS[[set_name]], top_n),
      error = function(e) {
        message("[ERROR] ", set_name, " top ", top_n, ": ", e$message)
        graphics.off()
        NULL
      }
    )
    if (!is.null(res)) summaries[[length(summaries) + 1]] <- res
  }
}

# Separate oncoplot for each cluster alone (top genes ranked within that cluster)
for (set_name in names(GENE_SETS)) {
  for (cl in cluster_levels_all) {
    for (top_n in PER_CLUSTER_TOP_N) {
      res <- tryCatch(
        make_cluster_oncoplot(set_name, GENE_SETS[[set_name]], top_n, cluster_name = cl),
        error = function(e) {
          message("[ERROR] ", set_name, " | ", cl, " | top ", top_n, ": ", e$message)
          graphics.off()
          NULL
        }
      )
      if (!is.null(res)) {
        res[, Cluster := cl]
        summaries[[length(summaries) + 1]] <- res
      }
    }
  }
}

if (length(summaries) > 0) print(rbindlist(summaries, fill = TRUE))

cat("\nAll done.\nOutputs in:\n", out_dir, "\n")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))