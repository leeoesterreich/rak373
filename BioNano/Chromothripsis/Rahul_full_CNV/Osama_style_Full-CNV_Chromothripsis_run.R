#!/usr/bin/env Rscript

# ============================================================
# Osama-style Chromothripsis / ShatterSeek reproduction 
# ============================================================
# Reproduces the ACTIVE Osama-style ShatterSeek workflow.
#
# SINGLE ANALYTICAL CHANGE RELATIVE TO THE CURRENT OSAMA-STYLE RUN:
#   CNV input is switched from Rick's filtered *_Annotated_CNV.txt files
#   to the fuller Rahul/CNV_Annotation/<SampleID>.txt files.
#
# Both originate from the same Bionano Solve CNV annotation
# pipeline, but the Rahul file contains the fuller CNV call set.
#
# Everything else is retained:
# 1) Same 62 original Bionano SMAP files from Rick's folder.
# 2) NON-DEDUPLICATED SV_df (same as Osama for ShatterSeek).
# 3) Commented Access 1.7 SV filters remain inactive.
# 4) Same SV type conversion and High/Low thresholds.
# 5) Same technical GRanges compatibility workaround.
#
# Validation:
#   Using Rick SMAP files + Rahul fuller CNV files reproduced
#   Osama's original chromothripsis heatmap.
# ============================================================

`%notin%` <- Negate(`%in%`)

suppressPackageStartupMessages({
  library(ShatterSeek)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(reshape2)
  library(ComplexHeatmap)
  library(grid)
})

# ============================================================
# PATHS
# ============================================================
SV_DIR <- "/ix1/alee/LO_LAB/Personal/Rick/Bionano/20231216_Liane_Osama_Bionano_DNA_SV"
CNV_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation"

# Separate output folder prevents mixing this run with the previous
# Rick-filtered-CNV Osama-style outputs.
OUT_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Chromothripsis/Osama_Script/Rahul_full_CNV"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("SCRIPT VERSION: V5_OSAMA_STYLE_RAHUL_FULL_CNV\n")
cat("SV_DIR :", SV_DIR, "\n")
cat("CNV_DIR:", CNV_DIR, "\n")
cat("OUT_DIR :", OUT_DIR, "\n\n")
cat("ANALYSIS CHANGE: CNV input only -> Rahul/CNV_Annotation full files\n\n")

# ============================================================
# 1. FIND SMAP + CNV FILES
# ============================================================
svFilesPaths <- list.files(
  SV_DIR,
  pattern = "\\.smap$",
  full.names = TRUE
)

cnvFilesPaths <- list.files(
  CNV_DIR,
  pattern = "\\.txt$",
  full.names = TRUE
)

cat("Number of SMAP files:", length(svFilesPaths), "\n")
cat("Number of CNV files :", length(cnvFilesPaths), "\n\n")

if (length(svFilesPaths) != 62) {
  stop("Expected 62 .smap files but found ", length(svFilesPaths))
}

# ============================================================
# 2. PAIR SMAPs WITH THE FULLER RAHUL CNV FILES
#
# Only CNV source is changed relative to Current Osama-style:
#   SVs  = same raw SMAP files from Rick folder
#   CNVs = Rahul/CNV_Annotation/<SampleID>.txt
# ============================================================

# Same sample-ID cleaning rule Osama used, adapted to basename().
SampleIDs <- gsub(
  x = basename(svFilesPaths),
  pattern = "Annotated_SV.smap|_12_16_|_12_17_|2023_[0-9]+_[0-9]+_[0-9]+_|_|_cr|_s|_b|_A",
  replacement = ""
)

# Rahul CNV files are named directly by SampleID, e.g. TP19-M6.txt.
cnv_ids <- sub("\\.txt$", "", basename(cnvFilesPaths))

file_map <- data.frame(
  SampleID = SampleIDs,
  smap = svFilesPaths,
  stringsAsFactors = FALSE
)

file_map$cnv <- cnvFilesPaths[
  match(file_map$SampleID, cnv_ids)
]

cat("Unique derived SampleIDs:", length(unique(file_map$SampleID)), "\n")
cat("SMAPs without matched Rahul CNV:", sum(is.na(file_map$cnv)), "\n\n")

if (anyDuplicated(file_map$SampleID)) {
  dup_ids <- unique(
    file_map$SampleID[
      duplicated(file_map$SampleID) |
      duplicated(file_map$SampleID, fromLast = TRUE)
    ]
  )
  stop(
    "Duplicate SampleIDs were derived from SMAP filenames: ",
    paste(dup_ids, collapse = ", ")
  )
}

if (anyDuplicated(cnv_ids)) {
  dup_ids <- unique(
    cnv_ids[
      duplicated(cnv_ids) |
      duplicated(cnv_ids, fromLast = TRUE)
    ]
  )
  stop(
    "Duplicate Rahul CNV SampleIDs were found: ",
    paste(dup_ids, collapse = ", ")
  )
}

if (any(is.na(file_map$cnv))) {
  print(file_map[is.na(file_map$cnv), , drop = FALSE])
  stop("At least one SMAP sample could not be paired with a Rahul CNV file.")
}

write.table(
  file_map,
  file = file.path(OUT_DIR, "SMAP_CNV_file_map.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# CNV input QC BEFORE ShatterSeek
# ------------------------------------------------------------

cnv_qc_list <- lapply(seq_len(nrow(file_map)), function(i) {

  f <- file_map$cnv[i]

  x <- read.delim(
    f,
    skip = 3,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c(
    "Chromosome",
    "Start",
    "End",
    "CopyNumber"
  )

  data.frame(
    SampleID = file_map$SampleID[i],
    CNV_file = basename(f),
    N_CNV_calls = nrow(x),
    Required_columns_present = all(required %in% names(x)),
    stringsAsFactors = FALSE
  )
})

cnv_qc <- do.call(rbind, cnv_qc_list)

write.table(
  cnv_qc,
  file = file.path(OUT_DIR, "CNV_input_QC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("CNV files with 0 calls:", sum(cnv_qc$N_CNV_calls == 0), "\n")
cat(
  "CNV files missing required columns:",
  sum(!cnv_qc$Required_columns_present),
  "\n\n"
)

if (any(!cnv_qc$Required_columns_present)) {
  print(cnv_qc[!cnv_qc$Required_columns_present, , drop = FALSE])
  stop("At least one Rahul CNV file is missing required ShatterSeek CNV columns.")
}

if (any(cnv_qc$N_CNV_calls == 0)) {
  print(cnv_qc[cnv_qc$N_CNV_calls == 0, , drop = FALSE])
  stop(
    "At least one Rahul CNV file contains zero CNV rows. ",
    "Stop here so the Osama-style 62-sample comparison remains clean."
  )
}

# ============================================================
# 3. LOAD SV FILES — SAME AS OSAMA
# ============================================================
SV_df <- lapply(svFilesPaths, function(x) {
  temp <- read.delim(x, skip = 6)
  temp <- temp[-1, ]       # remove #f row, same as Osama
  temp
})

names(SV_df) <- SampleIDs

cat("First derived sample IDs:\n")
print(head(names(SV_df), 10))
cat("\n")

# ============================================================
# 4. COMMON COLUMNS — SAME TWO SAMPLES OSAMA USED
# ============================================================
required_common_samples <- c("TP11-896", "TP21-M598")

if (!all(required_common_samples %in% names(SV_df))) {
  stop(
    "Cannot reproduce commonCols exactly. Missing: ",
    paste(setdiff(required_common_samples, names(SV_df)), collapse = ", ")
  )
}

commonCols <- intersect(
  colnames(SV_df[["TP11-896"]]),
  colnames(SV_df[["TP21-M598"]])
)

SV_df_ <- lapply(names(SV_df), function(x) {
  SV_df[[x]] <- SV_df[[x]][, commonCols, drop = FALSE]
  SV_df[[x]]$sample <- rep(x, nrow(SV_df[[x]]))
  SV_df[[x]]
})

SV_df_merge <- do.call(rbind, SV_df_)

cat("SV rows after raw merge:", nrow(SV_df_merge), "\n")
cat("Samples after merge:", length(unique(SV_df_merge$sample)), "\n\n")

# ============================================================
# 5. OSAMA'S ACTIVE SV PREPROCESSING
# ============================================================
# Intentionally keep the original OR condition.
SV_df_merge <- subset(
  SV_df_merge,
  RefcontigID1 %notin% "24" | RefcontigID2 %notin% "24"
)

SV_df_merge$Chr <- paste0("Chr", SV_df_merge$RefcontigID1)
SV_df_merge$Chr <- gsub("23", "X", SV_df_merge$Chr)

# Partial inversions
idx <- SV_df_merge$Type %in% "inversion_partial"
SV_df_merge[idx, "RefEndPos"] <- SV_df_merge[idx, "RefStartPos"]
SV_df_merge[idx, "Type"] <- "inversion_paired"

# Same positional renaming as Osama
names(SV_df_merge)[29] <- "Pct_in_BNG_control_samples"
names(SV_df_merge)[30] <- "Pct_in_BNG_control_samples_same_enzyme"

# Same coercions
mode(SV_df_merge$RefStartPos) <- "integer"
mode(SV_df_merge$RefEndPos) <- "integer"
mode(SV_df_merge$SVsize) <- "integer"
mode(SV_df_merge$VAF) <- "numeric"
mode(SV_df_merge$SVfreq) <- "numeric"
mode(SV_df_merge$Pct_in_BNG_control_samples) <- "numeric"
mode(SV_df_merge$Pct_in_BNG_control_samples_same_enzyme) <- "numeric"
mode(SV_df_merge$Self_molecule_count) <- "numeric"
mode(SV_df_merge$Confidence) <- "numeric"

# IMPORTANT: Osama's Access 1.7 filter was commented out.
SV_df_merge_filtered <- SV_df_merge

# ============================================================
# 6. SIMPLIFY SV TYPES — SAME AS OSAMA
# ============================================================
SV_df <- SV_df_merge_filtered

SV_df$RefcontigID1 <- paste0("chr", SV_df$RefcontigID1)
SV_df$RefcontigID1 <- gsub("23", "X", SV_df$RefcontigID1)
SV_df$RefcontigID2 <- paste0("chr", SV_df$RefcontigID2)
SV_df$RefcontigID2 <- gsub("23", "X", SV_df$RefcontigID2)

SV_df[SV_df$Type %in% "deletion", "Type"] <- "Del"
SV_df[SV_df$Type %in% "duplication", "Type"] <- "Dup"
SV_df[SV_df$Type %in% "duplication_inverted", "Type"] <- "Dup"
SV_df[SV_df$Type %in% "duplication_split", "Type"] <- "Dup"
SV_df[SV_df$Type %in% "insertion", "Type"] <- "Ins"
SV_df[SV_df$Type %in% "inversion", "Type"] <- "Inv"
SV_df[SV_df$Type %in% "inversion_paired", "Type"] <- "Inv"
SV_df[SV_df$Type %in% "inversion_partial", "Type"] <- "Inv"
SV_df[SV_df$Type %in% "translocation_interchr", "Type"] <- "Transloc-inter"
SV_df[SV_df$Type %in% "translocation_intrachr", "Type"] <- "Transloc-intra"

cat("Final NON-DEDUP SV_df dimensions:\n")
print(dim(SV_df))
cat("\nSV type counts:\n")
print(table(SV_df$Type, useNA = "ifany"))

write.csv(
  SV_df,
  file = file.path(OUT_DIR, "SV_df_used_for_ShatterSeek_NONDEDUP.csv"),
  row.names = FALSE
)

write.table(
  as.data.frame(table(sample = SV_df$sample, Type = SV_df$Type)),
  file = file.path(OUT_DIR, "SV_counts_by_sample_and_type.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ============================================================
# 7. CNV MERGING FUNCTION
# ============================================================
# Same coverage-based logic as Osama.
# Only technical change: explicit GRanges() construction.
mergeHomoCNregions <- function(Cell_CN_df) {

  Cell_CN_df$total_cn[Cell_CN_df$total_cn == 0] <- 150000
  Cell_CN_df$total_cn[is.na(Cell_CN_df$total_cn)] <- 0

  dd <- GenomicRanges::GRanges(
    seqnames = as.character(Cell_CN_df$chromosome),
    ranges = IRanges::IRanges(
      start = as.numeric(Cell_CN_df$start),
      end = as.numeric(Cell_CN_df$end)
    )
  )

  S4Vectors::mcols(dd)$total_cn <- as.numeric(Cell_CN_df$total_cn)

  cov <- GenomicRanges::coverage(
    dd,
    weight = S4Vectors::mcols(dd)$total_cn
  )

  # Compatibility replacement for Osama's:
  #   dd1 <- as(cov, "GRanges")
  # Some current Bioconductor versions fail on zero-length Rle elements.
  # Extract the same coverage runs manually; biological logic is unchanged.
  dd1_list <- lapply(names(cov), function(chr) {

    x <- cov[[chr]]

    if (length(x) == 0) {
      return(NULL)
    }

    run_lengths <- as.numeric(S4Vectors::runLength(x))
    run_values  <- as.numeric(S4Vectors::runValue(x))

    if (length(run_lengths) == 0) {
      return(NULL)
    }

    run_end <- cumsum(run_lengths)
    run_start <- c(1, head(run_end, -1) + 1)

    z <- data.frame(
      chromosome = chr,
      start = run_start,
      end = run_end,
      total_cn = run_values,
      stringsAsFactors = FALSE
    )

    # Same as Osama's dd1[dd1$score != 0, ]
    z <- z[z$total_cn != 0, , drop = FALSE]

    if (nrow(z) == 0) {
      return(NULL)
    }

    z
  })

  dd1_list <- dd1_list[!vapply(dd1_list, is.null, logical(1))]

  if (length(dd1_list) == 0) {
    return(
      data.frame(
        chromosome = character(0),
        start = numeric(0),
        end = numeric(0),
        total_cn = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }

  dd1 <- do.call(rbind, dd1_list)
  rownames(dd1) <- NULL

  dd1$total_cn[dd1$total_cn == 150000] <- 0

  dd1$chromosome <- as.character(dd1$chromosome)
  dd1$start <- as.numeric(dd1$start)
  dd1$end <- as.numeric(dd1$end)
  dd1$total_cn <- as.numeric(dd1$total_cn)

  dd1
}

# ============================================================
# 8. RUN SHATTERSEEK — SAME ACTIVE OSAMA WORKFLOW
# ============================================================
chromoThrp_df_list <- list()

AllsampleIDs <- unique(SV_df$sample)

cat("\nTotal samples in SV_df:", length(AllsampleIDs), "\n")
cat("Samples entering ShatterSeek:", length(AllsampleIDs), "\n\n")

if (length(AllsampleIDs) != 62) {
  stop(
    "Expected 62 samples entering ShatterSeek but found ",
    length(AllsampleIDs)
  )
}

for (sampleID in AllsampleIDs) {

  cat("====================================================\n")
  cat("Processing:", sampleID, "\n")

  # Paired fuller Rahul CNV file for this SampleID
  cnv_file <- file_map$cnv[match(sampleID, file_map$SampleID)]

  if (length(cnv_file) != 1 || is.na(cnv_file) || !file.exists(cnv_file)) {
    stop("CNV file not resolved for sample: ", sampleID)
  }

  cat("CNV file:", basename(cnv_file), "\n")

  # ---------------- CNV ----------------
  cnv <- read.delim(cnv_file, skip = 3)
  cnv <- cnv[, c("Chromosome", "Start", "End", "CopyNumber")]

  cnv$Chromosome[cnv$Chromosome %in% 23] <- "X"
  cnv <- subset(cnv, Chromosome %notin% c("Y", "24"))

  colnames(cnv) <- c("chromosome", "start", "end", "total_cn")
  cnv$start <- as.integer(cnv$start)
  cnv$end <- as.integer(cnv$end)
  cnv$total_cn <- as.integer(cnv$total_cn)

  cnv_merged <- mergeHomoCNregions(cnv)

  if (nrow(cnv_merged) == 0) {
    stop("cnv_merged contains zero rows for ", sampleID)
  }

  cat("CNV rows before merge:", nrow(cnv), "\n")
  cat("CNV rows after merge :", nrow(cnv_merged), "\n")

  CN_data <- CNVsegs(
    chrom = as.character(cnv_merged$chromosome),
    start = as.numeric(cnv_merged$start),
    end = as.numeric(cnv_merged$end),
    total_cn = as.numeric(cnv_merged$total_cn)
  )

  # ---------------- SV ----------------
  SV <- subset(SV_df, sample %in% sampleID)

  cat("SV rows before ShatterSeek conversion:", nrow(SV), "\n")
  print(table(SV$Type, useNA = "ifany"))

  # EXACT Osama conversion block
  SV$Type[SV$Type %in% "Del"] <- "DEL"
  SV$Type[SV$Type %in% c("Transloc-inter", "Transloc-intra")] <- "TRA"
  SV$Type[SV$Type %in% "Dup"] <- "DUP"
  SV$Type[SV$Type %in% "INS"] <- "INS"

  SV$RefcontigID1 <- gsub("chr", "", SV$RefcontigID1)
  SV$RefcontigID2 <- gsub("chr", "", SV$RefcontigID2)

  # Same Osama filter for unknown second breakpoint
  SV <- subset(SV, RefcontigID2 %notin% "-1")

  SV <- SV[, c(
    "RefcontigID1",
    "RefcontigID2",
    "RefStartPos",
    "RefEndPos",
    "Type"
  )]

  cat("SV rows entering SVs():", nrow(SV), "\n")
  print(table(SV$Type, useNA = "ifany"))

  # Audit only; does not alter computation
  write.table(
    SV,
    file = file.path(OUT_DIR, paste0(sampleID, "_SVs_entering_ShatterSeek.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  SV_data <- SVs(
    chrom1 = as.character(SV$RefcontigID1),
    pos1 = as.numeric(SV$RefStartPos),
    chrom2 = as.character(SV$RefcontigID2),
    pos2 = as.numeric(SV$RefEndPos),
    SVtype = as.character(SV$Type),
    strand1 = rep("+", length(SV$Type)),
    strand2 = rep("+", length(SV$Type))
  )

  chromothripsis <- shatterseek(
    SV.sample = SV_data,
    seg.sample = CN_data,
    genome = "hg38"
  )

  chromoThrp_df <- chromothripsis@chromSummary[, c(
    "chrom",
    "start",
    "end",
    "number_DEL",
    "number_DUP",
    "number_TRA",
    "number_h2hINV",
    "number_t2tINV",
    "clusterSize_including_TRA",
    "max_number_oscillating_CN_segments_2_states",
    "chr_breakpoint_enrichment",
    "pval_exp_cluster",
    "pval_fragment_joins",
    "inter_pval_fragment_joins",
    "inter_other_chroms_coords_all"
  )]

  chromoThrp_df$total_IntraSVs <- rowSums(
    chromothripsis@chromSummary[, c(
      "number_DEL",
      "number_DUP",
      "number_h2hINV",
      "number_t2tINV"
    )]
  )

  # EXACT active Osama High / Low thresholds
  chromoThrp_df$confidence <- ifelse(
    (
      chromoThrp_df$clusterSize_including_TRA >= 12 &
      chromoThrp_df$max_number_oscillating_CN_segments_2_states >= 4 &
      chromoThrp_df$pval_fragment_joins <= 0.001 &
      (
        chromoThrp_df$chr_breakpoint_enrichment <= 0.001 |
        chromoThrp_df$pval_exp_cluster <= 0.001
      )
    ),
    "High",
    ifelse(
      (
        chromoThrp_df$total_IntraSVs >= 8 &
        chromoThrp_df$max_number_oscillating_CN_segments_2_states >= 2 &
        chromoThrp_df$pval_fragment_joins <= 0.05 &
        (
          chromoThrp_df$chr_breakpoint_enrichment <= 0.05 |
          chromoThrp_df$pval_exp_cluster <= 0.05
        )
      ),
      "Low",
      ""
    )
  )

  chromoThrp_df$sample <- sampleID
  chromoThrp_df_list[[sampleID]] <- chromoThrp_df

  cat("High calls:", sum(chromoThrp_df$confidence == "High", na.rm = TRUE), "\n")
  cat("Low calls :", sum(chromoThrp_df$confidence == "Low", na.rm = TRUE), "\n\n")
}

# ============================================================
# 9. COMBINE + SAVE SHATTERSEEK RESULTS
# ============================================================
chromoThrp_df_list_df <- do.call(rbind, chromoThrp_df_list)

write.table(
  chromoThrp_df_list_df,
  file = file.path(OUT_DIR, "ShatterSeek_df.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

save(
  chromoThrp_df_list_df,
  file = file.path(OUT_DIR, "ShatterSeek_df.Rdata")
)

# ============================================================
# 10. HEATMAP MATRIX — SAME OSAMA LOGIC
# ============================================================
thripsisMat <- reshape2::dcast(
  chromoThrp_df_list_df,
  chrom ~ sample,
  value.var = "confidence"
)

rownames(thripsisMat) <- thripsisMat$chrom
thripsisMat <- thripsisMat[, -1, drop = FALSE]
thripsisMat <- as.matrix(thripsisMat)
thripsisMat[is.na(thripsisMat)] <- ""
thripsisMat[thripsisMat == ""] <- "ns"

write.table(
  cbind(chrom = rownames(thripsisMat), thripsisMat),
  file = file.path(OUT_DIR, "Chromothripsis_heatmap_matrix_FULL.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nFull matrix:", nrow(thripsisMat), "chromosomes x", ncol(thripsisMat), "samples\n")

# Keep Osama's original hard-coded 62 / 23 filtering exactly.
thripsisMat_ <- thripsisMat[
  rowSums(thripsisMat == "ns") != 62,
  colSums(thripsisMat == "ns") != 23,
  drop = FALSE
]

write.table(
  cbind(chrom = rownames(thripsisMat_), thripsisMat_),
  file = file.path(OUT_DIR, "Chromothripsis_heatmap_matrix_Osama_filtered.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Filtered matrix:", nrow(thripsisMat_), "chromosomes x", ncol(thripsisMat_), "samples\n")

# ============================================================
# 11. OSAMA-STYLE HEATMAP
# ============================================================
# Histology annotation is omitted here because the old cohort-master-sheet
# path is not part of the current input folder. This does not change calls.
ht <- Heatmap(
  thripsisMat_,
  col = c(
    "High" = "#eded2a",
    "Low" = "#a1d99b",
    "ns" = "white"
  ),
  border = TRUE,
  name = "Event Confidence",
  row_order = rownames(thripsisMat_)[
    order(rowSums(thripsisMat_ != "ns"), decreasing = TRUE)
  ],
  column_order = colnames(thripsisMat_)[
    order(colSums(thripsisMat_ != "ns"), decreasing = TRUE)
  ],
  border_gp = gpar(col = "gray"),
  rect_gp = gpar(col = "gray", lwd = 1)
)

pdf(
  file.path(OUT_DIR, "Chromothripsis_Osama_style_heatmap.pdf"),
  width = 10,
  height = 7
)
draw(ht, merge_legends = TRUE)
dev.off()

png(
  file.path(OUT_DIR, "Chromothripsis_Osama_style_heatmap.png"),
  width = 2600,
  height = 1800,
  res = 250
)
draw(ht, merge_legends = TRUE)
dev.off()

# ============================================================
# 12. SAMPLE-LEVEL SUMMARY
# ============================================================
sample_summary <- data.frame(
  sample = AllsampleIDs,
  Total_SVs_in_SV_df = sapply(
    AllsampleIDs,
    function(x) sum(SV_df$sample == x)
  ),
  High_chromosomes = sapply(
    AllsampleIDs,
    function(x) sum(
      chromoThrp_df_list_df$sample == x &
      chromoThrp_df_list_df$confidence == "High",
      na.rm = TRUE
    )
  ),
  Low_chromosomes = sapply(
    AllsampleIDs,
    function(x) sum(
      chromoThrp_df_list_df$sample == x &
      chromoThrp_df_list_df$confidence == "Low",
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

sample_summary$Any_High_or_Low <-
  sample_summary$High_chromosomes > 0 |
  sample_summary$Low_chromosomes > 0

write.table(
  sample_summary,
  file = file.path(OUT_DIR, "Osama_style_sample_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n====================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("====================================================\n")
cat("Samples analyzed:", nrow(sample_summary), "\n")
cat("Samples with >=1 High:", sum(sample_summary$High_chromosomes > 0), "\n")
cat("Samples with >=1 Low :", sum(sample_summary$Low_chromosomes > 0), "\n")
cat("Samples with >=1 High OR Low:", sum(sample_summary$Any_High_or_Low), "\n\n")

print(sample_summary[sample_summary$Any_High_or_Low, , drop = FALSE])

cat("\nOutputs saved to:\n", OUT_DIR, "\n")

# Save package/R versions
sink(file.path(OUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")