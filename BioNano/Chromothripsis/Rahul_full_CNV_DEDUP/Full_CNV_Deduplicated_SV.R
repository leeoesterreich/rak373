#!/usr/bin/env Rscript

# ============================================================
# Osama-style Full-CNV Chromothripsis
# DEDUPLICATED SV
# ============================================================
#
# PURPOSE
# -------
# This script is a controlled comparison against the successful:
#
#   Osama_style_Full-CNV_Chromothripsis_run.R
#
# The successful reproduction used:
#
#   ORIGINAL SV INPUT:
#     Rick *_Annotated_SV.smap files
#       -> Osama active preprocessing
#       -> NON-DEDUPLICATED SV_df
#
#   CNV INPUT:
#     Rahul fuller CNV files
#
# and reproduced Osama's original chromothripsis result.
#
# ============================================================
# SCIENTIFIC QUESTION
# ============================================================
#
# Does deduplication of the structural-variant dataset alter
# ShatterSeek chromothripsis calls?
#
# To test that, this script keeps the successful Full-CNV
# Osama-style workflow and changes the SV input to:
#
#   /ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation/
#   SV_df_unqiue.csv
#
# The all-sample QC independently reconstructed Osama's
# deduplication from the original SMAP files and confirmed that
# SV_df_unqiue.csv contains the expected deduplicated SV events.
#
# ============================================================
# ONLY TWO ANALYTICAL CHANGES RELATIVE TO THE SUCCESSFUL RUN
# ============================================================
#
# CHANGE 1
# --------
# Instead of rebuilding NON-DEDUPLICATED SV_df from Rick *.smap,
# directly load:
#
#   SV_df_unqiue.csv
#
# This is the deduplicated SV representation.
#
# CHANGE 2
# --------
# SV_df_unqiue.csv stores translocations as:
#
#   Transloc-Inter
#   Transloc-Intra
#
# whereas the successful SMAP-derived table used:
#
#   Transloc-inter
#   Transloc-intra
#
# R matching is case-sensitive. Therefore the ShatterSeek conversion
# explicitly recognizes either capitalization and converts both to:
#
#   TRA
#
# using:
#
#   tolower(SV$Type) %in% c("transloc-inter", "transloc-intra")
#
# This is a label-compatibility correction; it does not add or remove
# any SV rows.
#
# ============================================================
# EVERYTHING ELSE IS HELD CONSTANT
# ============================================================
#
#   - SAME Rahul fuller CNV files
#   - SAME CNV processing
#   - SAME coverage-based homogeneous-CN merging
#   - SAME RefcontigID2 != "-1" filter immediately before SVs()
#   - NO additional strict both-breakpoint 1-22/X filter
#   - SAME ShatterSeek genome = "hg38"
#   - SAME strand1 / strand2 = "+"
#   - SAME treatment of Del / Dup / Ins / Inv
#   - SAME High-confidence thresholds
#   - SAME Low-confidence thresholds
#   - SAME Osama heatmap construction
#
# IMPORTANT:
# Do NOT add extra SV filtering to this script. The intention is to
# isolate the effect of deduplication while otherwise reproducing the
# successful Osama-style Full-CNV analysis.
#
# ============================================================
# SAMPLE-NAME COMPATIBILITY
# ============================================================
#
# One technical naming difference is known:
#
#   SV_df_unqiue.csv : TP14-079
#   fuller CNV file  : TP14-79.txt
#
# These are the same biological sample.
#
# For file matching and output comparison only, sample IDs are
# normalized so TP14-079 -> TP14-79. This is NOT an analytical
# modification and does not alter any SV event.
#
# ============================================================
# INTERPRETATION
# ============================================================
#
# Compare this run directly with the successful NON-DEDUPLICATED
# Osama-style Full-CNV run.
#
# If High/Low calls are unchanged:
#   -> deduplication does not materially affect final ShatterSeek calls.
#
# If High/Low calls change:
#   -> after the validated SV QC and corrected TRA mapping, the
#      differences can be attributed to the deduplicated SV input.
#
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

SV_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation/SV_df_unqiue.csv"

CNV_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation"

OUT_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Chromothripsis/Osama_Script/Rahul_full_CNV_DEDUP"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("SCRIPT: OSAMA-STYLE FULL-CNV + DEDUP SV + CORRECT TRA\n")
cat("SV_FILE:", SV_FILE, "\n")
cat("CNV_DIR :", CNV_DIR, "\n")
cat("OUT_DIR :", OUT_DIR, "\n\n")


# ============================================================
# 1. SAMPLE-ID NORMALIZATION FOR MATCHING ONLY
# ============================================================

normalize_sample_id <- function(x) {

  x <- trimws(
    as.character(x)
  )

  # TP14-079 -> TP14-79
  x <- sub(
    "^(TP[0-9]+)-0+([0-9]+)$",
    "\\1-\\2",
    x
  )

  x
}


# ============================================================
# 2. LOAD THE VALIDATED DEDUPLICATED SV TABLE
# ============================================================
#
# This block replaces the SMAP loading + merging + preprocessing
# sections of the successful Osama-style script.
#
# SV_df_unqiue.csv is already the deduplicated, simplified SV table.
# No additional deduplication or analytical filtering is performed.
# ============================================================

SV_df <- read.csv(
  SV_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Preserve the original value for provenance.
SV_df$sample_original <- as.character(
  SV_df$sample
)

# Normalize only the sample ID used by this controlled comparison.
SV_df$sample <- normalize_sample_id(
  SV_df$sample
)

SV_df$RefcontigID1 <- as.character(
  SV_df$RefcontigID1
)

SV_df$RefcontigID2 <- as.character(
  SV_df$RefcontigID2
)

SV_df$RefStartPos <- as.integer(
  SV_df$RefStartPos
)

SV_df$RefEndPos <- as.integer(
  SV_df$RefEndPos
)

SV_df$Type <- as.character(
  SV_df$Type
)

cat("DEDUP SV_df dimensions:\n")
print(dim(SV_df))

cat("\nDEDUP SV type counts before ShatterSeek conversion:\n")
print(
  table(
    SV_df$Type,
    useNA = "ifany"
  )
)

cat(
  "\nNormalized biological samples in SV_df:",
  length(unique(SV_df$sample)),
  "\n\n"
)

if (length(unique(SV_df$sample)) != 62) {
  stop(
    "Expected 62 normalized biological samples in SV_df_unqiue.csv but found ",
    length(unique(SV_df$sample))
  )
}

write.csv(
  SV_df,
  file = file.path(
    OUT_DIR,
    "SV_df_used_for_ShatterSeek_DEDUP.csv"
  ),
  row.names = FALSE
)

write.table(
  as.data.frame(
    table(
      sample = SV_df$sample,
      Type = SV_df$Type
    )
  ),
  file = file.path(
    OUT_DIR,
    "SV_counts_by_sample_and_type_DEDUP.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 3. PAIR THE SAME FULLER RAHUL CNV FILES
# ============================================================
#
# CNV source is identical to the successful Full-CNV run.
# Sample-name normalization only resolves TP14-079 vs TP14-79.
# ============================================================

cnvFilesPaths <- list.files(
  CNV_DIR,
  pattern = "\\.txt$",
  full.names = TRUE
)

cnv_ids_original <- sub(
  "\\.txt$",
  "",
  basename(cnvFilesPaths)
)

cnv_ids_norm <- normalize_sample_id(
  cnv_ids_original
)

AllsampleIDs <- unique(
  SV_df$sample
)

file_map <- data.frame(
  SampleID = AllsampleIDs,
  stringsAsFactors = FALSE
)

file_map$cnv <- cnvFilesPaths[
  match(
    file_map$SampleID,
    cnv_ids_norm
  )
]

if (any(is.na(file_map$cnv))) {

  print(
    file_map[
      is.na(file_map$cnv),
      ,
      drop = FALSE
    ]
  )

  stop(
    "At least one deduplicated SV sample could not be paired with the same fuller Rahul CNV file."
  )
}

write.table(
  file_map,
  file = file.path(
    OUT_DIR,
    "DEDUP_SV_Full_CNV_file_map.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Samples paired to fuller CNV files:",
  nrow(file_map),
  "\n\n"
)

# ============================================================
# 4. CNV MERGING FUNCTION
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
# 5. RUN SHATTERSEEK — SAME ACTIVE OSAMA WORKFLOW
# ============================================================
chromoThrp_df_list <- list()

AllsampleIDs <- unique(SV_df$sample)

cat("\nTotal samples in DEDUP SV_df:", length(AllsampleIDs), "\n")
cat("Samples entering ShatterSeek:", length(AllsampleIDs), "\n\n")

if (length(AllsampleIDs) != 62) {
  stop(
    "Expected 62 normalized biological samples entering ShatterSeek but found ",
    length(AllsampleIDs)
  )
}

for (sampleID in AllsampleIDs) {

  cat("====================================================\n")
  cat("Processing:", sampleID, "\n")

  # SAME fuller Rahul CNV file used in the successful Osama-style run
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

  cat("DEDUP SV rows before ShatterSeek conversion:", nrow(SV), "\n")
  print(table(SV$Type, useNA = "ifany"))

  # EXACT Osama conversion logic, with one necessary label-compatibility fix:
  # SV_df_unqiue.csv contains Transloc-Inter / Transloc-Intra,
  # so normalize case ONLY for recognizing translocations.
  SV$Type[SV$Type %in% "Del"] <- "DEL"
  SV$Type[tolower(SV$Type) %in% c("transloc-inter", "transloc-intra")] <- "TRA"
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
    file = file.path(OUT_DIR, paste0(sampleID, "_DEDUP_SVs_entering_ShatterSeek.tsv")),
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
# 6. COMBINE + SAVE SHATTERSEEK RESULTS
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
# 7. HEATMAP MATRIX — SAME OSAMA LOGIC
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
# 8. OSAMA-STYLE HEATMAP
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
# 9. SAMPLE-LEVEL SUMMARY
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
  file = file.path(OUT_DIR, "Osama_style_DEDUP_correct_TRA_sample_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n====================================================\n")
cat("DEDUPLICATED-SV + CORRECT TRA ANALYSIS COMPLETE\n")
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