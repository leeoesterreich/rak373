#!/usr/bin/env Rscript

# ============================================================
# COMPLETE QC: Validate SV_df_unqiue.csv against original SMAP files
# ============================================================
#
# PURPOSE
# -------
# Before interpreting the difference between the Osama-style
# NON-DEDUPLICATED ShatterSeek result and the DEDUPLICATED result,
# verify that:
#
#   /ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation/
#   SV_df_unqiue.csv
#
# is truly the expected deduplicated version of the original
# Bionano SMAP-derived SV table.
#
# This QC is performed for the ENTIRE cohort, not just a few
# randomly selected samples.
#
# ============================================================
# WHAT THIS SCRIPT CHECKS
# ============================================================
#
# For every original Rick *_Annotated_SV.smap sample:
#
#   1. Read the original SMAP file.
#
#   2. Apply the same ACTIVE Osama preprocessing that was used
#      before deduplication:
#
#        - read SMAP with skip = 6
#        - remove the extra first row
#        - restrict to the same common SMAP columns
#        - add sample ID
#        - apply Osama's chromosome-24 OR filter
#        - create Chr from RefcontigID1
#        - convert chromosome 23 -> X
#        - convert inversion_partial to inversion_paired and set
#          RefEndPos = RefStartPos
#        - keep Access 1.7 filtering INACTIVE
#        - simplify SV types to:
#             Del
#             Dup
#             Ins
#             Inv
#             Transloc-inter
#             Transloc-intra
#        - add "chr" prefixes to RefcontigID1 / RefcontigID2
#
#   3. Recreate the Osama deduplication key:
#
#        paste0(
#          ISCN,
#          sample,
#          Chr,
#          PutativeGeneFusion,
#          OverlapGenes
#        )
#
#   4. Retain the first row for each unique key.
#
#   5. Compare the independently recomputed deduplicated rows
#      against the corresponding rows in SV_df_unqiue.csv.
#
# ============================================================
# SAMPLE-NAME NORMALIZATION
# ============================================================
#
# One known naming difference exists:
#
#      original SMAP ID      : TP14-79
#      SV_df_unqiue.csv ID   : TP14-079
#
# These are the same biological sample.
#
# For matching and dedup-key comparison only, the script normalizes:
#
#      TP14-079 -> TP14-79
#
# Original sample names are retained in the output tables so that
# this normalization remains transparent.
#
# ============================================================
# EXPECTED PASS CRITERIA
# ============================================================
#
# For every sample:
#
#   Recomputed_dedup_rows == Rows_in_SV_df_unqiue
#
#   All_recomputed_keys_in_saved == TRUE
#
#   All_saved_keys_in_recomputed == TRUE
#
#   Key_sets_identical == TRUE
#
# The script also checks the cohort-level totals:
#
#   total recomputed deduplicated rows
#          ==
#   total rows in SV_df_unqiue.csv
#
# If all samples pass, SV_df_unqiue.csv can be considered a valid
# deduplicated counterpart of the original SMAP-derived SV table
# for the purpose of the ShatterSeek comparison.
#
# IMPORTANT
# ---------
# This is a QC script only.
# It does NOT run ShatterSeek and does NOT use CNV files.
#
# ============================================================


# ============================================================
# 1. PATHS
# ============================================================

SV_DIR <- "/ix1/alee/LO_LAB/Personal/Rick/Bionano/20231216_Liane_Osama_Bionano_DNA_SV"

UNIQUE_FILE <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/CNV_Annotation/SV_df_unqiue.csv"

OUT_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Chromothripsis/SV_dedup_QC_All62"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

`%notin%` <- Negate(`%in%`)


# ============================================================
# 2. SAMPLE-ID NORMALIZATION
# ============================================================

normalize_sample_id <- function(x) {

  x <- trimws(
    as.character(x)
  )

  # Example:
  # TP14-079 -> TP14-79
  x <- sub(
    "^(TP[0-9]+)-0+([0-9]+)$",
    "\\1-\\2",
    x
  )

  x
}


# ============================================================
# 3. IDENTIFY ORIGINAL SMAP FILES
# ============================================================

smap_files <- list.files(
  SV_DIR,
  pattern = "\\.smap$",
  full.names = TRUE
)

smap_ids_original <- gsub(
  x = basename(smap_files),
  pattern = "Annotated_SV.smap|_12_16_|_12_17_|2023_[0-9]+_[0-9]+_[0-9]+_|_|_cr|_s|_b|_A",
  replacement = ""
)

smap_ids_norm <- normalize_sample_id(
  smap_ids_original
)

smap_map <- setNames(
  smap_files,
  smap_ids_norm
)

cat(
  "Original SMAP files:",
  length(smap_files),
  "\n"
)

cat(
  "Unique normalized SMAP samples:",
  length(unique(smap_ids_norm)),
  "\n\n"
)


# ============================================================
# 4. LOAD SAVED DEDUPLICATED SV TABLE
# ============================================================

saved_unique <- read.csv(
  UNIQUE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

saved_unique$sample_original <- saved_unique$sample

saved_unique$sample_norm <- normalize_sample_id(
  saved_unique$sample
)

cat(
  "Rows in SV_df_unqiue.csv:",
  nrow(saved_unique),
  "\n"
)

cat(
  "Unique normalized samples in SV_df_unqiue.csv:",
  length(unique(saved_unique$sample_norm)),
  "\n\n"
)


# ============================================================
# 5. COHORT MATCHING
# ============================================================

smap_samples <- sort(
  unique(smap_ids_norm)
)

saved_samples <- sort(
  unique(saved_unique$sample_norm)
)

missing_from_saved <- setdiff(
  smap_samples,
  saved_samples
)

extra_in_saved <- setdiff(
  saved_samples,
  smap_samples
)

cat(
  "Samples present in SMAP but missing from SV_df_unqiue:",
  length(missing_from_saved),
  "\n"
)

if (length(missing_from_saved) > 0) {
  print(missing_from_saved)
}

cat(
  "Samples present in SV_df_unqiue but not SMAP:",
  length(extra_in_saved),
  "\n"
)

if (length(extra_in_saved) > 0) {
  print(extra_in_saved)
}

cat("\n")

if (
  length(missing_from_saved) > 0 ||
  length(extra_in_saved) > 0
) {
  stop(
    "Sample sets do not match after normalization."
  )
}

test_samples <- smap_samples


# ============================================================
# 6. READ RAW SMAP
# ============================================================

read_smap_raw <- function(file) {

  x <- read.delim(
    file,
    skip = 6,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Same as Osama
  x <- x[-1, , drop = FALSE]

  x
}


# ============================================================
# 7. DEFINE THE SAME COMMON SMAP COLUMNS AS OSAMA
# ============================================================
#
# Osama used TP11-896 and TP21-M598 to establish commonCols.
# ============================================================

ref1 <- read_smap_raw(
  smap_map[["TP11-896"]]
)

ref2 <- read_smap_raw(
  smap_map[["TP21-M598"]]
)

commonCols <- intersect(
  colnames(ref1),
  colnames(ref2)
)

cat(
  "Common SMAP columns:",
  length(commonCols),
  "\n\n"
)


# ============================================================
# 8. RECREATE ACTIVE OSAMA SV PREPROCESSING
# ============================================================

prepare_osama_sv <- function(sampleID_norm) {

  file <- smap_map[[sampleID_norm]]

  x <- read_smap_raw(
    file
  )

  raw_n <- nrow(x)

  # Same common-column restriction.
  x <- x[
    ,
    commonCols,
    drop = FALSE
  ]

  # Use the SMAP-derived sample ID first.
  smap_original_id <- smap_ids_original[
    match(
      sampleID_norm,
      smap_ids_norm
    )
  ]

  x$sample <- smap_original_id
  x$sample_norm <- normalize_sample_id(
    x$sample
  )

  # ----------------------------------------------------------
  # Same active Osama chromosome-24 filter:
  #
  #   RefcontigID1 != 24 OR RefcontigID2 != 24
  #
  # Thus only rows where BOTH are chromosome 24 are removed.
  # ----------------------------------------------------------

  x <- subset(
    x,
    RefcontigID1 %notin% "24" |
    RefcontigID2 %notin% "24"
  )

  # Same Chr field.
  x$Chr <- paste0(
    "Chr",
    x$RefcontigID1
  )

  x$Chr <- gsub(
    "23",
    "X",
    x$Chr
  )

  # Same partial-inversion handling.
  idx_partial <- x$Type %in% "inversion_partial"

  x[
    idx_partial,
    "RefEndPos"
  ] <- x[
    idx_partial,
    "RefStartPos"
  ]

  x[
    idx_partial,
    "Type"
  ] <- "inversion_paired"

  # Same positional renaming present in Osama's workflow.
  if (ncol(x) >= 30) {

    names(x)[29] <-
      "Pct_in_BNG_control_samples"

    names(x)[30] <-
      "Pct_in_BNG_control_samples_same_enzyme"
  }

  # Same coordinate coercion.
  x$RefStartPos <- as.integer(
    x$RefStartPos
  )

  x$RefEndPos <- as.integer(
    x$RefEndPos
  )

  # ----------------------------------------------------------
  # IMPORTANT:
  # Access 1.7 filtering was commented out in Osama's script.
  # Therefore no additional SV filtering is applied here.
  # ----------------------------------------------------------

  # Same Refcontig representation.
  x$RefcontigID1 <- paste0(
    "chr",
    x$RefcontigID1
  )

  x$RefcontigID1 <- gsub(
    "23",
    "X",
    x$RefcontigID1
  )

  x$RefcontigID2 <- paste0(
    "chr",
    x$RefcontigID2
  )

  x$RefcontigID2 <- gsub(
    "23",
    "X",
    x$RefcontigID2
  )

  # ----------------------------------------------------------
  # Same simplified SV types
  # ----------------------------------------------------------

  x[
    x$Type == "deletion",
    "Type"
  ] <- "Del"

  x[
    x$Type %in% c(
      "duplication",
      "duplication_inverted",
      "duplication_split"
    ),
    "Type"
  ] <- "Dup"

  x[
    x$Type == "insertion",
    "Type"
  ] <- "Ins"

  x[
    x$Type %in% c(
      "inversion",
      "inversion_paired",
      "inversion_partial"
    ),
    "Type"
  ] <- "Inv"

  x[
    x$Type == "translocation_interchr",
    "Type"
  ] <- "Transloc-inter"

  x[
    x$Type == "translocation_intrachr",
    "Type"
  ] <- "Transloc-intra"

  attr(
    x,
    "raw_n"
  ) <- raw_n

  x
}


# ============================================================
# 9. DEDUPLICATION KEY
# ============================================================
#
# The original Osama key used:
#
#   ISCN + sample + Chr + PutativeGeneFusion + OverlapGenes
#
# For this QC only, sample_norm is used instead of the literal
# sample string so TP14-79 and TP14-079 are treated as the same
# biological sample.
# ============================================================

make_dedup_key_norm <- function(x) {

  paste0(
    x$ISCN,
    x$sample_norm,
    x$Chr,
    x$PutativeGeneFusion,
    x$OverlapGenes
  )
}


# ============================================================
# 10. RUN QC FOR ALL SAMPLES
# ============================================================

summary_list <- list()

all_recomputed_unique <- list()

for (sampleID in test_samples) {

  cat(
    "====================================================\n"
  )

  cat(
    "Checking:",
    sampleID,
    "\n"
  )

  processed <- prepare_osama_sv(
    sampleID
  )

  raw_n <- attr(
    processed,
    "raw_n"
  )

  processed$key_dedup_norm <-
    make_dedup_key_norm(
      processed
    )

  recomputed_unique <- processed[
    !duplicated(
      processed$key_dedup_norm
    ),
    ,
    drop = FALSE
  ]

  # Saved rows for the same normalized biological sample.
  saved_sample <- saved_unique[
    saved_unique$sample_norm == sampleID,
    ,
    drop = FALSE
  ]

  saved_sample$key_dedup_norm <-
    make_dedup_key_norm(
      saved_sample
    )

  recomputed_keys <- recomputed_unique$key_dedup_norm

  saved_keys <- saved_sample$key_dedup_norm

  all_recomputed_in_saved <- all(
    recomputed_keys %in% saved_keys
  )

  all_saved_in_recomputed <- all(
    saved_keys %in% recomputed_keys
  )

  key_sets_identical <- setequal(
    recomputed_keys,
    saved_keys
  )

  duplicate_rows_removed <-
    nrow(processed) -
    nrow(recomputed_unique)

  summary_list[[sampleID]] <- data.frame(

    SampleID_normalized = sampleID,

    SMAP_SampleID = unique(
      processed$sample
    )[1],

    Saved_SampleID = if (
      nrow(saved_sample) > 0
    ) {
      unique(
        saved_sample$sample_original
      )[1]
    } else {
      NA_character_
    },

    Raw_SMAP_rows = raw_n,

    Osama_processed_rows_before_dedup =
      nrow(processed),

    Duplicate_rows_removed =
      duplicate_rows_removed,

    Recomputed_dedup_rows =
      nrow(recomputed_unique),

    Rows_in_SV_df_unqiue =
      nrow(saved_sample),

    Row_counts_identical =
      nrow(recomputed_unique) ==
      nrow(saved_sample),

    All_recomputed_keys_in_saved =
      all_recomputed_in_saved,

    All_saved_keys_in_recomputed =
      all_saved_in_recomputed,

    Key_sets_identical =
      key_sets_identical,

    stringsAsFactors = FALSE
  )

  all_recomputed_unique[[sampleID]] <-
    recomputed_unique

  cat(
    "Raw SMAP rows                :",
    raw_n,
    "\n"
  )

  cat(
    "Processed rows before dedup  :",
    nrow(processed),
    "\n"
  )

  cat(
    "Duplicate rows removed       :",
    duplicate_rows_removed,
    "\n"
  )

  cat(
    "Recomputed dedup rows        :",
    nrow(recomputed_unique),
    "\n"
  )

  cat(
    "Rows in SV_df_unqiue.csv     :",
    nrow(saved_sample),
    "\n"
  )

  cat(
    "Row counts identical         :",
    nrow(recomputed_unique) ==
      nrow(saved_sample),
    "\n"
  )

  cat(
    "All recomputed keys in saved :",
    all_recomputed_in_saved,
    "\n"
  )

  cat(
    "All saved keys in recomputed :",
    all_saved_in_recomputed,
    "\n"
  )

  cat(
    "Key sets identical           :",
    key_sets_identical,
    "\n\n"
  )


  # ==========================================================
  # SAVE MISMATCHES ONLY IF PRESENT
  # ==========================================================

  recomputed_not_saved <- recomputed_unique[
    !(
      recomputed_unique$key_dedup_norm
      %in%
      saved_keys
    ),
    ,
    drop = FALSE
  ]

  saved_not_recomputed <- saved_sample[
    !(
      saved_sample$key_dedup_norm
      %in%
      recomputed_keys
    ),
    ,
    drop = FALSE
  ]

  if (nrow(recomputed_not_saved) > 0) {

    write.csv(
      recomputed_not_saved,
      file.path(
        OUT_DIR,
        paste0(
          sampleID,
          "_RECOMPUTED_NOT_IN_SAVED.csv"
        )
      ),
      row.names = FALSE
    )
  }

  if (nrow(saved_not_recomputed) > 0) {

    write.csv(
      saved_not_recomputed,
      file.path(
        OUT_DIR,
        paste0(
          sampleID,
          "_SAVED_NOT_IN_RECOMPUTED.csv"
        )
      ),
      row.names = FALSE
    )
  }
}


# ============================================================
# 11. SAMPLE-LEVEL QC SUMMARY
# ============================================================

qc_summary <- do.call(
  rbind,
  summary_list
)

rownames(qc_summary) <- NULL

write.csv(
  qc_summary,
  file.path(
    OUT_DIR,
    "SV_dedup_All62_QC_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. COHORT-LEVEL COUNTS
# ============================================================

recomputed_all <- do.call(
  rbind,
  all_recomputed_unique
)

rownames(recomputed_all) <- NULL

total_raw_smap_rows <- sum(
  qc_summary$Raw_SMAP_rows
)

total_processed_rows <- sum(
  qc_summary$Osama_processed_rows_before_dedup
)

total_duplicates_removed <- sum(
  qc_summary$Duplicate_rows_removed
)

total_recomputed_dedup <- sum(
  qc_summary$Recomputed_dedup_rows
)

total_saved_unique <- nrow(
  saved_unique
)

cohort_summary <- data.frame(

  Number_of_samples =
    nrow(qc_summary),

  Total_raw_SMAP_rows =
    total_raw_smap_rows,

  Total_Osama_processed_rows_before_dedup =
    total_processed_rows,

  Total_duplicate_rows_removed =
    total_duplicates_removed,

  Total_recomputed_dedup_rows =
    total_recomputed_dedup,

  Total_rows_in_SV_df_unqiue =
    total_saved_unique,

  Total_row_counts_identical =
    total_recomputed_dedup ==
    total_saved_unique,

  All_samples_row_counts_identical =
    all(
      qc_summary$Row_counts_identical
    ),

  All_samples_key_sets_identical =
    all(
      qc_summary$Key_sets_identical
    ),

  stringsAsFactors = FALSE
)

write.csv(
  cohort_summary,
  file.path(
    OUT_DIR,
    "SV_dedup_All62_COHORT_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. TYPE-LEVEL COUNT COMPARISON
# ============================================================
#
# This is an additional cohort-level sanity check.
# ============================================================

type_recomputed <- as.data.frame(
  table(
    recomputed_all$Type
  ),
  stringsAsFactors = FALSE
)

names(type_recomputed) <- c(
  "Type",
  "Recomputed_dedup_count"
)

type_saved <- as.data.frame(
  table(
    saved_unique$Type
  ),
  stringsAsFactors = FALSE
)

names(type_saved) <- c(
  "Type",
  "Saved_SV_df_unqiue_count"
)

type_comparison <- merge(
  type_recomputed,
  type_saved,
  by = "Type",
  all = TRUE
)

type_comparison[
  is.na(type_comparison)
] <- 0

type_comparison$Counts_identical <-
  type_comparison$Recomputed_dedup_count ==
  type_comparison$Saved_SV_df_unqiue_count

write.csv(
  type_comparison,
  file.path(
    OUT_DIR,
    "SV_dedup_All62_Type_count_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. FINAL REPORT
# ============================================================

cat(
  "\n====================================================\n"
)

cat(
  "FINAL ALL-SAMPLE SV DEDUP QC\n"
)

cat(
  "====================================================\n\n"
)

print(
  cohort_summary,
  row.names = FALSE
)

cat(
  "\nSV type comparison:\n"
)

print(
  type_comparison,
  row.names = FALSE
)

cat(
  "\nSamples failing row-count comparison:\n"
)

failed_row_counts <- qc_summary[
  !qc_summary$Row_counts_identical,
  ,
  drop = FALSE
]

if (nrow(failed_row_counts) == 0) {

  cat(
    "None\n"
  )

} else {

  print(
    failed_row_counts,
    row.names = FALSE
  )
}

cat(
  "\nSamples failing dedup-key comparison:\n"
)

failed_keys <- qc_summary[
  !qc_summary$Key_sets_identical,
  ,
  drop = FALSE
]

if (nrow(failed_keys) == 0) {

  cat(
    "None\n"
  )

} else {

  print(
    failed_keys,
    row.names = FALSE
  )
}


cat(
  "\n====================================================\n"
)

if (
  nrow(qc_summary) == 62 &&
  cohort_summary$Total_row_counts_identical &&
  cohort_summary$All_samples_row_counts_identical &&
  cohort_summary$All_samples_key_sets_identical &&
  all(type_comparison$Counts_identical)
) {

  cat(
    "PASS\n"
  )

  cat(
    "SV_df_unqiue.csv exactly matches the independently\n"
  )

  cat(
    "recomputed Osama-style deduplicated SV table across\n"
  )

  cat(
    "all 62 normalized biological samples.\n"
  )

} else {

  cat(
    "QC NOT FULLY PASSED\n"
  )

  cat(
    "Inspect the sample-level summary and mismatch files.\n"
  )
}

cat(
  "====================================================\n"
)

cat(
  "\nOutput directory:\n",
  OUT_DIR,
  "\n"
)