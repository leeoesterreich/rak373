```r
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# -----------------------------
# Inputs / outputs
# -----------------------------
in_rds  <- "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif/combined_motif/endmotif_4bp_gc_summary.rds"
outdir  <- "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif/motif_analysis_no_annot"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(outdir, "endmotif_4bp_freq_per_sample.rds")
out_tsv <- file.path(outdir, "endmotif_4bp_freq_per_sample.tsv.gz")

# -----------------------------
# Load combined tibble
# -----------------------------
df <- readRDS(in_rds)

# Expected minimum columns: id, end_motif, count, count.corrected
required <- c("id", "end_motif", "count", "count.corrected")
missing <- setdiff(required, colnames(df))
if (length(missing) > 0) {
  stop("Missing required columns in input RDS: ", paste(missing, collapse = ", "))
}

# If you have 'frac' column keep it; otherwise assume single fraction
has_frac <- "frac" %in% colnames(df)

# -----------------------------
# 1) Collapse across GC bins (and other per-bin fields)
# -----------------------------
if (has_frac) {
  df_sum <- df %>%
    group_by(id, frac, end_motif) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      count.corrected = sum(count.corrected, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  df_sum <- df %>%
    group_by(id, end_motif) %>%
    summarise(
      count = sum(count, na.rm = TRUE),
      count.corrected = sum(count.corrected, na.rm = TRUE),
      .groups = "drop"
    )
}

# -----------------------------
# 2) Compute per-sample frequencies
# -----------------------------
if (has_frac) {
  df_freq <- df_sum %>%
    group_by(id, frac) %>%
    mutate(
      total = sum(count, na.rm = TRUE),
      total.corrected = sum(count.corrected, na.rm = TRUE),
      freq = ifelse(total > 0, count / total, NA_real_),
      freq.corrected = ifelse(total.corrected > 0, count.corrected / total.corrected, NA_real_)
    ) %>%
    ungroup()
} else {
  df_freq <- df_sum %>%
    group_by(id) %>%
    mutate(
      total = sum(count, na.rm = TRUE),
      total.corrected = sum(count.corrected, na.rm = TRUE),
      freq = ifelse(total > 0, count / total, NA_real_),
      freq.corrected = ifelse(total.corrected > 0, count.corrected / total.corrected, NA_real_)
    ) %>%
    ungroup()
}

# -----------------------------
# Save outputs
# -----------------------------
saveRDS(df_freq, out_rds)
write_tsv(df_freq, out_tsv)

message("Wrote RDS: ", out_rds)
message("Wrote TSV: ", out_tsv)

# -----------------------------
# Optional quick sanity prints (motif-level)
# -----------------------------
message("Samples: ", dplyr::n_distinct(df_freq$id))

motifs_per_sample <- df_freq %>%
  group_by(id) %>%
  summarise(n_motifs = n_distinct(end_motif), .groups = "drop")

message("Motifs per sample (median): ",
        median(motifs_per_sample$n_motifs))

# Debug info (optional)
message("Has 'frac' column: ", "frac" %in% colnames(df_freq))
if ("frac" %in% colnames(df_freq)) {
  message("Number of frac levels: ", length(unique(df_freq$frac)))
}
```
