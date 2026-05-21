suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# -----------------------------
# Arg parsing
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  
  while (i <= length(args)) {
    key <- args[[i]]
    
    if (!startsWith(key, "--")) {
      stop("Invalid arg: ", key)
    }
    
    key <- sub("^--", "", key)
    
    # allow --key=value
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- kv[[2]]
      i <- i + 1
      next
    }
    
    if (i == length(args)) {
      stop("Missing value for --", key)
    }
    
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  
  out
}

opt <- parse_args(args)

required <- c("in_rds", "outdir", "prefix")
missing <- setdiff(required, names(opt))

if (length(missing) > 0) {
  stop("Missing required args: ", paste(missing, collapse = ", "))
}

in_rds <- opt$in_rds
outdir <- opt$outdir
prefix <- opt$prefix

if (!file.exists(in_rds)) {
  stop("Input RDS does not exist: ", in_rds)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(outdir, paste0(prefix, ".rds"))
out_tsv <- file.path(outdir, paste0(prefix, ".tsv.gz"))

message("Input RDS: ", in_rds)
message("Output dir: ", outdir)
message("Output prefix: ", prefix)

# -----------------------------
# Load combined tibble
# -----------------------------
df <- readRDS(in_rds)

# Expected minimum columns
required_cols <- c("id", "end_motif", "count", "count.corrected")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in input RDS: ",
    paste(missing_cols, collapse = ", ")
  )
}

has_frac <- "frac" %in% colnames(df)

message("Rows in input: ", nrow(df))
message("Samples in input: ", dplyr::n_distinct(df$id))
message("Unique motifs in input: ", dplyr::n_distinct(df$end_motif))
message("Has 'frac' column: ", has_frac)

# -----------------------------
# 1) Collapse across GC bins and other per-bin fields
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
      freq.corrected = ifelse(
        total.corrected > 0,
        count.corrected / total.corrected,
        NA_real_
      )
    ) %>%
    ungroup()
  
} else {
  
  df_freq <- df_sum %>%
    group_by(id) %>%
    mutate(
      total = sum(count, na.rm = TRUE),
      total.corrected = sum(count.corrected, na.rm = TRUE),
      freq = ifelse(total > 0, count / total, NA_real_),
      freq.corrected = ifelse(
        total.corrected > 0,
        count.corrected / total.corrected,
        NA_real_
      )
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
# Sanity checks
# -----------------------------
message("Samples in output: ", dplyr::n_distinct(df_freq$id))
message("Unique motifs in output: ", dplyr::n_distinct(df_freq$end_motif))

motifs_per_sample <- df_freq %>%
  group_by(id) %>%
  summarise(n_motifs = n_distinct(end_motif), .groups = "drop")

message("Motifs per sample median: ", median(motifs_per_sample$n_motifs))

if ("frac" %in% colnames(df_freq)) {
  message("Number of frac levels: ", length(unique(df_freq$frac)))
}

message("DONE")