suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
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

required <- c("motifdir", "outdir")
missing <- setdiff(required, names(opt))

if (length(missing) > 0) {
  stop("Missing required args: ", paste(missing, collapse = ", "))
}

motifdir <- opt$motifdir
outdir   <- opt$outdir

message("motifdir: ", motifdir)
message("outdir:   ", outdir)

if (!dir.exists(motifdir)) {
  stop("motifdir does not exist: ", motifdir)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read all per-sample GC RDS files
# -----------------------------
files <- list.files(
  motifdir,
  pattern = "_4bpmotif_gc\\.rds$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *_4bpmotif_gc.rds files found in: ", motifdir)
}

ids <- basename(files)
ids <- sub("_4bpmotif_gc\\.rds$", "", ids)

message("Combining ", length(files), " samples...")

files.list <- lapply(files, readRDS)

# Combine into one tibble with sample id
tib_all <- map2(files.list, ids, ~ as_tibble(.x) %>% mutate(id = .y)) %>%
  bind_rows() %>%
  select(id, everything())

if (!("end_motif" %in% colnames(tib_all))) {
  stop(
    "Missing 'end_motif' column. Available columns: ",
    paste(colnames(tib_all), collapse = ", ")
  )
}

message("All motif table dimensions:")
print(dim(tib_all))

message("Number of samples combined: ", length(unique(tib_all$id)))
message("Unique motifs before ATGC filtering: ", n_distinct(tib_all$end_motif))

# -----------------------------
# Write all motif outputs
# -----------------------------
out_all_rds <- file.path(outdir, "endmotif_4bp_gc_summary.rds")
out_all_tsv <- file.path(outdir, "endmotif_4bp_gc_summary.tsv.gz")

saveRDS(tib_all, out_all_rds)
write_tsv(tib_all, out_all_tsv)

message("Wrote all motif RDS: ", out_all_rds)
message("Wrote all motif TSV: ", out_all_tsv)

# -----------------------------
# ATGC-only + 4bp filtering
# -----------------------------
tib_atgc <- tib_all %>%
  mutate(end_motif = toupper(str_trim(end_motif))) %>%
  filter(!is.na(end_motif)) %>%
  filter(str_length(end_motif) == 4) %>%
  filter(!grepl("[^ATGC]", end_motif))

message("ATGC-only table dimensions:")
print(dim(tib_atgc))

message("Unique motifs after ATGC filtering: ", n_distinct(tib_atgc$end_motif))
message("ATGC 4-mers present globally: ", n_distinct(tib_atgc$end_motif), " expected up to 256")
message("Number of samples in ATGC-only table: ", length(unique(tib_atgc$id)))

# -----------------------------
# Write ATGC-only outputs
# -----------------------------
out_atgc_rds <- file.path(outdir, "endmotif_4bp_gc_summary_ATGC_only.rds")
out_atgc_tsv <- file.path(outdir, "endmotif_4bp_gc_summary_ATGC_only.tsv.gz")

saveRDS(tib_atgc, out_atgc_rds)
write_tsv(tib_atgc, out_atgc_tsv)

message("Wrote ATGC-only RDS: ", out_atgc_rds)
message("Wrote ATGC-only TSV: ", out_atgc_tsv)

# -----------------------------
# Final output check
# -----------------------------
expected_files <- c(
  out_all_rds,
  out_all_tsv,
  out_atgc_rds,
  out_atgc_tsv
)

message("Checking output files:")

for (ff in expected_files) {
  if (file.exists(ff) && file.info(ff)$size > 0) {
    message("FOUND: ", ff, " | size: ", round(file.info(ff)$size / 1024^2, 3), " MB")
  } else {
    warning("MISSING or EMPTY: ", ff)
  }
}

message("DONE: combined motif outputs written to: ", outdir)