```r
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
})

# -----------------------------
# Arg parsing (robust)
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Invalid arg: ", key)
    key <- sub("^--", "", key)

    # allow --key=value
    if (grepl("=", key)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      out[[kv[[1]]]] <- kv[[2]]
      i <- i + 1
      next
    }

    if (i == length(args)) stop("Missing value for --", key)
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opt <- parse_args(args)

required <- c("motifdir", "outdir")
missing <- setdiff(required, names(opt))
if (length(missing) > 0)
  stop("Missing required args: ", paste(missing, collapse = ", "))

motifdir <- opt$motifdir
outdir   <- opt$outdir

if (!dir.exists(motifdir))
  stop("motifdir does not exist: ", motifdir)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read all per-sample GC RDS files
# -----------------------------
files <- list.files(motifdir, pattern = "_4bpmotif_gc\\.rds$", full.names = TRUE)
if (length(files) == 0)
  stop("No *_4bpmotif_gc.rds files found in: ", motifdir)

ids <- basename(files)
ids <- sub("_4bpmotif_gc\\.rds$", "", ids)

message("Combining ", length(files), " samples...")

files.list <- lapply(files, readRDS)

# Combine into one tibble with sample id
tib.list <- map2(files.list, ids, ~ as_tibble(.x) %>% mutate(id = .y)) %>%
  bind_rows() %>%
  select(id, everything())

# -----------------------------
# Write outputs
# -----------------------------
out.file <- file.path(outdir, "endmotif_4bp_gc_summary.rds")
saveRDS(tib.list, out.file)

tsv.file <- file.path(outdir, "endmotif_4bp_gc_summary.tsv.gz")
write_tsv(tib.list, tsv.file)

message("Wrote RDS: ", out.file)
message("Wrote TSV: ", tsv.file)

```
