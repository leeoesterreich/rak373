```r
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
})

# ---- set paths here ----
motifdir <- "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif/motif_gc"
outdir   <- "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/Broad/NST/End_motif/combined_motif"

stopifnot(dir.exists(motifdir))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read all per-sample GC RDS files
# -----------------------------
files <- list.files(motifdir, pattern = "_4bpmotif_gc\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No *_4bpmotif_gc.rds files found in: ", motifdir)

ids <- basename(files) |> sub("_4bpmotif_gc\\.rds$", "", x = _)

message("Combining ", length(files), " samples from: ", motifdir)

tib_all <- map2_dfr(files, ids, function(f, id) {
  x <- readRDS(f)
  tb <- as_tibble(x)

  if (!("end_motif" %in% names(tb))) {
    stop("Missing 'end_motif' in file: ", f, " | columns: ", paste(names(tb), collapse = ", "))
  }

  tb %>% mutate(id = id)
}) %>%
  select(id, everything()) %>%
  mutate(end_motif = toupper(str_trim(end_motif)))

n_before <- n_distinct(tib_all$end_motif)

# -----------------------------
# ATGC-only + 4bp filtering
# -----------------------------
tib_atgc <- tib_all %>%
  filter(!is.na(end_motif)) %>%
  filter(str_length(end_motif) == 4) %>%
  filter(!grepl("[^ATGC]", end_motif))

n_after <- n_distinct(tib_atgc$end_motif)

message("Unique motifs BEFORE filter: ", n_before)
message("Unique motifs AFTER  filter: ", n_after)
message("Rows BEFORE: ", nrow(tib_all), " | Rows AFTER: ", nrow(tib_atgc))
message("ATGC 4-mers present (global): ", n_after, " (expected up to 256)")

# -----------------------------
# Write outputs
# -----------------------------
out_rds <- file.path(outdir, "endmotif_4bp_gc_summary_ATGC_only.rds")
out_tsv <- file.path(outdir, "endmotif_4bp_gc_summary_ATGC_only.tsv.gz")

saveRDS(tib_atgc, out_rds)
write_tsv(tib_atgc, out_tsv)

message("Wrote RDS: ", out_rds)
message("Wrote TSV: ", out_tsv)

dim(tib_atgc)
head(tib_atgc)
length(unique(tib_atgc$end_motif))   # should be <= 256
length(unique(tib_atgc$id))         # number of samples combined

```
