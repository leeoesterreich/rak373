```r
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
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

required <- c("file", "outdir", "plotdir", "statdir")
missing <- setdiff(required, names(opt))
if (length(missing) > 0) stop("Missing required args: ", paste(missing, collapse = ", "))

file     <- opt$file
outdir   <- opt$outdir
plotdir  <- opt$plotdir
statdir  <- opt$statdir

nmer     <- if (!is.null(opt$nmer)) as.integer(opt$nmer) else 4L
min_frag <- if (!is.null(opt$min_frag)) as.integer(opt$min_frag) else 90L
max_frag <- if (!is.null(opt$max_frag)) as.integer(opt$max_frag) else 650L

if (!file.exists(file)) stop("Input file not found: ", file)
dir.create(outdir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)
dir.create(statdir, recursive = TRUE, showWarnings = FALSE)

id <- basename(file)
id <- sub(paste0("_", nmer, "bp_motif\\.bed$"), "", id)
id <- sub("_4bp_motif\\.bed$", "", id)  # fallback

# ensure plot subdir exists
dir.create(file.path(plotdir, id), recursive = TRUE, showWarnings = FALSE)

out_file  <- file.path(outdir,  paste0(id, "_", nmer, "bpmotif_gc.rds"))
stat_file <- file.path(statdir, paste0(id, "_", nmer, "bpmotif_stat.rds"))

message("Processing sample: ", id)
message("Input: ", file)

# -----------------------------
# Functions (GC correction)
# -----------------------------
seqlast <- function(from, to, by) {
  vec <- seq(from = from, to = to, by = by)
  if (tail(vec, 1) != to) c(vec, to) else vec
}

gc.correct <- function(coverage, bias) {
  # coverage, bias are numeric vectors (same length), bias in [0,1]
  # Returns numeric vector of corrected coverage
  ok <- is.finite(coverage) & is.finite(bias)
  coverage <- coverage[ok]
  bias <- bias[ok]

  if (length(coverage) < 10) return(rep(NA_real_, length(ok)))  # too few points

  i <- seqlast(min(bias, na.rm = TRUE), max(bias, na.rm = TRUE), by = 0.01)

  # loess can fail for small/degenerate data; guard it
  corrected <- tryCatch({
    coverage.trend <- loess(coverage ~ bias, span = 0.75)
    coverage.model <- loess(predict(coverage.trend, i) ~ i, span = 0.75)
    coverage.pred  <- predict(coverage.model, bias)
    coverage - coverage.pred + median(coverage, na.rm = TRUE)
  }, error = function(e) {
    rep(NA_real_, length(coverage))
  })

  # re-inflate to original length with NAs where excluded
  out <- rep(NA_real_, length(ok))
  out[ok] <- corrected
  out
}

# -----------------------------
# Read motif bed
# Expected columns in your merged file:
# 1 chr, 2 start, 3 end, 4 strand, 5 frag_size, 6 gc, 7 end_motif
# -----------------------------
df <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")

if (ncol(df) < 7) {
  msg <- paste0("Skipping ", basename(file), ": expected >=7 columns, got ", ncol(df))
  warning(msg)
  saveRDS(list(sample=id, input_file=file, ncol=ncol(df), note=msg), stat_file)
  quit(save="no", status=0)
}

# Ensure standard column names exist
colnames(df)[1:7] <- paste0("V", 1:7)

df <- df %>%
  dplyr::transmute(
    chr       = V1,
    start     = as.integer(V2),
    end       = as.integer(V3),
    strand    = V4,
    frag_size = as.integer(V5),
    gc        = as.numeric(V6),
    end_motif = toupper(V7)
  )

# -----------------------------
# Filtering + nucleosome fraction
# -----------------------------
n_in <- nrow(df)

df_filt <- df %>%
  filter(!is.na(frag_size), !is.na(gc)) %>%
  filter(frag_size >= min_frag, frag_size <= max_frag) %>%
  mutate(
    frac = case_when(
      frag_size >=  90 & frag_size <= 250 ~ "mono",
      frag_size >= 251 & frag_size <= 450 ~ "di",
      frag_size >= 451 & frag_size <= 650 ~ "tri",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(frac)) %>%
  mutate(gc = round(pmax(pmin(gc, 1), 0), 2))  # clamp into [0,1] then round

n_after_frag <- nrow(df_filt)

# If nothing remains, still save stats for debugging
if (n_after_frag == 0) {
  df_stats <- list(
    sample = id,
    input_rows = n_in,
    rows_after_frag_filter = n_after_frag,
    note = "No rows after fragment-size filtering"
  )
  saveRDS(df_stats, stat_file)
  stop("No motifs remaining after fragment-size filter for sample: ", id)
}

# -----------------------------
# Count motifs by (end_motif, frac, gc)
# -----------------------------
df_motif <- df_filt %>%
  select(end_motif, frag_size, gc, frac) %>%
  group_by(end_motif, frac, gc) %>%
  summarise(count = n(), .groups = "drop")

# -----------------------------
# GC scalar per frac/gc with 5-95 percentile trimming
# -----------------------------
df_gc <- df_motif %>%
  group_by(frac, gc) %>%
  summarise(frac.count = sum(count), .groups = "drop") %>%
  arrange(frac, gc) %>%
  group_by(frac) %>%
  mutate(
    frac.total = sum(frac.count),
    cum_count  = cumsum(frac.count),
    percentile = (cum_count / frac.total) * 100
  ) %>%
  ungroup() %>%
  filter(percentile >= 5.00 & percentile <= 95.00)

# apply loess correction within each fraction
df_gc <- df_gc %>%
  group_by(frac) %>%
  mutate(
    frac.count.corrected = gc.correct(frac.count, gc)
  ) %>%
  ungroup() %>%
  group_by(frac) %>%
  mutate(
    total.corrected = sum(frac.count.corrected, na.rm = TRUE),
    gc.scale = frac.count.corrected / frac.count
  ) %>%
  ungroup() %>%
  select(frac, gc, frac.count, frac.count.corrected, gc.scale)

# if loess failed badly, protect against all-NA
if (all(!is.finite(df_gc$gc.scale))) {
  warning("GC correction produced all NA scales; writing uncorrected counts only for: ", id)
  df_gc$gc.scale <- 1.0
  df_gc$frac.count.corrected <- df_gc$frac.count
}

# -----------------------------
# Apply scalar to each motif
# -----------------------------
df_out <- df_motif %>%
  inner_join(df_gc, by = c("frac", "gc")) %>%
  mutate(count.corrected = gc.scale * count) %>%
  group_by(frac) %>%
  mutate(
    frac.total = sum(count),
    frac.total.corrected = sum(count.corrected)
  ) %>%
  ungroup() %>%
  mutate(
    total.count = sum(count),
    total.corrected = sum(count.corrected)
  )

# -----------------------------
# Stats (minimal but useful)
# -----------------------------
df_stats <- list(
  sample = id,
  input_rows = n_in,
  rows_after_frag_filter = n_after_frag,
  unique_motifs = n_distinct(df_out$end_motif),
  gc_bins_used = n_distinct(df_out$gc),
  frac_levels = sort(unique(df_out$frac))
)

# -----------------------------
# Plots (overall GC bias per fraction)
# -----------------------------
# Make a plot per fraction: uncorrected vs corrected coverage trend
fracs <- sort(unique(df_gc$frac))

for (flevel in fracs) {
  dsub <- df_gc %>% filter(frac == flevel)

  p <- ggplot(dsub, aes(x = gc, y = frac.count)) +
    geom_line() +
    geom_line(aes(y = frac.count.corrected), color = "red") +
    labs(
      title = paste0(id, " - ", flevel, " (", nmer, "-mer)"),
      x = "GC",
      y = "Counts",
      color = "Series"
    )

  out_png <- file.path(plotdir, id, paste0(id, "_", flevel, "_", nmer, "bpmotif_gc.png"))
  ggsave(filename = out_png, plot = p, width = 7, height = 5)
}

# -----------------------------
# Save outputs
# -----------------------------
saveRDS(df_out, out_file)
saveRDS(df_stats, stat_file)

message("Wrote: ", out_file)
message("Wrote: ", stat_file)
message("Plots: ", file.path(plotdir, id))
```
