# Template provided by Brent
# Converts multiple .cns files to GISTIC-compatible segmentation format

library(tidyverse)

# ==== CONFIGURATION ====
cns_directory <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/GISTIC/cns_file"
output_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/GISTIC"

intermediate_seg <- file.path(output_dir, "gistic_input.seg")
final_output_seg <- file.path(output_dir, "gistic_input_fixed.seg")

cns_file_pattern <- "\\.bqsr\\.cns$"  # Match .bqsr.cns files
merge_similar_segments_flag <- FALSE
cn_diff_threshold <- 0.1

# ==== STEP 1: Convert .cns files to .seg ====

process_cns_file <- function(file_path) {
  sample_name <- basename(file_path) %>%
    str_remove("\\.bqsr\\.cns$")

  cat("Processing:", sample_name, "\n")

  tryCatch({
    data <- read_tsv(file_path, col_types = cols(
      chromosome = col_character(),
      start = col_double(),
      end = col_double(),
      gene = col_character(),
      log2 = col_double(),
      depth = col_double(),
      probes = col_double(),
      weight = col_double()
    ))

    gistic_data <- data %>%
      select(chromosome, start, end, probes, log2) %>%
      mutate(
        Sample = sample_name,
        Chromosome = str_remove(chromosome, "^chr"),
        `Start Position` = start,
        `End Position` = end,
        `Num Markers` = as.integer(probes),
        `Seg.CN` = log2
      ) %>%
      select(Sample, Chromosome, `Start Position`, `End Position`, `Num Markers`, `Seg.CN`)

    return(gistic_data)
  }, error = function(e) {
    cat("Error processing", sample_name, ":", e$message, "\n")
    return(NULL)
  })
}

convert_cns_to_seg <- function(cns_directory, output_file) {
  cns_files <- list.files(cns_directory, pattern = cns_file_pattern, full.names = TRUE)

  if (length(cns_files) == 0) stop("No .bqsr.cns files found in directory: ", cns_directory)

  cat("Found", length(cns_files), ".bqsr.cns files\n")

  all_segments <- map_dfr(cns_files, process_cns_file)

  if (nrow(all_segments) == 0) stop("No segments processed successfully.")

  all_segments <- all_segments %>%
    arrange(Sample, Chromosome, `Start Position`)

  write_tsv(all_segments, output_file, na = "")
  cat("Wrote initial segmentation file:", output_file, "\n")

  return(output_file)
}

# ==== STEP 2: Fix Overlapping Segments ====

fix_overlaps_per_sample_chr <- function(sample_chr_data) {
  sample_chr_data <- sample_chr_data %>% arrange(`Start Position`)
  if (nrow(sample_chr_data) > 1) {
    for (i in 1:(nrow(sample_chr_data) - 1)) {
      if (sample_chr_data$`End Position`[i] >= sample_chr_data$`Start Position`[i + 1]) {
        sample_chr_data$`End Position`[i] <- sample_chr_data$`Start Position`[i + 1] - 1
      }
    }
  }
  sample_chr_data %>% filter(`Start Position` < `End Position`)
}

merge_similar_segments <- function(sample_chr_data, cn_threshold = 0.1) {
  if (nrow(sample_chr_data) <= 1) return(sample_chr_data)

  sample_chr_data <- sample_chr_data %>%
    arrange(`Start Position`) %>%
    mutate(
      next_cn = lead(`Seg.CN`),
      next_start = lead(`Start Position`),
      cn_diff = abs(`Seg.CN` - next_cn),
      adjacent = (`End Position` + 1) == next_start,
      should_merge = adjacent & !is.na(cn_diff) & cn_diff <= cn_threshold,
      merge_group = cumsum(!should_merge | is.na(should_merge))
    ) %>%
    group_by(merge_group) %>%
    summarise(
      Sample = first(Sample),
      Chromosome = first(Chromosome),
      `Start Position` = min(`Start Position`),
      `End Position` = max(`End Position`),
      `Num Markers` = sum(`Num Markers`),
      `Seg.CN` = mean(`Seg.CN`, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(`Start Position`)

  return(sample_chr_data)
}

check_overlaps <- function(seg_file) {
  seg_data <- read_tsv(seg_file, col_types = cols(
    Sample = col_character(),
    Chromosome = col_character(),
    `Start Position` = col_double(),
    `End Position` = col_double(),
    `Num Markers` = col_double(),
    `Seg.CN` = col_double()
  ))

  overlap_count <- seg_data %>%
    group_by(Sample, Chromosome) %>%
    arrange(`Start Position`) %>%
    mutate(
      next_start = lead(`Start Position`, default = Inf),
      overlap = `End Position` >= next_start
    ) %>%
    filter(overlap) %>%
    nrow()

  return(overlap_count)
}

fix_segmentation_file <- function(input_file, output_file, merge_similar = FALSE, cn_threshold = 0.1) {
  seg_data <- read_tsv(input_file, col_types = cols(
    Sample = col_character(),
    Chromosome = col_character(),
    `Start Position` = col_double(),
    `End Position` = col_double(),
    `Num Markers` = col_double(),
    `Seg.CN` = col_double()
  ))

  seg_data <- seg_data %>%
    mutate(
      Chromosome = as.character(Chromosome),
      chr_sort = case_when(
        Chromosome == "X" ~ 23,
        Chromosome == "Y" ~ 24,
        TRUE ~ as.numeric(Chromosome)
      )
    ) %>%
    arrange(Sample, chr_sort, `Start Position`) %>%
    select(-chr_sort)

  # Fix overlaps
  cat("Fixing overlaps...\n")
  fixed_data <- seg_data %>%
    group_by(Sample, Chromosome) %>%
    group_split() %>%
    map_dfr(fix_overlaps_per_sample_chr)

  # Optional merge
  if (merge_similar) {
    cat("Merging similar segments...\n")
    fixed_data <- fixed_data %>%
      group_by(Sample, Chromosome) %>%
      group_split() %>%
      map_dfr(~merge_similar_segments(.x, cn_threshold))
  }

  fixed_data <- fixed_data %>%
    mutate(
      chr_sort = case_when(
        Chromosome == "X" ~ 23,
        Chromosome == "Y" ~ 24,
        TRUE ~ as.numeric(Chromosome)
      )
    ) %>%
    arrange(Sample, chr_sort, `Start Position`) %>%
    select(-chr_sort)

  write_tsv(fixed_data, output_file, na = "")
  cat("Wrote fixed segmentation file:", output_file, "\n")
}

# ==== RUNNING THE PIPELINE ====
cat("=== CNVkit to GISTIC Pipeline ===\n")

# Step 1: Convert .cns ? .seg
cat("\n[1] Converting .cns files to GISTIC .seg format...\n")
convert_cns_to_seg(cns_directory, intermediate_seg)

# Step 2: Check for overlaps
cat("\n[2] Checking for overlapping segments...\n")
overlaps <- check_overlaps(intermediate_seg)

if (overlaps > 0) {
  cat("Overlaps found:", overlaps, "\n")
  cat("\n[3] Fixing overlaps and writing final file...\n")
  fix_segmentation_file(
    input_file = intermediate_seg,
    output_file = final_output_seg,
    merge_similar = merge_similar_segments_flag,
    cn_threshold = cn_diff_threshold
  )

  # Recheck
  final_overlaps <- check_overlaps(final_output_seg)
  if (final_overlaps == 0) {
    cat("Final segmentation file is clean: ", final_output_seg, "\n")
  } else {
    cat("Overlaps still exist after fixing!\n")
  }

} else {
  cat("No overlaps found. Using initial file for GISTIC:", intermediate_seg, "\n")
}

cat("\nPipeline complete!\n")
