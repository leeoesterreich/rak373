#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(facetsSuite)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript extract_cna_from_rds.R <input_rds> <output_tsv>")
}

input_rds <- args[1]
output_tsv <- args[2]

# Load RDS
out <- readRDS(input_rds)

# Set NA lcn.em to 0 and compute major/minor CN
out$segs$lcn.em[is.na(out$segs$lcn.em)] <- 0
minor_cn <- out$segs$lcn.em
total_cn <- out$segs$tcn.em
major_cn <- total_cn - minor_cn

# Make data.frame
cn_df <- data.frame(
  Chromosome = out$segs$chrom,
  Start = out$segs$start,
  End = out$segs$end,
  nMajor = major_cn,
  nMinor = minor_cn
)

# Write output
write_tsv(cn_df, output_tsv)

