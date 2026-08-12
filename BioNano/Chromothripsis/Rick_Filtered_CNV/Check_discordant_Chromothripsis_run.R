# ============================================================
# Check discordant ShatterSeek calls: Ours vs Osama
# ============================================================

OUT_DIR <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Chromothripsis/Osama_Script"

ss <- read.delim(
  file.path(OUT_DIR, "ShatterSeek_df.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Make chromosome format consistent
ss$chrom_clean <- gsub("^chr", "", as.character(ss$chrom))


# ------------------------------------------------------------
# Discordant loci identified from heatmap comparison
# ------------------------------------------------------------

targets <- data.frame(
  sample = c(
    "IPM-BO-046",
    "IPM-BO-053",
    "TP05-232",
    "IPM-BO-084",
    "IPM-BO-084",
    "TP04-104",
    "TP11-2004",
    "TP11-2004",
    "TP11-2004",
    "TP11-2004"
  ),

  chrom_clean = c(
    "10",
    "4",
    "10",
    "4",
    "7",
    "5",
    "12",
    "17",
    "20",
    "21"
  ),

  Osama_call = c(
    "Low",
    "Low",
    "Low",
    "High",
    "High",
    "Low",
    "High",
    "High",
    "High",
    "High"
  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Merge target loci with our ShatterSeek result
# ------------------------------------------------------------

x <- merge(
  targets,
  ss,
  by = c("sample", "chrom_clean"),
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# Recalculate Osama High criteria
# ------------------------------------------------------------

x$H_cluster <- !is.na(x$clusterSize_including_TRA) &
  x$clusterSize_including_TRA >= 12

x$H_osc <- !is.na(x$max_number_oscillating_CN_segments_2_states) &
  x$max_number_oscillating_CN_segments_2_states >= 4

x$H_frag <- !is.na(x$pval_fragment_joins) &
  x$pval_fragment_joins <= 0.001

x$H_enrich <- (
  (!is.na(x$chr_breakpoint_enrichment) &
     x$chr_breakpoint_enrichment <= 0.001) |
  (!is.na(x$pval_exp_cluster) &
     x$pval_exp_cluster <= 0.001)
)

x$High_pass <- (
  x$H_cluster &
  x$H_osc &
  x$H_frag &
  x$H_enrich
)


# ------------------------------------------------------------
# Recalculate Osama Low criteria
# ------------------------------------------------------------

x$L_intra <- !is.na(x$total_IntraSVs) &
  x$total_IntraSVs >= 8

x$L_osc <- !is.na(x$max_number_oscillating_CN_segments_2_states) &
  x$max_number_oscillating_CN_segments_2_states >= 2

x$L_frag <- !is.na(x$pval_fragment_joins) &
  x$pval_fragment_joins <= 0.05

x$L_enrich <- (
  (!is.na(x$chr_breakpoint_enrichment) &
     x$chr_breakpoint_enrichment <= 0.05) |
  (!is.na(x$pval_exp_cluster) &
     x$pval_exp_cluster <= 0.05)
)

x$Low_pass <- (
  x$L_intra &
  x$L_osc &
  x$L_frag &
  x$L_enrich
)


# ------------------------------------------------------------
# What our criteria predict
# ------------------------------------------------------------

x$Recalculated_call <- ifelse(
  x$High_pass,
  "High",
  ifelse(
    x$Low_pass,
    "Low",
    "ns"
  )
)


# Original confidence from our ShatterSeek_df
x$Our_call <- x$confidence
x$Our_call[
  is.na(x$Our_call) |
  x$Our_call == ""
] <- "ns"


# ------------------------------------------------------------
# Compact output
# ------------------------------------------------------------

result <- x[, c(
  "sample",
  "chrom_clean",
  "Osama_call",
  "Our_call",
  "Recalculated_call",

  "clusterSize_including_TRA",
  "total_IntraSVs",
  "max_number_oscillating_CN_segments_2_states",

  "pval_fragment_joins",
  "chr_breakpoint_enrichment",
  "pval_exp_cluster",

  "H_cluster",
  "H_osc",
  "H_frag",
  "H_enrich",

  "L_intra",
  "L_osc",
  "L_frag",
  "L_enrich"
)]


# Arrange in target order
result$order <- match(
  paste(result$sample, result$chrom_clean),
  paste(targets$sample, targets$chrom_clean)
)

result <- result[
  order(result$order),
]

result$order <- NULL


# ------------------------------------------------------------
# Print everything clearly
# ------------------------------------------------------------

options(
  width = 250,
  digits = 6
)

print(
  result,
  row.names = FALSE
)


# ------------------------------------------------------------
# Save it
# ------------------------------------------------------------

write.table(
  result,
  file = file.path(
    OUT_DIR,
    "Osama_vs_Our_discordant_ShatterSeek_metrics.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)