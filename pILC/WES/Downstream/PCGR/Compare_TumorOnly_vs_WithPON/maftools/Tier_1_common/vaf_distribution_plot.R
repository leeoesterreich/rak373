# ============================================================
# VAF distribution across all nonsilent somatic mutations
# Shows how many mutations fall below the 1% / 3% / 5% cutoffs
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)
graphics.off()

# ------------------------------------------------------------
# Paths (same inputs as the oncoplot scripts)
# ------------------------------------------------------------

combined_maf_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Tier_1_common_all_samples.combined.maf"

vaf_audit_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/Combined_MAF/Full_MAF_nonsilent_variant_PCGR_VAF_match_audit_2026-09-23.tsv"

out_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/maftools/Tier_1_common/plots/VAF_distribution"
table_dir <- file.path(out_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Cutoffs to mark on the plot (percent)
VAF_CUTOFFS <- c(1, 3, 5)

nonsilent_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Splice_Site", "Splice_Region", "Translation_Start_Site"
)

# ------------------------------------------------------------
# Helper functions (identical to the oncoplot scripts, so
# Match_key agrees with the VAF audit table)
# ------------------------------------------------------------

vaf_clean_sample <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("(_T|-T|\\.T)$", "", x)
  x[x %in% c("", ".", "NA", "N/A")] <- NA_character_
  x
}

vaf_clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", ".", "NA", "N/A")] <- NA_character_
  x
}

edit_key <- function(chr, pos, ref, alt, maf_format = FALSE) {
  chr <- sub("^CHR", "", toupper(trimws(as.character(chr))))
  pos <- suppressWarnings(as.integer(pos))
  ref <- toupper(trimws(as.character(ref)))
  alt <- toupper(trimws(as.character(alt)))

  if (anyNA(c(chr, pos, ref, alt)) ||
      !nzchar(chr) || pos < 1L ||
      !grepl("^([ACGTN]+|-)$", ref) ||
      !grepl("^([ACGTN]+|-)$", alt)) {
    return(NA_character_)
  }

  if (maf_format && ref == "-") pos <- pos + 1L
  if (ref == "-") ref <- ""
  if (alt == "-") alt <- ""

  while (nzchar(ref) && nzchar(alt) &&
         substr(ref, nchar(ref), nchar(ref)) == substr(alt, nchar(alt), nchar(alt))) {
    ref <- substr(ref, 1, nchar(ref) - 1L)
    alt <- substr(alt, 1, nchar(alt) - 1L)
  }
  while (nzchar(ref) && nzchar(alt) && substr(ref, 1, 1) == substr(alt, 1, 1)) {
    ref <- substring(ref, 2)
    alt <- substring(alt, 2)
    pos <- pos + 1L
  }
  if (!nzchar(ref) && !nzchar(alt)) return(NA_character_)

  paste(chr, pos, ifelse(nzchar(ref), ref, "-"), ifelse(nzchar(alt), alt, "-"), sep = ":")
}

# ------------------------------------------------------------
# Read MAF, rebuild Match_key, attach VAF (same logic as the
# oncoplot scripts, but WITHOUT applying a VAF threshold here —
# we need the full distribution, not a filtered subset)
# ------------------------------------------------------------

cat("Reading combined MAF:\n", combined_maf_file, "\n\n")

maf_full <- fread(
  combined_maf_file, sep = "\t", header = TRUE, data.table = TRUE,
  na.strings = c("", ".", "NA", "N/A")
)

maf_full[, Variant_key := mapply(
  edit_key, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
  MoreArgs = list(maf_format = TRUE), USE.NAMES = FALSE
)]
maf_full[, Match_key := paste(
  vaf_clean_sample(Tumor_Sample_Barcode), vaf_clean_gene(Hugo_Symbol), Variant_key, sep = "|"
)]

vaf_audit <- fread(vaf_audit_file, sep = "\t", header = TRUE, data.table = TRUE)
vaf_by_key <- unique(
  vaf_audit[, .(Match_key, VAF_TUMOR_fraction = suppressWarnings(as.numeric(VAF_TUMOR_fraction)))],
  by = "Match_key"
)
maf_full <- merge(maf_full, vaf_by_key, by = "Match_key", all.x = TRUE, sort = FALSE)

maf_full[, Is_nonsilent := Variant_Classification %in% nonsilent_classes]

# Nonsilent mutations only, from here on
nonsilent_dt <- maf_full[Is_nonsilent == TRUE]
nonsilent_dt[, VAF_percent := 100 * VAF_TUMOR_fraction]

n_total_nonsilent <- nrow(nonsilent_dt)
n_missing_vaf <- nonsilent_dt[is.na(VAF_percent), .N]
vaf_dt <- nonsilent_dt[!is.na(VAF_percent)]

cat("Total nonsilent mutations:", n_total_nonsilent, "\n")
cat("Missing a matched PCGR VAF (excluded from the plot):", n_missing_vaf, "\n")
cat("Nonsilent mutations with a VAF value (plotted):", nrow(vaf_dt), "\n\n")

# ------------------------------------------------------------
# Counts at each cutoff (for the PPT text / table)
# ------------------------------------------------------------

cutoff_summary <- rbindlist(lapply(VAF_CUTOFFS, function(cutoff) {
  data.table(
    VAF_cutoff_percent = cutoff,
    Mutations_retained_at_or_above_cutoff = vaf_dt[VAF_percent >= cutoff, .N],
    Mutations_excluded_below_cutoff = vaf_dt[VAF_percent < cutoff, .N]
  )
}))

cutoff_summary <- rbind(
  data.table(
    VAF_cutoff_percent = 0,
    Mutations_retained_at_or_above_cutoff = nrow(vaf_dt),
    Mutations_excluded_below_cutoff = 0L
  ),
  cutoff_summary
)

cat("Mutation counts by VAF cutoff:\n")
print(cutoff_summary)

fwrite(
  cutoff_summary,
  file.path(table_dir, "VAF_cutoff_mutation_counts.tsv"),
  sep = "\t", quote = FALSE, na = ""
)

# Bin the low-VAF range finely so reviewers can see where mass sits
range_bins_dt <- data.table(
  VAF_range = c("<1%", "1-3%", "3-5%", ">=5%"),
  N_mutations = c(
    vaf_dt[VAF_percent < 1, .N],
    vaf_dt[VAF_percent >= 1 & VAF_percent < 3, .N],
    vaf_dt[VAF_percent >= 3 & VAF_percent < 5, .N],
    vaf_dt[VAF_percent >= 5, .N]
  )
)
fwrite(
  range_bins_dt,
  file.path(table_dir, "VAF_range_bin_counts.tsv"),
  sep = "\t", quote = FALSE, na = ""
)
cat("\nMutations by VAF range bin:\n")
print(range_bins_dt)

# ------------------------------------------------------------
# Histogram with cutoff lines (publication style)
#
# Key changes from the earlier version:
# - center = 0 (instead of boundary = 0) so each axis number sits in
#   the MIDDLE of its bar, not on the edge between two bars.
# - The "1% / 3% / 5%" tags on the dashed lines are now short labels
#   only; the exclusion counts moved into a figure caption below the
#   plot, so nothing overlaps regardless of how close the cutoffs are.
# - Comma-formatted axis numbers, bold titles, cleaner theme.
# ------------------------------------------------------------

cutoff_caption <- paste0(
  "Excluded below 1%: ", format(cutoff_summary$Mutations_excluded_below_cutoff[2], big.mark = ","),
  "    |    below 3%: ", format(cutoff_summary$Mutations_excluded_below_cutoff[3], big.mark = ","),
  "    |    below 5%: ", format(cutoff_summary$Mutations_excluded_below_cutoff[4], big.mark = ",")
)

subtitle_base <- paste0(
  "n = ", format(nrow(vaf_dt), big.mark = ","), " mutations with a matched VAF  (",
  format(n_missing_vaf, big.mark = ","), " excluded for missing VAF)"
)

# Shared look for both the full-range and zoomed plots
p_base <- ggplot(vaf_dt, aes(x = VAF_percent)) +
  geom_histogram(
    binwidth = 1, center = 0,
    fill = "#3B6FA0", color = "white", linewidth = 0.15
  ) +
  geom_vline(xintercept = VAF_CUTOFFS, linetype = "dashed", color = "#C0392B", linewidth = 0.6) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "VAF distribution across nonsilent somatic mutations",
    subtitle = subtitle_base,
    x = "Variant allele frequency (%)",
    y = "Number of mutations",
    caption = cutoff_caption
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey30", size = 11),
    plot.caption = element_text(color = "#C0392B", size = 10, hjust = 0, margin = margin(t = 8)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(clip = "off")

# ------------------------------------------------------------
# Full-range version (0-100%)
# ------------------------------------------------------------

p_full <- p_base +
  scale_x_continuous(limits = c(-0.5, 100.5), breaks = seq(0, 100, by = 10), expand = c(0, 0))

ggsave(file.path(out_dir, "VAF_distribution_nonsilent_mutations.pdf"), p_full, width = 9, height = 5.5)
ggsave(file.path(out_dir, "VAF_distribution_nonsilent_mutations.png"), p_full, width = 9, height = 5.5, dpi = 300)

# ------------------------------------------------------------
# Zoomed-in version (0-15%): enough room per bar to add short
# "1% / 3% / 5%" tags above each dashed line without overlap.
# ------------------------------------------------------------

# Compute the tallest bar within the 0-15% zoom window directly
# (avoids hist()'s "breaks do not span range of x" error, since
# vaf_dt still holds the full 0-100% range).
max_label_y <- vaf_dt[
  VAF_percent >= -0.5 & VAF_percent < 15.5,
  .N,
  by = .(bin = floor(VAF_percent + 0.5))   # same binning as center = 0, binwidth = 1
][, max(N)]

p_zoom <- p_base +
  scale_x_continuous(limits = c(-0.5, 15.5), breaks = seq(0, 15, by = 1), expand = c(0, 0)) +
  annotate(
    "text",
    x = VAF_CUTOFFS,
    y = max_label_y * 1.08,
    label = paste0(VAF_CUTOFFS, "%"),
    color = "#C0392B",
    fontface = "bold",
    size = 4,
    vjust = 0
  ) +
  labs(subtitle = paste0(subtitle_base, "  |  zoomed to 0-15% VAF"))

ggsave(file.path(out_dir, "VAF_distribution_nonsilent_mutations_zoom0to15.pdf"), p_zoom, width = 9, height = 5.5)
ggsave(file.path(out_dir, "VAF_distribution_nonsilent_mutations_zoom0to15.png"), p_zoom, width = 9, height = 5.5, dpi = 300)

cat("\nSaved:\n")
cat(" ", file.path(out_dir, "VAF_distribution_nonsilent_mutations.pdf/.png"), "\n")
cat(" ", file.path(out_dir, "VAF_distribution_nonsilent_mutations_zoom0to15.pdf/.png"), "\n")
cat(" ", file.path(table_dir, "VAF_cutoff_mutation_counts.tsv"), "\n")
cat(" ", file.path(table_dir, "VAF_range_bin_counts.tsv"), "\n")
