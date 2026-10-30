# ==============================================================================
# Count fusions per tool (Arriba, FusionCatcher, STAR-Fusion) from nf-core/rnafusion
#
# Steps for each tool:
#   1) all calls                 -> raw output of the tool
#   2) tool threshold            -> Arriba high/medium | FusionCatcher label filter | STAR-Fusion FFPM >= 0.1
#   3) after artifact removal    -> likely artifacts removed from the calls kept in step 2
# Counts reported: total calls, unique fusions (A--B == B--A), samples with >= 1 call
#
# How the 5'/3' orientation is decided (section 13; details in fusion_orientation_details.tsv):
#   1) Take the order reported by each caller (gene_a = 5' partner, gene_b = 3' partner):
#        Arriba         gene1 (5')                         / gene2 (3')
#        FusionCatcher  Gene_1_symbol(5end_fusion_partner) / Gene_2_symbol(3end_fusion_partner)
#        STAR-Fusion    LeftGene (5', by convention)       / RightGene (3')
#      The fusion name (A--B) is sorted alphabetically, so it carries NO orientation.
#   2) For each sample + gene pair, compare the 5'>3' order across callers
#      (orientation_consistent = TRUE when all callers give the same order).
#   3) Final orientation = the order supported by the MOST callers; if tied, the order with
#      the most calls. A remaining tie is flagged (orientation_tie) for manual review.
#   4) The other order is kept (orientation_other) with its callers and breakpoints,
#      so a possible reciprocal fusion or strand misprediction can still be checked.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(VennDiagram)   # Venn diagrams (section 14)
  library(grid)
})

# ---- 1. Paths ----------------------------------------------------------------
BASE  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/nfcore_rnafusion"
RES   <- file.path(BASE, "results")                    # nf-core output folder
SHEET <- file.path(BASE, "samplesheet_fusion.csv")     # used to find samples without output
OUT   <- file.path(BASE, "fusion_counts")              # all tables and figures go here
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Tool thresholds (same as the previous analysis) -----------------------
ARRIBA_CONF <- c("high", "medium")    # Arriba confidence levels kept
STARF_FFPM  <- 0.1                    # STAR-Fusion: FFPM >= 0.1

# FusionCatcher: a call is dropped if ANY of these labels is in its Fusion_description.
# Labels are FusionCatcher's own false-positive flags (checked against its manual):
#   healthy-sample evidence : banned, hpa, bodymap2, 1000genomes, gtex, healthy, cortex, non_tumor_cells,
#                             non_cancer_tissues, cacg, conjoing
#   sequence similarity     : paralogs, duplicates, similar_reads, similar_symbols, pair_pseudo_genes
#   gene type / position    : mt, ribosomal, pseudogene, readthrough, adjacent
#   gene overlap            : *_fully_overlapping, *_partially_overlapping (Ensembl, RefSeq, UCSC)
# Not used: mX labels (read counts) and labels I could not verify (ambiguous, rt_circ_rna, mitelman, ccle, lncrna).
# Compare the resulting FusionCatcher counts with your earlier run (721 calls / 554 unique) and edit if needed;
# fc_label_tally.csv lists every label that occurs in your data.
FC_DROP_LABELS <- c("banned", "hpa", "bodymap2", "1000genomes", "gtex", "healthy", "cortex",
                    "non_tumor_cells", "non_cancer_tissues", "cacg", "conjoing",
                    "paralogs", "duplicates", "similar_reads", "similar_symbols", "pair_pseudo_genes",
                    "mt", "ribosomal", "pseudogene", "readthrough", "adjacent",
                    "ensembl_fully_overlapping", "refseq_fully_overlapping", "ucsc_fully_overlapping",
                    "ensembl_partially_overlapping", "refseq_partially_overlapping", "ucsc_partially_overlapping")

# ---- 3. Artifact rules (edit freely; each rule is reported separately) --------
# 3a. gene-name patterns: mitochondrial, rRNA/snRNA/snoRNA, miRNA, lncRNA (LINC, -AS1),
#     clone-based names (RP11-..., CTD-..., AC012345.1 ...)
#     HLA and immunoglobulin / T-cell receptor genes are deliberately KEPT (no evidence found for removing them).
ARTIFACT_GENE_REGEX <- paste0(
  "^MT-|^RNA5S|^RN7S|^RNU[0-9]|^SNOR|^MIR[0-9]|^LINC[0-9]|-AS[0-9]$|",
  "^(RP[0-9]+|CT[ABCD]|XXbac)-|^(AC|AL|AP|AF|AJ|BX|CR|FP|Z)[0-9]{5,6}\\.[0-9]+$"
)
# 3b. partner without a gene symbol (bare Ensembl ID such as ENSG00000282885)
DROP_ENSG_ONLY <- TRUE
# 3c. adjacent / read-through-like events: both breakpoints on the same chromosome and closer than this
READTHROUGH_BP <- 100000
# 3d. STAR-Fusion annotation tags that mark likely artifacts (checked against annots, case-insensitive).
#     Check the tags that really occur in starfusion_annots_tally.csv and edit this vector.
STARF_DROP_ANNOTS <- c("GTEx_recurrent", "BodyMap", "Greger_Normal", "Babiceanu_Normal", "ConjoinG", "NEIGHBORS")

# ---- 4. Helper functions -----------------------------------------------------
# keep one clean gene symbol: "ELF1(41522),WBP4(671)" -> "ELF1" ; "TP53^ENSG0000..." -> "TP53"
clean_gene <- function(x) {
  x <- sub(",.*$", "", as.character(x))
  x <- sub("\\^.*$", "", x)
  x <- sub("\\(.*\\)$", "", x)
  x <- trimws(x)
  x[x %in% c("", ".")] <- NA_character_
  x
}

# "chr15:39600845" or "14:49862595:-" -> chromosome (without "chr") and position
parse_bp <- function(x) {
  p <- tstrsplit(as.character(x), ":", fixed = TRUE, keep = 1:2)
  list(chr = sub("^chr", "", p[[1]]), pos = suppressWarnings(as.numeric(p[[2]])))
}

check_cols <- function(dt, need, f) {
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing column(s) ", paste(miss, collapse = ", "), " in ", f)
}

# ---- 5. Readers: one function per tool, each returns a standard table --------
read_arriba <- function(f) {
  dt <- fread(f, sep = "\t", header = TRUE, colClasses = "character", fill = TRUE, showProgress = FALSE)
  if (nrow(dt) == 0) return(NULL)
  setnames(dt, sub("^#", "", names(dt)))                      # header starts with "#gene1"
  check_cols(dt, c("gene1", "gene2", "breakpoint1", "breakpoint2", "confidence"), f)
  a <- parse_bp(dt$breakpoint1); b <- parse_bp(dt$breakpoint2)
  data.table(gene_a_raw = dt$gene1, gene_b_raw = dt$gene2,
             chr_a = a$chr, pos_a = a$pos, chr_b = b$chr, pos_b = b$pos,
             pass_tool = tolower(dt$confidence) %in% ARRIBA_CONF,
             r_tag = FALSE)
}

read_fusioncatcher <- function(f) {
  dt <- fread(f, sep = "\t", header = TRUE, colClasses = "character", fill = TRUE, showProgress = FALSE)
  if (nrow(dt) == 0) return(NULL)
  g1 <- "Gene_1_symbol(5end_fusion_partner)"; g2 <- "Gene_2_symbol(3end_fusion_partner)"
  p1 <- "Fusion_point_for_gene_1(5end_fusion_partner)"; p2 <- "Fusion_point_for_gene_2(3end_fusion_partner)"
  check_cols(dt, c(g1, g2, p1, p2, "Fusion_description"), f)
  a <- parse_bp(dt[[p1]]); b <- parse_bp(dt[[p2]])
  labels <- strsplit(tolower(dt$Fusion_description), ",", fixed = TRUE)
  data.table(gene_a_raw = dt[[g1]], gene_b_raw = dt[[g2]],
             chr_a = a$chr, pos_a = a$pos, chr_b = b$chr, pos_b = b$pos,
             pass_tool = !vapply(labels, function(l) any(trimws(l) %in% FC_DROP_LABELS), logical(1)),
             r_tag = FALSE,
             labels_raw = dt$Fusion_description)
}

read_starfusion <- function(f) {
  dt <- fread(f, sep = "\t", header = TRUE, colClasses = "character", fill = TRUE, showProgress = FALSE)
  if (nrow(dt) == 0) return(NULL)                              # header only = no fusions in this sample
  setnames(dt, sub("^#", "", names(dt)))
  check_cols(dt, c("LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint", "FFPM", "annots"), f)
  a <- parse_bp(dt$LeftBreakpoint); b <- parse_bp(dt$RightBreakpoint)
  tag_pat <- paste(STARF_DROP_ANNOTS, collapse = "|")
  data.table(gene_a_raw = dt$LeftGene, gene_b_raw = dt$RightGene,
             chr_a = a$chr, pos_a = a$pos, chr_b = b$chr, pos_b = b$pos,
             pass_tool = suppressWarnings(as.numeric(dt$FFPM)) >= STARF_FFPM,
             r_tag = grepl(tag_pat, dt$annots, ignore.case = TRUE),
             annots_raw = dt$annots)
}

# list the per-sample files of one tool, read them, and return calls + sample names found
run_tool <- function(tool, folder, pattern, reader) {
  files   <- list.files(file.path(RES, folder), pattern = pattern, full.names = TRUE)
  samples <- sub(pattern, "", basename(files))
  calls   <- rbindlist(Map(function(f, s) {
    d <- reader(f)
    if (!is.null(d)) d[, sample := s]
    d
  }, files, samples), fill = TRUE)
  if (nrow(calls) > 0) calls[, tool := tool]
  list(calls = calls, samples = samples)
}

# ---- 6. Read the three tools --------------------------------------------------
res_list <- list(
  Arriba        = run_tool("Arriba",        "arriba",        "\\.arriba\\.fusions\\.tsv$",              read_arriba),
  FusionCatcher = run_tool("FusionCatcher", "fusioncatcher", "\\.fusion-genes\\.txt$",                  read_fusioncatcher),
  `STAR-Fusion` = run_tool("STAR-Fusion",   "starfusion",    "\\.starfusion\\.abridged\\.tsv$",         read_starfusion)
)

calls <- rbindlist(lapply(res_list, `[[`, "calls"), fill = TRUE)

# samples expected (samplesheet) vs samples with an output file, per tool
expected <- fread(SHEET)$sample
for (tl in names(res_list)) {
  miss <- setdiff(expected, res_list[[tl]]$samples)
  message(tl, ": output files for ", length(res_list[[tl]]$samples), " of ", length(expected),
          " samples", if (length(miss)) paste0("  | MISSING: ", paste(miss, collapse = ", ")) else "")
}

# ---- 7. Clean gene names, build the fusion key, flag artifacts ---------------
calls[, gene_a := clean_gene(gene_a_raw)]
calls[, gene_b := clean_gene(gene_b_raw)]
# orientation-independent key: A--B and B--A are the same fusion
calls[, fusion := fifelse(is.na(gene_a) | is.na(gene_b), NA_character_,
                          fifelse(gene_a <= gene_b, paste0(gene_a, "--", gene_b),
                                                    paste0(gene_b, "--", gene_a)))]
calls <- calls[!is.na(fusion)]

calls[, r_name     := grepl(ARTIFACT_GENE_REGEX, gene_a) | grepl(ARTIFACT_GENE_REGEX, gene_b)]
calls[, r_ensg     := DROP_ENSG_ONLY & (grepl("^ENSG[0-9]+", gene_a) | grepl("^ENSG[0-9]+", gene_b))]
calls[, r_adjacent := (chr_a == chr_b & abs(pos_a - pos_b) < READTHROUGH_BP) %in% TRUE]
calls[, artifact   := r_name | r_ensg | r_adjacent | r_tag]

# ---- 8. Counts per tool and step ----------------------------------------------
count_row <- function(dt, step) {
  data.table(step = step, calls = nrow(dt),
             unique_fusions = uniqueN(dt$fusion), samples_with_calls = uniqueN(dt$sample))
}
summ <- rbindlist(lapply(names(res_list), function(tl) {
  d <- calls[tool == tl]
  rbindlist(list(
    count_row(d,                                "1_all_calls"),
    count_row(d[pass_tool == TRUE],             "2_tool_threshold"),
    count_row(d[pass_tool == TRUE & !artifact], "3_after_artifact_removal")
  ))[, tool := tl]
}))
setcolorder(summ, c("tool", "step", "calls", "unique_fusions", "samples_with_calls"))
fwrite(summ, file.path(OUT, "fusion_counts_per_tool.csv"))
print(summ)

# how many threshold-passing calls each artifact rule flags (rules overlap, so they do not add up)
rule_hits <- rbindlist(lapply(names(res_list), function(tl) {
  d <- calls[tool == tl & pass_tool == TRUE]
  data.table(tool = tl, calls_passing_threshold = nrow(d),
             gene_name_pattern = sum(d$r_name), no_gene_symbol = sum(d$r_ensg),
             adjacent_same_chr = sum(d$r_adjacent), database_tag = sum(d$r_tag),
             any_rule = sum(d$artifact))
}))
fwrite(rule_hits, file.path(OUT, "artifact_rule_hits.csv"))
print(rule_hits)

# ---- 9. Counts per sample (NA = no output file for that tool and sample) ------
grid <- CJ(sample = expected, tool = names(res_list))
grid[, has_file := mapply(function(s, tl) s %in% res_list[[tl]]$samples, sample, tool)]
n_thr <- calls[pass_tool == TRUE,               .(calls_tool_threshold = .N),       by = .(tool, sample)]
n_fin <- calls[pass_tool == TRUE & !artifact,   .(calls_after_artifacts = .N),      by = .(tool, sample)]
per_sample <- merge(merge(grid, n_thr, by = c("tool", "sample"), all.x = TRUE),
                    n_fin, by = c("tool", "sample"), all.x = TRUE)
per_sample[has_file == TRUE & is.na(calls_tool_threshold), calls_tool_threshold := 0L]
per_sample[has_file == TRUE & is.na(calls_after_artifacts), calls_after_artifacts := 0L]
per_sample[has_file == FALSE, c("calls_tool_threshold", "calls_after_artifacts") := NA_integer_]
setorder(per_sample, tool, sample)
fwrite(per_sample[, .(tool, sample, calls_tool_threshold, calls_after_artifacts)],
       file.path(OUT, "fusion_counts_per_sample.csv"))

# ---- 10. Tables to check the filters against the real labels / tags ----------
fc_labels <- calls[tool == "FusionCatcher", unlist(strsplit(tolower(labels_raw), ",", fixed = TRUE))]
fwrite(as.data.table(sort(table(trimws(fc_labels)), decreasing = TRUE))[, .(label = V1, calls = N)],
       file.path(OUT, "fc_label_tally.csv"))

sf_tokens <- calls[tool == "STAR-Fusion", unlist(strsplit(annots_raw, "[^A-Za-z0-9_]+"))]
sf_tokens <- sf_tokens[nchar(sf_tokens) >= 4 & !grepl("^[0-9]+$|^chr|Mb$", sf_tokens)]
fwrite(as.data.table(sort(table(sf_tokens), decreasing = TRUE))[, .(tag = sf_tokens, calls = N)],
       file.path(OUT, "starfusion_annots_tally.csv"))

# ---- 11. Filtered calls, one TSV per tool (same file names and columns sample/caller/fusion as your
#          earlier tables, so RNA_Fusion.Rmd only needs REPORT_ROOT changed to OUT) ----------------
keep <- calls[pass_tool == TRUE & !artifact,
              .(sample, caller = tool, fusion, gene_a, gene_b, chr_a, pos_a, chr_b, pos_b)]
out_names <- c(Arriba = "arriba.calls.high_medium.tsv",
               FusionCatcher = "fusioncatcher.calls.kept.tsv",
               `STAR-Fusion` = "starfusion.calls.filtered.tsv")
for (tl in names(out_names)) {
  fwrite(keep[caller == tl], file.path(OUT, out_names[[tl]]), sep = "\t")
}

# ---- 12. Figure (PNG + PDF) --------------------------------------------------
plot_dt <- melt(summ, id.vars = c("tool", "step"), measure.vars = c("calls", "unique_fusions"),
                variable.name = "metric", value.name = "n")
plot_dt[, step := factor(step, levels = c("1_all_calls", "2_tool_threshold", "3_after_artifact_removal"),
                         labels = c("All calls", "Tool threshold", "After artifact removal"))]
plot_dt[, metric := factor(metric, levels = c("calls", "unique_fusions"),
                           labels = c("Total calls", "Unique fusions"))]
# criteria shown under the plot; built from the settings above so the text always matches the filters used
wrap_txt <- function(x) paste(strwrap(x, width = 125), collapse = "\n")
crit_thr <- wrap_txt(paste0(
  "Tool threshold: Arriba = confidence ", paste(ARRIBA_CONF, collapse = " or "),
  "; FusionCatcher = calls without any FusionCatcher false-positive label (healthy-sample databases, ",
  "paralogs / similar reads, read-through / adjacent genes, pseudogene / ribosomal / mitochondrial genes, ",
  "overlapping genes); STAR-Fusion = FFPM >= ", STARF_FFPM, "."))
crit_art <- wrap_txt(paste0(
  "After artifact removal: tool-threshold calls without (i) a mitochondrial, small-RNA, ",
  "lncRNA or clone-style (e.g. RP11-) partner, (ii) a partner without a gene symbol (Ensembl ID only), ",
  "(iii) both breakpoints on the same chromosome and < ", format(READTHROUGH_BP / 1000), " kb apart, ",
  "(iv) STAR-Fusion normal-tissue / neighbour tags (", paste(STARF_DROP_ANNOTS, collapse = ", "), ")."))

p <- ggplot(plot_dt, aes(x = tool, y = n, fill = step)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 3) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("#B0B0B0", "#4C9F70", "#2C6E9B"), name = NULL) +
  labs(x = NULL, y = "Number of fusions",
       subtitle = paste0(length(expected), " samples. Unique fusion = gene pair, regardless of orientation."),
       caption = paste(crit_thr, crit_art, sep = "\n\n")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 7.5, lineheight = 1.1),
        plot.caption.position = "plot")
ggsave(file.path(OUT, "fusion_counts_per_tool.png"), p, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(OUT, "fusion_counts_per_tool.pdf"), p, width = 10, height = 6.5)

# ---- 13. Concordance across the three callers (uses the filtered calls `keep` from section 11) ----------
# Two levels:
#   cohort     : gene pair pooled over all samples (same logic as your earlier Venn diagram)
#   per sample : sample + gene pair must be called by the callers in the SAME sample (stricter)
# Agreement is at gene-pair level (A--B == B--A); breakpoints are listed in the table but not compared.
CALLERS   <- c("Arriba", "FusionCatcher", "STAR-Fusion")
CAT_ORDER <- c("Arriba+FusionCatcher+STAR-Fusion", "Arriba+FusionCatcher", "Arriba+STAR-Fusion",
               "FusionCatcher+STAR-Fusion", "Arriba", "FusionCatcher", "STAR-Fusion")
LV_COHORT <- "cohort (gene pair pooled over samples)"
LV_SAMPLE <- "per sample (sample + gene pair)"

# breakpoint text per call; one row per sample + fusion + caller (a caller may report a fusion more than once)
keep[, bp := sprintf("%s:%.0f-%s:%.0f", chr_a, pos_a, chr_b, pos_b)]
per <- keep[, .(bp = paste(unique(bp), collapse = ";")), by = .(sample, fusion, caller)]

# one row per sample + fusion, one column per caller (NA = not called by that caller)
wide <- dcast(per, sample + fusion ~ caller, value.var = "bp")
for (cl in CALLERS) if (!cl %in% names(wide)) wide[, (cl) := NA_character_]
wide[, n_callers := rowSums(!is.na(.SD)), .SDcols = CALLERS]
wide[, callers := apply(!is.na(.SD), 1, function(x) paste(CALLERS[x], collapse = "+")), .SDcols = CALLERS]
wide[, c("gene1", "gene2") := tstrsplit(fusion, "--", fixed = TRUE)]   # alphabetical order, NOT 5' / 3'
setnames(wide, CALLERS, paste0("breakpoints_", gsub("-", "_", CALLERS)))

# 5' > 3' order as reported by each caller (Arriba gene1 > gene2, FusionCatcher Gene_1 > Gene_2,
# STAR-Fusion LeftGene > RightGene). The breakpoints_* text above is in the same order (5' partner first).
keep[, orient := paste0(gene_a, ">", gene_b)]
per_o  <- keep[, .(orient = paste(unique(orient), collapse = ";")), by = .(sample, fusion, caller)]
wide_o <- dcast(per_o, sample + fusion ~ caller, value.var = "orient")
for (cl in CALLERS) if (!cl %in% names(wide_o)) wide_o[, (cl) := NA_character_]
O_COLS <- paste0("five_to_three_", gsub("-", "_", CALLERS))
setnames(wide_o, CALLERS, O_COLS)
wide <- merge(wide, wide_o[, c("sample", "fusion", O_COLS), with = FALSE], by = c("sample", "fusion"), all.x = TRUE)
# TRUE when every call of this gene pair in this sample has the same 5' > 3' order
wide <- merge(wide, keep[, .(orientation_consistent = uniqueN(orient) == 1L), by = .(sample, fusion)],
              by = c("sample", "fusion"), all.x = TRUE)
# One 5' > 3' order per fusion: the order supported by the MOST callers (for these data: all 3 callers where possible).
# The other order, if any, stays visible in orientation_other. orientation_tie = TRUE when two orders have the same
# number of callers (the first one listed is then arbitrary: check the breakpoints before using it).
ori_support <- keep[, .(n_callers_ori = uniqueN(caller), callers_ori = paste(sort(unique(caller)), collapse = "+")),
                    by = .(sample, fusion, orient)]
setorder(ori_support, sample, fusion, -n_callers_ori, orient)
ori_final <- ori_support[, .(
  orientation_final         = orient[1],
  orientation_final_callers = callers_ori[1],
  orientation_other         = if (.N > 1) paste0(orient[-1], " (", callers_ori[-1], ")", collapse = " ; ") else NA_character_,
  orientation_tie           = .N > 1 && n_callers_ori[1] == n_callers_ori[2]
), by = .(sample, fusion)]
wide <- merge(wide, ori_final, by = c("sample", "fusion"), all.x = TRUE)
setcolorder(wide, c("sample", "fusion", "gene1", "gene2", "n_callers", "callers", "orientation_consistent",
                    "orientation_final", "orientation_final_callers", "orientation_other", "orientation_tie", O_COLS))
setorder(wide, -n_callers, fusion, sample)
fwrite(wide, file.path(OUT, "fusion_concordance_per_sample.tsv"), sep = "\t")

# --- Orientation details: one row per sample + gene pair + order (5'>3') -------------------------
# Shows which callers support each order, how many calls, and the breakpoints, so the decision
# (FINAL vs alternative) can be checked. Rule: most callers first, then most calls; ties are flagged.
ori_det <- keep[, .(callers     = paste(sort(unique(caller)), collapse = "+"),
                    n_callers   = uniqueN(caller),
                    n_calls     = .N,
                    breakpoints = paste(unique(bp), collapse = " ; ")),
                by = .(sample, fusion, orient)]
setorder(ori_det, sample, fusion, -n_callers, -n_calls)
ori_det[, decision := fifelse(seq_len(.N) == 1, "FINAL", "alternative"), by = .(sample, fusion)]
ori_det[, tie := .N > 1 && n_callers[1] == n_callers[2] && n_calls[1] == n_calls[2], by = .(sample, fusion)]
fwrite(ori_det, file.path(OUT, "fusion_orientation_details.tsv"), sep = "\t")

# --- Recurrent genes among the fusions found by all 3 callers (one row per gene + fusion) -------
# A gene is recurrent if it is part of > 1 fusion. role = position of that gene in the final
# 5'>3' order; breakpoints are per caller, 5' partner first (chr:pos-chr:pos).
tri <- wide[n_callers == 3]
tri[, c("g5", "g3") := tstrsplit(orientation_final, ">", fixed = TRUE)]
bp_cols <- grep("^breakpoints_", names(tri), value = TRUE)
mk <- function(gcol, partner_col, role) {
  tri[, c(list(gene = get(gcol), role = role, partner = get(partner_col),
              sample = sample, fusion = fusion, orientation_final = orientation_final),
          lapply(.SD, identity)), .SDcols = bp_cols]
}
rec_long <- rbind(mk("g5", "g3", "5prime"), mk("g3", "g5", "3prime"))
rec_long[, `:=`(n_fusions = uniqueN(fusion), n_samples = uniqueN(sample)), by = gene]
rec <- rec_long[n_fusions > 1]
setcolorder(rec, c("gene", "n_fusions", "n_samples", "role", "partner", "sample", "fusion", "orientation_final", bp_cols))
setorder(rec, -n_fusions, gene, sample, fusion)
fwrite(rec, file.path(OUT, "recurrent_gene_breakpoints.tsv"), sep = "\t")

# cohort level: which callers reported each gene pair anywhere in the cohort
coh <- per[, .(callers = paste(CALLERS[CALLERS %in% unique(caller)], collapse = "+")), by = fusion]

# summary counts for every caller combination (zeros included)
counts <- rbind(coh[,  .(n = .N), by = callers][, level := LV_COHORT],
                wide[, .(n = .N), by = callers][, level := LV_SAMPLE])
summ_conc <- merge(CJ(level = c(LV_COHORT, LV_SAMPLE), callers = CAT_ORDER),
                   counts, by = c("level", "callers"), all.x = TRUE)
summ_conc[is.na(n), n := 0L]
summ_conc[, pct_of_level := round(100 * n / sum(n), 1), by = level]
summ_conc[, callers := factor(callers, levels = CAT_ORDER)]
setorder(summ_conc, level, callers)
fwrite(summ_conc, file.path(OUT, "fusion_concordance_summary.tsv"), sep = "\t")
print(summ_conc)

# list of gene pairs called by all 3 callers (cohort level), same file name / column as your earlier table
all3 <- coh[callers == CAT_ORDER[1], sort(fusion)]
fwrite(data.table(Fusion = all3), file.path(OUT, "fusions_in_all_3_callers.tsv"), sep = "\t")

message("Called by all 3 callers: ", length(all3), " gene pairs (cohort level), ",
        sum(wide$n_callers == 3), " sample-specific calls (per sample).")

# ---- 14. Venn diagrams (same style as RNA_Fusion.Rmd; PNG + PDF) ---------------------------------------
draw_venn <- function(sets, title, file_stub) {
  v <- venn.diagram(x = sets, filename = NULL, main = title, main.cex = 1.4,
                    fill = c("#1b9e77", "#F8766D", "#00BFC4"), alpha = 0.5,
                    cex = 1.1, cat.cex = 1.1, cat.col = "black", margin = 0.08,
                    euler.d = FALSE, scaled = FALSE,   # standard 3-circle layout, even when a region is empty
                    disable.logging = TRUE)
  png(file.path(OUT, paste0(file_stub, ".png")), width = 1800, height = 1400, res = 200)
  grid.newpage(); grid.draw(v); invisible(dev.off())
  pdf(file.path(OUT, paste0(file_stub, ".pdf")), width = 8, height = 6)
  grid.newpage(); grid.draw(v); invisible(dev.off())
}

# cohort level: gene pairs pooled over samples (same as the Venn on your slide)
draw_venn(lapply(setNames(CALLERS, CALLERS), function(cl) per[caller == cl, unique(fusion)]),
          "Fusion overlap across callers (cohort-level)", "fusion_venn_cohort")

message("Done. Tables and figures saved in: ", OUT)