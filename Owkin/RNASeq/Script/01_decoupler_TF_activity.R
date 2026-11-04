# ==============================================================================
# Script 1: ESR1 activity with decoupleR, 8 Primary vs 8 Recurrence (paired)
# Output: scores, paired stats, paired plots (PDF + PNG)
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(decoupleR)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(readr)
})
select <- dplyr::select   # AnnotationDbi masks dplyr::select

SEED <- 42
set.seed(SEED)

# ---- Settings (CHANGE THESE) --------------------------------------------------
BASE        <- "/ix1/alee/LO_LAB/Personal/Rahul/Okwin/RNASeq"
COUNTS_FILE <- file.path(BASE, "counts/matrix/gene_counts.tsv")
META_FILE   <- file.path(BASE, "sample_metadata.csv")   # columns: block, patient, sample_type
OUT_DIR     <- file.path(BASE, "Decoupler")   # script 1 (decoupler) results only
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# CollecTRI already downloaded/converted in the pILC project (reused, no new download)
NET_DIR  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Decoupler"
NET_GZ   <- file.path(NET_DIR, "collectri-human.tsv.gz")   # downloaded table
NET_FILE <- file.path(NET_DIR, "collectri_human.tsv")      # converted network (source, target, mor)

# Save a ggplot as both PDF and PNG
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# ---- 1. Metadata: keep complete Primary/Recurrence pairs ----------------------
meta <- read.csv(META_FILE, stringsAsFactors = FALSE) %>%
  mutate(
    sample    = paste0(block, "_mRNA"),
    timepoint = factor(ifelse(sample_type == "Primary", "Primary", "Recurrence"),
                       levels = c("Primary", "Recurrence"))
  )

counts <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE)

meta <- meta %>%
  filter(sample %in% colnames(counts)) %>%
  group_by(patient) %>%
  filter(any(timepoint == "Primary"), any(timepoint == "Recurrence")) %>%
  ungroup() %>%
  mutate(patient = factor(patient))

stopifnot(all(table(meta$patient) == 2))   # exactly one Primary + one Recurrence per patient
cat("Pairs used:", nlevels(meta$patient), "\n")

counts <- counts[, meta$sample]

# ---- 2. Gene IDs -> symbols (keep highest-expressed row per symbol) -----------
ens    <- sub("\\..*$", "", rownames(counts))
symbol <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep_g <- !is.na(symbol)
counts <- counts[keep_g, ]
symbol <- unname(symbol[keep_g])
ord    <- order(rowSums(counts), decreasing = TRUE)
counts <- counts[ord, ]
symbol <- symbol[ord]
dedup  <- !duplicated(symbol)
counts <- counts[dedup, ]
rownames(counts) <- symbol[dedup]

# ---- 3. Filter, TMM normalize, log-CPM ----------------------------------------
y <- DGEList(counts = counts, group = meta$timepoint)
y <- y[filterByExpr(y, group = meta$timepoint), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
logcpm <- cpm(y, log = TRUE, prior.count = 2)
colnames(logcpm) <- meta$sample

# ---- 4. TF activity: CollecTRI + ulm ------------------------------------------
# 4a. Network: reuse the converted pILC file; otherwise convert the downloaded table
net_path <- if (file.exists(NET_FILE)) NET_FILE else file.path(OUT_DIR, "collectri_human.tsv")
if (!file.exists(net_path)) {
  if (!file.exists(NET_GZ)) stop("Downloaded network not found: ", NET_GZ)
  raw <- read_tsv(NET_GZ, show_col_types = FALSE, col_types = cols(.default = "c"))
  stopifnot(all(c("source_genesymbol", "target_genesymbol",
                  "is_stimulation", "is_inhibition") %in% colnames(raw)))
  is_true <- function(x) x %in% c("True", "TRUE", "true")
  net0 <- raw %>%
    transmute(source = source_genesymbol, target = target_genesymbol,
              stim = is_true(is_stimulation), inhib = is_true(is_inhibition)) %>%
    mutate(mor = case_when(stim & !inhib ~ 1, inhib & !stim ~ -1, TRUE ~ NA_real_))
  cat("Unsigned/ambiguous interactions removed:", sum(is.na(net0$mor)), "\n")
  net0 <- net0 %>% filter(!is.na(mor), !is.na(source), !is.na(target))
  # Same TF-target pair with conflicting signs is ambiguous: remove it
  conflict <- net0 %>% group_by(source, target) %>%
    filter(n_distinct(mor) > 1) %>% ungroup() %>% distinct(source, target)
  cat("Pairs with conflicting signs removed:", nrow(conflict), "\n")
  net0 <- net0 %>% anti_join(conflict, by = c("source", "target")) %>% distinct(source, target, mor)
  write_tsv(net0, net_path)
}
net <- read_tsv(net_path, show_col_types = FALSE)
stopifnot(all(c("source", "target", "mor") %in% colnames(net)), "ESR1" %in% net$source)
cat("Network:", nrow(net), "interactions |", length(unique(net$source)), "TFs\n")

# 4b. z-score each gene across the 16 samples, so scores reflect differences
#     between samples and not baseline gene expression
keep   <- apply(logcpm, 1, sd) > 0
expr_z <- t(scale(t(logcpm[keep, ])))
cat("Network targets found in expression matrix:", length(intersect(net$target, rownames(expr_z))), "\n")

acts_tf <- run_ulm(mat = expr_z, network = net, .source = source, .target = target,
                   .mor = mor, minsize = 5) %>%
  filter(statistic == "ulm", source %in% c("ESR1", "FOXA1", "GATA3")) %>%
  select(feature = source, sample = condition, value = score)

# ---- 5. Hallmark estrogen response (cross-check) ------------------------------
msig <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
net_er <- msig %>%
  filter(gs_name %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE")) %>%
  transmute(source = gs_name, target = gene_symbol, mor = 1)

acts_er <- run_ulm(mat = expr_z, network = net_er, .source = source, .target = target,
                   .mor = mor, minsize = 5) %>%
  filter(statistic == "ulm") %>%
  select(feature = source, sample = condition, value = score)

acts <- bind_rows(acts_tf, acts_er) %>% left_join(meta, by = "sample")
write.csv(acts, file.path(OUT_DIR, "decoupler_activity_scores.csv"), row.names = FALSE)

# ---- 6. Paired statistics (Recurrence - Primary) ------------------------------
paired_test <- function(df) {
  w <- df %>%
    select(patient, timepoint, value) %>%
    pivot_wider(names_from = timepoint, values_from = value) %>%
    drop_na()
  tibble(
    n_pairs      = nrow(w),
    median_delta = median(w$Recurrence - w$Primary),
    mean_delta   = mean(w$Recurrence - w$Primary),
    n_up         = sum(w$Recurrence > w$Primary),
    p_wilcoxon   = suppressWarnings(wilcox.test(w$Recurrence, w$Primary, paired = TRUE, exact = TRUE)$p.value),
    p_ttest      = t.test(w$Recurrence, w$Primary, paired = TRUE)$p.value
  )
}

act_stats <- acts %>%
  group_by(feature) %>% group_modify(~ paired_test(.x)) %>% ungroup() %>%
  mutate(
    role = ifelse(feature == "ESR1", "primary", "secondary"),
    FDR_secondary = ifelse(role == "secondary", p.adjust(p_wilcoxon, "BH"), NA_real_)
  )
# ESR1 is the single primary hypothesis (no adjustment); the other features are secondary.
write.csv(act_stats, file.path(OUT_DIR, "decoupler_activity_paired_stats.csv"), row.names = FALSE)
print(act_stats)

# ---- 7. Figures ---------------------------------------------------------------
p_pair <- ggplot(acts, aes(timepoint, value, group = patient, colour = patient)) +
  geom_line(alpha = 0.8) + geom_point(size = 2) +
  facet_wrap(~feature, scales = "free_y",
             labeller = as_labeller(c(
               ESR1  = "ESR1", FOXA1 = "FOXA1", GATA3 = "GATA3",
               HALLMARK_ESTROGEN_RESPONSE_EARLY = "Hallmark estrogen\nresponse (early)",
               HALLMARK_ESTROGEN_RESPONSE_LATE  = "Hallmark estrogen\nresponse (late)"))) +
  labs(x = NULL, y = "Activity score (ulm)", colour = "Patient") +
  theme_bw(base_size = 12)
save_fig(p_pair, "decoupler_activity_paired", w = 9, h = 6)

p_esr1 <- acts %>% filter(feature == "ESR1") %>%
  ggplot(aes(timepoint, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = patient), colour = "grey60") +
  geom_point(aes(colour = patient), size = 3) +
  labs(x = NULL, y = "ESR1 activity (ulm)", colour = "Patient",
       title = sprintf("ESR1 activity, Wilcoxon p = %.3f",
                       act_stats$p_wilcoxon[act_stats$feature == "ESR1"])) +
  theme_bw(base_size = 12)
save_fig(p_esr1, "decoupler_ESR1_paired", w = 5, h = 5)

sessionInfo()