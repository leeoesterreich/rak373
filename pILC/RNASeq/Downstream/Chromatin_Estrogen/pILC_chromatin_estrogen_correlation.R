############################################################
## PURPOSE
##   Is there a correlation between estrogen-related and chromatin-related
##   pathway activity in pILC RNA-seq? Uses the saved VST matrix.
##
## DESIGN (why each step exists)
##   1. Saved VST matrix     -> made before clustering, so labels cannot shape the scores.
##   2. Database gene sets   -> no hand-typed gene lists; every selection
##                              decision is saved in an audit table.
##   3. Shared genes removed -> a pair that shares genes is correlated by
##                              construction (sensitivity run keeps them).
##   4. One summary score per family (PC1) -> ONE primary test, not 400.
##   5. Three models: unadjusted / cluster-adjusted / cluster + cellularity
##      + TILs adjusted. Permutation P-values (inside clusters when the
##      cluster is adjusted for).
##   6. Random-family null   -> is the observed rho larger than the
##                              background correlation of unrelated genes?
##   7. Per-cluster rho      -> descriptive consistency check (n ~ 10 each).
##   8. Pathway-by-pathway matrix -> exploratory only, BH FDR.
##
## NOTE
##   Clusters were derived from the same expression data, so the
##   cluster-adjusted result is the meaningful one. Pathway scores from RNA
##   show coordinated expression, not a causal or mechanistic link.
############################################################

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ----------------------------
# Paths (fixed names, no date or time stamps)
# ----------------------------
input_dir  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/4_Count_file_HTSeq_stranded/VST_data_GENCODEv47"
vst_file   <- file.path(input_dir, "htseq_VST_matrix_all_expression_filtered_genes.rds")   # saved VST, all expression-filtered genes
annot_file <- file.path(input_dir, "gene_annotation_from_GTF.tsv")

cluster_gct_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/NMF/GenePattern_NMFConsensus/vst_transformed_data.consensus.k.3.gct"
c5_gmt_file      <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Pathway_analysis/Dataset/c5.all.v2025.1.Hs.symbols.gmt"
h_gmt_file       <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Pathway_analysis/Dataset/h.all.v2025.1.Hs.symbols.gmt"
clinical_file    <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/Clinical_correlation/pILC_ClinicoPathological_details.csv"

out_dir   <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Downstream/Chromatin_Estrogen_correlation"
table_dir <- file.path(out_dir, "tables")
plot_dir  <- file.path(out_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)

# Cluster membership: (0) Cluster_membership_k3_GPlabels.csv, (1) the NMF GCT at the default path, (2) any k=3 consensus GCT found under
# the NMF folder, (3) cluster_metadata.tsv saved by the earlier pathway-score run.
nmf_root         <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq"
meta_fallback    <- file.path(nmf_root, "Downstream/Chromatin_Estrogen_pathway_scores/tables/cluster_metadata.tsv")
cluster_csv_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3_GPlabels.csv"
cluster_source   <- NULL
if (file.exists(cluster_csv_file)) {
  cluster_source <- "csv"
} else if (file.exists(cluster_gct_file)) {
  cluster_source <- "gct"
} else {
  hits <- list.files(nmf_root, pattern = "consensus\\.k\\.3\\.gct$", recursive = TRUE, full.names = TRUE)
  if (length(hits) >= 1) {
    cluster_gct_file <- hits[1]; cluster_source <- "gct"
    message("Default cluster GCT not found; using: ", cluster_gct_file)
    if (length(hits) > 1) message("Other k=3 GCT files found: ", paste(hits[-1], collapse = " | "))
  } else if (file.exists(meta_fallback)) {
    cluster_source <- "metadata"
    message("No k=3 consensus GCT found; using the saved cluster_metadata.tsv: ", meta_fallback)
  } else {
    stop("No cluster file found. Set cluster_gct_file to your k=3 NMF consensus GCT.")
  }
}

for (f in c(vst_file, annot_file, c5_gmt_file, h_gmt_file, clinical_file)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

# ----------------------------
# Settings
# ----------------------------
min_set_size   <- 15      # GSEA default minimum
max_set_size   <- 500     # GSEA default maximum
jaccard_cutoff <- 0.5     # near-duplicate sets inside a family are collapsed
n_perm         <- 10000   # permutations for the primary test
n_random       <- 1000    # random family pairs for the null
top_n_bubble   <- 10      # pathways shown in each top-pathway bubble plot
n_max_leading_genes <- 60  # genes drawn in each leading-edge heatmap
top_rank_by    <- "correlation"   # "correlation" = |rho| with the other family's score (adjusted); "cluster" = Kruskal-Wallis across clusters
cluster_colors <- c("Cluster 1" = "#1B9E77", "Cluster 2" = "#D95F02", "Cluster 3" = "#7570B3")

# ----------------------------
# Helper functions
# ----------------------------
clean_sample <- function(x) sub("_S[0-9]+$", "", trimws(as.character(x)))

clinical_key <- function(x) {   # join key only: TP03-117, TP03_117, TP03.117, XTP03.117_S5 and TP03-117-T all -> TP03117
  x <- toupper(trimws(as.character(x)))
  x <- sub("^X(?=TP)", "", x, perl = TRUE)
  x <- sub("_S[0-9]+$", "", x)
  x <- sub("^JC287-", "", x)
  x <- gsub("[._]", "-", x)
  x <- sub("-T$", "", x)
  x <- sub("(?<=[0-9])T$", "", x, perl = TRUE)
  x <- sub("-0+(?=[0-9])", "-", x, perl = TRUE)      # TP09-09 = TP09-9
  gsub("[^A-Z0-9]", "", x)
}

parse_pct <- function(x) {      # "<5" is coded as 2.5 (half the bound)
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NaN", "nan")] <- NA_character_
  v  <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", x)))
  lt <- grepl("^<", x) & !is.na(v)
  v[lt] <- v[lt] / 2
  v
}

read_gmt <- function(file) {
  lines <- strsplit(readLines(file), "\t")
  setNames(lapply(lines, function(v) toupper(v[-(1:2)])), vapply(lines, `[`, character(1), 1))
}

save_plot <- function(p, stem, w, h) {
  ggsave(file.path(plot_dir, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(plot_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
}

# ----------------------------
# 1. Cluster membership
# ----------------------------
if (cluster_source == "csv") {
  cc <- fread(cluster_csv_file, data.table = TRUE)
  message("Cluster CSV columns: ", paste(names(cc), collapse = ", "))
  samp_col <- intersect(names(cc), c("Sample", "sample", "Name", "name", "Sample_ID", "SampleID", "ID"))
  samp_col <- if (length(samp_col)) samp_col[1] else names(cc)[1]
  clus_col <- grep("cluster|membership|group", names(cc), ignore.case = TRUE, value = TRUE)
  clus_col <- setdiff(clus_col, samp_col)
  if (!length(clus_col)) clus_col <- setdiff(names(cc), samp_col)[1]
  clus_col <- clus_col[1]
  message("Using sample column '", samp_col, "' and cluster column '", clus_col, "'")
  k_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(cc[[clus_col]]))))
  meta <- data.table(Sample = clean_sample(cc[[samp_col]]),
                     Cluster = factor(paste0("Cluster ", k_num), levels = names(cluster_colors)))
} else if (cluster_source == "gct") {
  cl <- fread(cluster_gct_file, skip = 2, data.table = TRUE)
  membership_col <- grep("^membership", names(cl), value = TRUE)
  if (length(membership_col) != 1) stop("Expected exactly one membership column.")
  meta <- cl[, .(Sample = clean_sample(Name),
                 Cluster = factor(paste0("Cluster ", as.integer(get(membership_col))),
                                  levels = names(cluster_colors)))]
} else {
  mf <- fread(meta_fallback)
  meta <- mf[, .(Sample = clean_sample(Sample),
                 Cluster = factor(as.character(Cluster), levels = names(cluster_colors)))]
}
if (anyNA(meta$Cluster)) stop("Cluster membership must be 1, 2 or 3.")
if (anyDuplicated(meta$Sample)) stop("Duplicate samples in cluster file.")
setorder(meta, Cluster)
n <- nrow(meta)
print(table(meta$Cluster))
fwrite(meta, file.path(table_dir, "cluster_metadata.tsv"), sep = "\t")

# ----------------------------
# 2. Covariates: tumor cellularity and TILs (purity / immune infiltration)
# ----------------------------
clin <- as.data.table(read.csv(clinical_file, check.names = FALSE, strip.white = TRUE,
                               na.strings = c("", "NA", "NaN", "N/A"), stringsAsFactors = FALSE,
                               fileEncoding = "UTF-8-BOM"))
setnames(clin, gsub("^\\s+|\\s+$", "", names(clin)))   # stray spaces in column names
message("Clinical file: ", nrow(clin), " rows; columns: ", paste(names(clin), collapse = " | "))
need <- c("Sample", "tumor_cellularity (%)", "TILs (%)")
if (!all(need %in% names(clin))) stop("Clinical file is missing: ", paste(setdiff(need, names(clin)), collapse = ", "))
clin_dt <- unique(clin[, .(Key = clinical_key(Sample),
                           Cellularity = parse_pct(`tumor_cellularity (%)`),
                           TILs = parse_pct(`TILs (%)`))], by = "Key")

key_rna <- clinical_key(meta$Sample)
cov_df <- data.frame(Cluster = meta$Cluster,
                     Cellularity = clin_dt$Cellularity[match(key_rna, clin_dt$Key)],
                     TILs        = clin_dt$TILs[match(key_rna, clin_dt$Key)],
                     row.names = meta$Sample)
cov_missing <- data.table(Sample = meta$Sample,
                          Cellularity_missing = is.na(cov_df$Cellularity),
                          TILs_missing = is.na(cov_df$TILs))
fwrite(cov_missing, file.path(table_dir, "covariate_missingness.tsv"), sep = "\t")
message("Samples matched to the clinical file: ", sum(key_rna %in% clin_dt$Key), " of ", length(key_rna))
if (sum(key_rna %in% clin_dt$Key) < length(key_rna))
  message("Unmatched RNA samples: ", paste(meta$Sample[!key_rna %in% clin_dt$Key], collapse = ", "))
message("Missing cellularity: ", sum(is.na(cov_df$Cellularity)),
        " | missing TILs: ", sum(is.na(cov_df$TILs)), " (filled with the cohort median)")
for (v in c("Cellularity", "TILs")) cov_df[[v]][is.na(cov_df[[v]])] <- median(cov_df[[v]], na.rm = TRUE)

# ----------------------------
# 3. Saved VST matrix -> gene symbols
#    The VST was made before the clustering (the NMF clusters were built from
#    vst_transformed_data.gct in this folder), so the labels cannot have shaped it.
# ----------------------------
vst_all <- readRDS(vst_file)
vst_all <- as.matrix(vst_all)
if (is.null(rownames(vst_all)) || is.null(colnames(vst_all))) stop("VST matrix needs gene row names and sample column names.")
vst_key  <- clinical_key(colnames(vst_all))      # TP09.9, TP09-09, TP09_9_S3 ... all -> same key
meta_key <- clinical_key(meta$Sample)
if (anyDuplicated(vst_key))  stop("Two VST columns map to the same sample key.")
if (anyDuplicated(meta_key)) stop("Two cluster samples map to the same sample key.")
idx <- match(meta_key, vst_key)
if (anyNA(idx)) stop("Cluster samples missing in the VST matrix: ", paste(meta$Sample[is.na(idx)], collapse = ", "))
extra <- colnames(vst_all)[!vst_key %in% meta_key]
if (length(extra)) message("VST samples not in the cluster file (dropped): ", paste(extra, collapse = ", "))
vst_mat <- vst_all[, idx, drop = FALSE]
colnames(vst_mat) <- meta$Sample                 # cluster-file names are used from here on
storage.mode(vst_mat) <- "double"
if (anyNA(vst_mat)) stop("Missing values in the VST matrix.")
message("VST matrix: ", nrow(vst_mat), " genes x ", ncol(vst_mat), " samples")

# Row names may be symbols, Ensembl IDs, or a mix: map Ensembl IDs through the annotation file
annot <- fread(annot_file)
rn   <- rownames(vst_mat)
sym  <- toupper(rn)
is_e <- grepl("^ENSG[0-9]+", sym)
if (any(is_e)) {
  ann_key <- sub("\\.[0-9]+$", "", annot$Gene)
  sym[is_e] <- toupper(annot$gene_symbol[match(sub("\\.[0-9]+$", "", rn[is_e]), ann_key)])
}
sym_dt <- data.table(Row = rn, Symbol = sym, Mean_VST = rowMeans(vst_mat))
sym_dt <- sym_dt[!is.na(Symbol) & Symbol != "" & !grepl("^ENSG[0-9]+", Symbol)]
setorder(sym_dt, -Mean_VST)
sym_dt <- sym_dt[!duplicated(Symbol)]      # shared symbol: keep the most expressed gene
expr <- vst_mat[sym_dt$Row, , drop = FALSE]
rownames(expr) <- sym_dt$Symbol
expressed_genes <- rownames(expr)
message("Symbols available for scoring: ", length(expressed_genes))

# ----------------------------
# 4. Gene sets from databases (GO and Hallmark local GMT; Reactome, KEGG Medicus via msigdbr)
# ----------------------------
get_msig <- function(collection, subcollection = NULL) {
  res <- tryCatch(
    msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection),
    error = function(e) NULL)
  if (is.null(res)) {   # older msigdbr versions
    res <- tryCatch(
      msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection),
      error = function(e) NULL)
  }
  if (is.null(res) || nrow(res) == 0) {
    message("Collection not available in msigdbr: ", collection, " ", subcollection)
    return(list())
  }
  res <- as.data.table(res)
  lapply(split(res$gene_symbol, res$gs_name), function(g) unique(toupper(g)))
}

go_all  <- read_gmt(c5_gmt_file)
go_sets <- go_all[grepl("^(GOBP|GOCC|GOMF)_", names(go_all))]

hallmark_sets <- read_gmt(h_gmt_file)   # local Hallmark GMT (same release as C5)
db_sets <- c(go_sets, hallmark_sets, get_msig("C2", "CP:REACTOME"), get_msig("C2", "CP:KEGG_MEDICUS"))
db_sets <- db_sets[!duplicated(names(db_sets))]

set_source <- fifelse(grepl("^(GOBP|GOCC|GOMF)_", names(db_sets)), "GO",
              fifelse(grepl("^HALLMARK_", names(db_sets)), "Hallmark",
              fifelse(grepl("^REACTOME_", names(db_sets)), "Reactome",
              fifelse(grepl("^KEGG_MEDICUS", names(db_sets)), "KEGG_Medicus", "Other"))))
names(set_source) <- names(db_sets)

# ----------------------------
# 5. Rule-based family assignment on set NAMES (saved in an audit table)
# ----------------------------
estrogen_include <- "ESTROGEN|ESTRADIOL|ESR_MEDIATED"
estrogen_exclude <- paste(c("BIOSYNTH", "METABOL", "CATABOL", "SECRETION", "TRANSPORT",
                            "GLUCURON", "SULFO", "AROMATASE", "HYDROXYLASE", "ANTIESTROGEN"), collapse = "|")
chromatin_include <- paste(c("CHROMATIN", "HETEROCHROMATIN", "EUCHROMATIN", "NUCLEOSOME",
                             "HISTONE.*(ACETYL|DEACETYL|METHYL|DEMETHYL|UBIQUITIN|DEUBIQUITIN|PHOSPH|SUMO)",
                             "(ACETYL|DEACETYL|METHYL|DEMETHYL|UBIQUITIN|DEUBIQUITIN|PHOSPH|SUMO).*HISTONE",
                             "DNA_METHYLATION", "DNA_DEMETHYLATION", "SWI.*SNF", "BAF_COMPLEX",
                             "PBAF", "POLYCOMB", "PRC1", "PRC2", "NURD", "NURF", "INO80", "SWR1"), collapse = "|")
chromatin_exclude <- paste(c("(^|_)RNA_", "MRNA", "RRNA", "TRNA", "RIBOSOM", "METHYLADENOSINE",
                             "N6_METHYL", "PSEUDOURID", "HISTONE_MRNA", "ZYGOTIC", "GERM_CELL",
                             "SPERM", "OOCYTE", "ANTIBODY", "IMMUNOGLOBULIN"), collapse = "|")

audit <- data.table(
  Set = names(db_sets), Source = unname(set_source),
  Estrogen_include  = grepl(estrogen_include,  names(db_sets)),
  Estrogen_exclude  = grepl(estrogen_exclude,  names(db_sets)),
  Chromatin_include = grepl(chromatin_include, names(db_sets)),
  Chromatin_exclude = grepl(chromatin_exclude, names(db_sets))
)
audit[, Family := fifelse(Estrogen_include & !Estrogen_exclude, "Estrogen",
                  fifelse(Chromatin_include & !Chromatin_exclude & !Estrogen_include, "Chromatin",
                          NA_character_))]
audit[, Expressed_size := vapply(db_sets[Set], function(g) sum(g %in% expressed_genes), integer(1))]
audit[, Size_ok := Expressed_size >= min_set_size & Expressed_size <= max_set_size]

cand <- audit[!is.na(Family) & Size_ok]
cat("Candidate sets before redundancy collapse:\n")
print(cand[, .N, by = .(Family, Source)])

# Collapse near-duplicates inside each family. Priority: Hallmark > Reactome > KEGG Medicus > GO.
priority <- c(Hallmark = 1, Reactome = 2, KEGG_Medicus = 3, GO = 4, Other = 5)
cand[, Priority := priority[Source]]
setorder(cand, Family, Priority, -Expressed_size)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
expressed_sets <- lapply(db_sets[cand$Set], intersect, y = expressed_genes)

cand[, Kept := FALSE]
cand[, Collapsed_into := NA_character_]
for (fam in unique(cand$Family)) {
  kept <- character()
  for (s in cand[Family == fam, Set]) {
    hit <- NULL
    for (k in kept) {
      if (jaccard(expressed_sets[[s]], expressed_sets[[k]]) >= jaccard_cutoff) { hit <- k; break }
    }
    if (is.null(hit)) {
      kept <- c(kept, s)
      cand[Set == s, Kept := TRUE]
    } else {
      cand[Set == s, Collapsed_into := hit]
    }
  }
}
audit <- merge(audit, cand[, .(Set, Kept, Collapsed_into)], by = "Set", all.x = TRUE, sort = FALSE)
audit[, Final_decision := fifelse(is.na(Family), "Not chromatin/estrogen by name rule",
                          fifelse(!Size_ok, "Excluded by gene-set size",
                          fifelse(Kept %in% TRUE, "Retained", "Collapsed (redundant)")))]
fwrite(audit, file.path(table_dir, "geneset_selection_audit.tsv"), sep = "\t")

final <- cand[Kept == TRUE]
if (!all(c("Chromatin", "Estrogen") %in% final$Family)) stop("A family has no retained sets.")
fwrite(final[, .(Set, Family, Source, Expressed_size)], file.path(table_dir, "retained_gene_sets.tsv"), sep = "\t")
cat("Final sets retained:\n"); print(final[, .N, by = .(Family, Source)])

sets0   <- expressed_sets[final$Set]
fam     <- setNames(final$Family, final$Set)
est_ids <- names(fam)[fam == "Estrogen"]
chr_ids <- names(fam)[fam == "Chromatin"]

# Genes present in both families: removed in the primary analysis
est_genes    <- unique(unlist(sets0[est_ids]))
chr_genes    <- unique(unlist(sets0[chr_ids]))
shared_genes <- intersect(est_genes, chr_genes)
message("Estrogen genes: ", length(est_genes), " | chromatin genes: ", length(chr_genes),
        " | shared (removed in primary): ", length(shared_genes))
fwrite(data.table(Shared_gene = sort(shared_genes)), file.path(table_dir, "shared_genes_removed.tsv"), sep = "\t")
sets_clean <- lapply(sets0, setdiff, y = shared_genes)

label_map <- function(ids) {
  x <- sub("^(GOBP|GOCC|GOMF|HALLMARK|REACTOME|KEGG_MEDICUS_REFERENCE|KEGG_MEDICUS)_", "", ids)
  x <- tolower(gsub("_", " ", x))
  x <- paste0(toupper(substr(x, 1, 1)), substring(x, 2))
  make.unique(paste0(x, " [", set_source[ids], "]"), sep = " #")
}

# ----------------------------
# 6. ssGSEA
# ----------------------------
run_ssgsea <- function(expr, sets) {
  sets <- sets[lengths(sets) >= min_set_size & lengths(sets) <= max_set_size]
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    sc <- GSVA::gsva(GSVA::ssgseaParam(exprData = expr, geneSets = sets,
                                       minSize = min_set_size, maxSize = max_set_size, normalize = TRUE),
                     verbose = FALSE)
  } else {
    sc <- GSVA::gsva(expr = expr, gset.idx.list = sets, method = "ssgsea", kcdf = "Gaussian",
                     min.sz = min_set_size, max.sz = max_set_size, abs.ranking = TRUE,
                     ssgsea.norm = TRUE, verbose = FALSE)
  }
  as.matrix(sc)[, meta$Sample, drop = FALSE]
}

scores_clean <- run_ssgsea(expr, sets_clean)   # primary: shared genes removed
scores_full  <- run_ssgsea(expr, sets0)        # sensitivity: shared genes kept
est_c <- intersect(est_ids, rownames(scores_clean)); chr_c <- intersect(chr_ids, rownames(scores_clean))
est_f <- intersect(est_ids, rownames(scores_full));  chr_f <- intersect(chr_ids, rownames(scores_full))
message("Sets scored (primary): estrogen ", length(est_c), " | chromatin ", length(chr_c))
if (length(est_c) < 2 || length(chr_c) < 2) stop("Too few sets left after shared-gene removal.")

# Which pathway sets were actually used (also lists the genes in each set)
used_tbl <- data.table(
  Set = final$Set, Family = final$Family, Source = final$Source, Label = label_map(final$Set),
  Expressed_size = final$Expressed_size,
  Size_after_shared_gene_removal = lengths(sets_clean[final$Set]),
  Used_in_primary_analysis = final$Set %in% c(est_c, chr_c),
  Genes_after_shared_gene_removal = vapply(sets_clean[final$Set], paste, character(1), collapse = ";"))
setorder(used_tbl, Family, -Used_in_primary_analysis, Source, Set)
fwrite(used_tbl, file.path(table_dir, "gene_sets_used.tsv"), sep = "\t")
for (f in c("Estrogen", "Chromatin")) {
  cat("\n", f, "sets used in the primary analysis (", used_tbl[Family == f & Used_in_primary_analysis, .N], "):\n", sep = "")
  print(used_tbl[Family == f & Used_in_primary_analysis, .(Label, Expressed_size, Size_after_shared_gene_removal)], row.names = FALSE)
}

fwrite(data.table(Set = rownames(scores_clean), Family = fam[rownames(scores_clean)],
                  Label = label_map(rownames(scores_clean)), scores_clean),
       file.path(table_dir, "ssGSEA_scores_shared_genes_removed.tsv"), sep = "\t")

# ----------------------------
# 7. Family summary score = PC1 of the z-scored pathway scores
# ----------------------------
family_pc1 <- function(scores, ids) {
  z  <- t(scale(t(scores[ids, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  pc <- prcomp(t(z))
  s  <- pc$x[, 1]
  if (mean(cor(s, t(z))) < 0) s <- -s            # higher = more expressed
  list(score = s, var_explained = pc$sdev[1]^2 / sum(pc$sdev^2))
}

# ----------------------------
# 8. Adjustment and permutation test
# ----------------------------
adjust <- function(M, vars) {                    # M: samples x features
  M <- as.matrix(M)
  if (length(vars) == 0) return(M)
  r <- resid(lm(M ~ ., data = cov_df[, vars, drop = FALSE]))
  if (is.null(dim(r))) r <- matrix(r, ncol = 1)
  r
}

# Spearman rho on (adjusted) residuals; permutation inside clusters when the cluster is adjusted for.
primary_test <- function(x, y, vars, B = n_perm) {
  R  <- adjust(cbind(x, y), vars)
  rx <- rank(R[, 1]); ry <- rank(R[, 2])
  obs <- cor(rx, ry)
  idx <- if ("Cluster" %in% vars) split(seq_len(n), cov_df$Cluster) else list(seq_len(n))
  null <- replicate(B, {
    yp <- ry
    for (i in idx) yp[i] <- ry[i][sample.int(length(i))]
    cor(rx, yp)
  })
  list(rho = obs, p = (1 + sum(abs(null) >= abs(obs))) / (B + 1))
}

models <- list(Unadjusted = character(),
               Cluster_adjusted = "Cluster",
               Cluster_cellularity_TILs_adjusted = c("Cluster", "Cellularity", "TILs"))

run_primary <- function(scores, est_ids, chr_ids, label) {
  e <- family_pc1(scores, est_ids)
  c <- family_pc1(scores, chr_ids)
  res <- rbindlist(lapply(names(models), function(m) {
    set.seed(42)
    t <- primary_test(e$score, c$score, models[[m]])
    data.table(Gene_sets = label, Model = m, Spearman_rho = t$rho, Permutation_P = t$p,
               Estrogen_PC1_variance = e$var_explained, Chromatin_PC1_variance = c$var_explained,
               N_estrogen_sets = length(est_ids), N_chromatin_sets = length(chr_ids))
  }))
  list(res = res, e = e$score, c = c$score)
}

p_clean <- run_primary(scores_clean, est_c, chr_c, "Shared genes removed (primary)")
p_full  <- run_primary(scores_full,  est_f, chr_f, "Shared genes kept (sensitivity)")
primary_tbl <- rbind(p_clean$res, p_full$res)
print(primary_tbl)
fwrite(primary_tbl, file.path(table_dir, "primary_family_correlation.tsv"), sep = "\t")

# Per-cluster rho (descriptive) and Fisher-z combination
within <- rbindlist(lapply(levels(meta$Cluster), function(k) {
  i <- cov_df$Cluster == k
  data.table(Cluster = k, N = sum(i),
             Spearman_rho = cor(p_clean$e[i], p_clean$c[i], method = "spearman"),
             Nominal_P = suppressWarnings(cor.test(p_clean$e[i], p_clean$c[i], method = "spearman", exact = FALSE)$p.value))
}))
z <- atanh(within$Spearman_rho); w <- within$N - 3
zbar <- sum(w * z) / sum(w); se <- 1 / sqrt(sum(w))
within <- rbind(within, data.table(Cluster = "Combined (Fisher z)", N = sum(within$N),
                                   Spearman_rho = tanh(zbar), Nominal_P = 2 * pnorm(-abs(zbar / se))), fill = TRUE)
print(within)
fwrite(within, file.path(table_dir, "per_cluster_family_correlation.tsv"), sep = "\t")

sc_dt <- data.table(Sample = meta$Sample, Cluster = meta$Cluster,
                    Estrogen_PC1 = p_clean$e, Chromatin_PC1 = p_clean$c)
fwrite(sc_dt, file.path(table_dir, "family_PC1_scores.tsv"), sep = "\t")

mc <- primary_tbl[Gene_sets == "Shared genes removed (primary)"]
get_m <- function(m, col) mc[Model == m][[col]]
subtitle_txt <- paste0(
  "Unadjusted rho = ", round(get_m("Unadjusted", "Spearman_rho"), 2),
  " (P = ", signif(get_m("Unadjusted", "Permutation_P"), 2), ")\n",
  "Adjusted for cluster + cellularity + TILs: rho = ",
  round(get_m("Cluster_cellularity_TILs_adjusted", "Spearman_rho"), 2),
  " (P = ", signif(get_m("Cluster_cellularity_TILs_adjusted", "Permutation_P"), 2), ")")

p_sc <- ggplot(sc_dt, aes(Chromatin_PC1, Estrogen_PC1, color = Cluster)) +
  geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "grey30", linewidth = 0.7) +
  scale_color_manual(values = cluster_colors) +
  labs(title = "Estrogen vs chromatin pathway summary scores (shared genes removed)",
       subtitle = subtitle_txt, x = "Chromatin score (PC1)", y = "Estrogen score (PC1)",
       caption = paste(strwrap(paste0(
         "Shared genes removed: the ", length(shared_genes), " genes found in both the estrogen and chromatin gene sets ",
         "were deleted from every set before ssGSEA. Clusters were derived from the same expression data; colours are for display."),
         110), collapse = "\n")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, size = 8))
save_plot(p_sc, "family_PC1_scatter", 8, 6.5)

# Per-cluster scatter (descriptive, about 10 tumors per cluster): does the trend hold inside each cluster?
cl_rows <- within[Cluster %in% levels(meta$Cluster)]
fmt_p2  <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
lab_map <- setNames(paste0(cl_rows$Cluster, " (n = ", cl_rows$N, ", rho = ", sprintf("%.2f", cl_rows$Spearman_rho),
                           ", P = ", fmt_p2(cl_rows$Nominal_P), ")"),
                    cl_rows$Cluster)
sc_dt[, Panel := factor(lab_map[as.character(Cluster)], levels = unname(lab_map))]
p_cl <- ggplot(sc_dt, aes(Chromatin_PC1, Estrogen_PC1, color = Cluster)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ Panel, nrow = 1, scales = "free") +
  scale_color_manual(values = cluster_colors) +
  labs(title = "Estrogen vs chromatin summary scores within each RNA cluster",
       subtitle = "Spearman rho and nominal P per cluster (descriptive; small n; P not corrected for multiple testing).",
       x = "Chromatin score (PC1)", y = "Estrogen score (PC1)",
       caption = paste(strwrap(paste0(
         "Shared genes removed: the ", length(shared_genes), " genes found in both the estrogen and chromatin gene sets ",
         "were deleted from every set before ssGSEA, so the two scores contain no common genes. ",
         "Clusters were derived from the same expression data."), 150), collapse = "\n")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        legend.position = "none", strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8))
save_plot(p_cl, "family_PC1_scatter_by_cluster", 13, 4.8)

# ----------------------------
# 8b. Heatmaps of pathway scores per tumor: one for estrogen, one for chromatin
#     Rows: pathways used in the primary analysis (shared genes removed), z-scored across tumors.
#     Columns: tumors grouped by RNA cluster. Cellularity / TILs shown as bars (missing = cohort median).
# ----------------------------
save_ht <- function(ht, stem, w, h) {
  for (ext in c("pdf", "png")) {
    f <- file.path(plot_dir, paste0(stem, ".", ext))
    if (ext == "pdf") pdf(f, width = w, height = h) else png(f, width = w, height = h, units = "in", res = 300)
    draw(ht, merge_legends = TRUE, padding = unit(c(4, 8, 6, 4), "mm"))
    dev.off()
  }
}
col_order <- meta$Sample                      # already ordered by cluster
cl_vec    <- meta$Cluster
top_ann <- HeatmapAnnotation(
  Cluster = cl_vec,
  `Tumor cellularity (%)` = cov_df[col_order, "Cellularity"],
  `TILs (%)` = cov_df[col_order, "TILs"],
  col = list(Cluster = cluster_colors,
             `Tumor cellularity (%)` = colorRamp2(c(0, 50, 100), c("#F7FBFF", "#6BAED6", "#08306B")),
             `TILs (%)` = colorRamp2(c(0, 30, 100), c("#FFF5F0", "#FB6A4A", "#67000D"))),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9),
  simple_anno_size = unit(0.4, "cm"))

family_heatmap <- function(ids, family, stem) {
  mat <- t(scale(t(scores_clean[ids, col_order, drop = FALSE])))
  mat[!is.finite(mat)] <- 0
  mat[mat > 2.5] <- 2.5; mat[mat < -2.5] <- -2.5
  ht <- Heatmap(mat, name = "Pathway score\n(Z across tumors)",
                col = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
                cluster_rows = nrow(mat) > 2, cluster_columns = FALSE,
                column_split = cl_vec, top_annotation = top_ann,
                row_labels = label_map(rownames(mat)), row_names_gp = gpar(fontsize = 8),
                row_names_max_width = unit(11, "cm"),
                show_column_names = FALSE, border = TRUE,
                column_title = paste0(family, " pathway scores across pILC tumors (shared genes removed)"))
  save_ht(ht, stem, w = 14, h = max(5, 0.22 * nrow(mat) + 3.5))
  fwrite(data.table(Set = rownames(scores_clean[ids, , drop = FALSE]), Label = label_map(ids),
                    scores_clean[ids, col_order, drop = FALSE]),
         file.path(table_dir, paste0(stem, "_scores.tsv")), sep = "\t")
}
family_heatmap(est_c, "Estrogen",  "estrogen_pathway_scores_heatmap")
family_heatmap(chr_c, "Chromatin", "chromatin_pathway_scores_heatmap")

# ----------------------------
# 8c. Top pathways: ONE bubble plot (top chromatin + top estrogen) and leading-edge heatmaps
#     Bubble: colour = median pathway score per cluster (Z across tumors); size = -log10 FDR of the
#     3-cluster Kruskal-Wallis test (FDR inside the family).
#     "Top" = ranked by |rho| with the other family's summary score (cluster + cellularity + TILs
#     adjusted); set top_rank_by <- "cluster" to rank by the cluster test instead.
#     Leading edge: fgsea of each top pathway on a gene ranking = adjusted Spearman rho of every gene
#     with the OTHER family's summary score (Cluster + cellularity + TILs adjusted). Genes in the
#     leading edge are the ones that drive the link between the two families.
# ----------------------------
vars_adj <- models[["Cluster_cellularity_TILs_adjusted"]]
k_cov_b  <- ncol(model.matrix(~ ., data = cov_df[, vars_adj, drop = FALSE])) - 1
df_b     <- n - 2 - k_cov_b

rank_family <- function(ids, other_pc1, family, stem) {
  sc  <- scores_clean[ids, meta$Sample, drop = FALSE]
  Rr  <- apply(adjust(cbind(t(sc), .other = other_pc1), vars_adj), 2, rank)
  colnames(Rr) <- c(ids, ".other")
  rho <- cor(Rr)[ids, ".other"]
  pv  <- 2 * pt(-abs(rho * sqrt(df_b / (1 - rho^2))), df_b)
  kw  <- apply(sc, 1, function(v) kruskal.test(v ~ meta$Cluster)$p.value)
  tab <- data.table(Family = family, Set = ids, Label = label_map(ids), Rho_adjusted = rho,
                    P_correlation = pv, FDR_correlation = p.adjust(pv, "BH"),
                    P_cluster = kw, FDR_cluster = p.adjust(kw, "BH"))
  tab <- tab[if (top_rank_by == "cluster") order(FDR_cluster, P_cluster) else order(-abs(Rho_adjusted))]
  fwrite(tab, file.path(table_dir, paste0(stem, "_pathway_ranking.tsv")), sep = "\t")
  tab
}
rk_chr <- rank_family(chr_c, p_clean$e, "Chromatin", "chromatin")
rk_est <- rank_family(est_c, p_clean$c, "Estrogen",  "estrogen")
top_chr <- head(rk_chr, top_n_bubble); top_est <- head(rk_est, top_n_bubble)
top_all <- rbind(top_chr, top_est)

# ---- single bubble plot
sc_top <- scores_clean[top_all$Set, meta$Sample, drop = FALSE]
z_top  <- t(scale(t(sc_top))); z_top[!is.finite(z_top)] <- 0
zl <- melt(as.data.table(z_top, keep.rownames = "Set"), id.vars = "Set", variable.name = "Sample", value.name = "Z")
zl <- merge(zl, meta, by = "Sample")
bub <- merge(zl[, .(Median_Z = median(Z)), by = .(Set, Cluster)],
             top_all[, .(Set, Family, Label, FDR_cluster)], by = "Set")
bub[, Family := factor(Family, levels = c("Chromatin", "Estrogen"))]
bub[, minus_log10_FDR := pmin(-log10(pmax(FDR_cluster, 1e-10)), 2.5)]
bub[, Label := factor(Label, levels = rev(top_all$Label))]
fwrite(bub[order(Family, Label, Cluster)], file.path(table_dir, "top_pathways_bubble_data.tsv"), sep = "\t")

rank_txt <- if (top_rank_by == "cluster") {
  "ranked by Kruskal-Wallis FDR across clusters"
} else {
  "ranked by |correlation| with the other family's score"
}
p_bub <- ggplot(bub, aes(Cluster, Label, fill = Median_Z, size = minus_log10_FDR)) +
  geom_point(shape = 21, color = "grey20", stroke = 0.4) +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Median score\n(Z across tumors)") +
  scale_size_continuous(name = "-log10 FDR\n(3-cluster test)", limits = c(0, 2.5),
                        breaks = c(0.5, 1.3, 2), labels = c("0.5", "1.3 (FDR 0.05)", "2"), range = c(2, 9)) +
  labs(title = "Top chromatin and estrogen pathway scores across pILC RNA clusters",
       subtitle = paste0("Top ", top_n_bubble, " per family, ", rank_txt,
                         ". Colour: relative score in the cluster; size: Kruskal-Wallis FDR."),
       x = NULL, y = NULL,
       caption = "ssGSEA on VST, shared genes removed. FDR computed within each family. Clusters were derived from the same expression data.") +
  theme_bw(base_size = 11) +
  theme(plot.title.position = "plot", plot.caption.position = "plot",
        plot.title = element_text(face = "bold"), strip.text.y = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 9), panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8))
save_plot(p_bub, "top_chromatin_estrogen_pathways_bubble", w = 12, h = max(6, 0.3 * nrow(top_all) + 3))

# ---- gene-level ranking: adjusted Spearman rho of every gene with the other family's score
Rg <- apply(adjust(t(expr), vars_adj), 2, rank)
colnames(Rg) <- rownames(expr)
gene_rho <- function(other_pc1) {
  ro <- rank(adjust(cbind(o = other_pc1), vars_adj)[, 1])
  r  <- drop(cor(Rg, ro))
  r[is.finite(r)]
}

leading_edge_heatmap <- function(top, other_pc1, family, other_name, stem) {
  ranks <- sort(gene_rho(other_pc1), decreasing = TRUE)
  pw <- lapply(sets_clean[top$Set], intersect, y = names(ranks))
  set.seed(42)
  fg <- as.data.table(fgseaMultilevel(pathways = pw, stats = ranks, minSize = min_set_size,
                                      maxSize = max_set_size, eps = 0, nproc = 1))
  if (nrow(fg) == 0) { message("fgsea returned no result for ", family); return(invisible(NULL)) }
  fg[, Label := label_map(top$Set)[match(pathway, top$Set)]]
  fwrite(fg[, .(pathway, Label, size, NES, pval, padj, leading_edge_size = lengths(leadingEdge),
                leading_edge = vapply(leadingEdge, paste, character(1), collapse = ";"))],
         file.path(table_dir, paste0(stem, "_top_pathways_fgsea.tsv")), sep = "\t")
  le <- rbindlist(lapply(seq_len(nrow(fg)), function(i)
    data.table(Set = fg$pathway[i], Gene = fg$leadingEdge[[i]])))
  if (nrow(le) == 0) { message("No leading-edge genes for ", family); return(invisible(NULL)) }
  gtab <- le[, .(N_top_pathways = uniqueN(Set), Pathways = paste(unique(label_map(top$Set)[match(Set, top$Set)]), collapse = "; ")), by = Gene]
  gtab[, Rho_with_other_score := ranks[Gene]]
  setorder(gtab, -N_top_pathways, -Rho_with_other_score)
  fwrite(gtab, file.path(table_dir, paste0(stem, "_leading_edge_genes.tsv")), sep = "\t")

  shown <- head(gtab[order(-abs(Rho_with_other_score))], n_max_leading_genes)
  shown <- shown[order(-Rho_with_other_score)]            # ordering only; most positive at the top
  if (nrow(shown) < 2) { message("Too few leading-edge genes to draw for ", family); return(invisible(NULL)) }
  mat <- t(scale(t(expr[shown$Gene, col_order, drop = FALSE])))
  mat[!is.finite(mat)] <- 0; mat[mat > 2.5] <- 2.5; mat[mat < -2.5] <- -2.5
  lim <- max(0.3, ceiling(max(abs(shown$Rho_with_other_score)) * 10) / 10)
  ra <- rowAnnotation(
    `Rho with score` = shown$Rho_with_other_score,
    `Top pathways (n)` = anno_barplot(shown$N_top_pathways, width = unit(1.6, "cm"),
                                      gp = gpar(fill = "grey50", col = NA), border = FALSE),
    col = list(`Rho with score` = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B"))),
    annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(0.4, "cm"))
  ht <- Heatmap(mat, name = "Gene Z-score",
                col = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
                cluster_rows = FALSE, cluster_columns = FALSE, column_split = cl_vec,
                top_annotation = top_ann, right_annotation = ra,
                row_names_gp = gpar(fontsize = 9), show_column_names = FALSE, border = TRUE,
                column_title = paste0(family, " top pathways: leading-edge genes (rho with the ", other_name,
                                      " score; ", nrow(mat), " of ", nrow(gtab), " genes shown)"))
  save_ht(ht, paste0(stem, "_top_pathways_leading_edge_heatmap"), w = 12, h = max(6, 0.18 * nrow(mat) + 3.5))
}
leading_edge_heatmap(top_chr, p_clean$e, "Chromatin", "estrogen",  "chromatin")
leading_edge_heatmap(top_est, p_clean$c, "Estrogen",  "chromatin", "estrogen")

# ----------------------------
# 9. Random-family null (same construction as the observed test)
#    Random gene sets, size-matched to the real sets, are grouped into random
#    "estrogen" and "chromatin" families; PC1 and the same adjustments are applied.
# ----------------------------
pool_genes <- setdiff(rownames(expr), c(est_genes, chr_genes))
est_sz <- lengths(sets_clean[est_c])
chr_sz <- lengths(sets_clean[chr_c])
k_e <- length(est_c); k_c <- length(chr_c)
pool_e <- max(200, 3 * k_e); pool_c <- max(200, 3 * k_c)

make_random <- function(prefix, n_sets, sizes) {
  out <- lapply(seq_len(n_sets), function(i) sample(pool_genes, sample(sizes, 1)))
  names(out) <- paste0(prefix, seq_len(n_sets))
  out
}
set.seed(42)
rand_scores <- run_ssgsea(expr, c(make_random("RE", pool_e, est_sz), make_random("RC", pool_c, chr_sz)))
re_ids <- grep("^RE", rownames(rand_scores), value = TRUE)
rc_ids <- grep("^RC", rownames(rand_scores), value = TRUE)

set.seed(42)
null_scores <- lapply(seq_len(n_random), function(i) {
  list(e = family_pc1(rand_scores, sample(re_ids, k_e))$score,
       c = family_pc1(rand_scores, sample(rc_ids, k_c))$score)
})
null_E <- do.call(cbind, lapply(null_scores, `[[`, "e"))
null_C <- do.call(cbind, lapply(null_scores, `[[`, "c"))

null_rho <- function(vars) {
  rE <- apply(adjust(null_E, vars), 2, rank)
  rC <- apply(adjust(null_C, vars), 2, rank)
  vapply(seq_len(n_random), function(i) cor(rE[, i], rC[, i]), numeric(1))
}
null_dt <- rbindlist(lapply(c("Unadjusted", "Cluster_cellularity_TILs_adjusted"), function(m) {
  data.table(Model = m, Rho = null_rho(models[[m]]))
}))
obs_dt <- data.table(Model = c("Unadjusted", "Cluster_cellularity_TILs_adjusted"),
                     Observed = c(get_m("Unadjusted", "Spearman_rho"),
                                  get_m("Cluster_cellularity_TILs_adjusted", "Spearman_rho")))
emp <- rbindlist(lapply(obs_dt$Model, function(m) {
  obs <- obs_dt[Model == m, Observed]
  r   <- null_dt[Model == m, Rho]
  data.table(Model = m, Observed_rho = obs,
             Null_median_rho = median(r), Null_95th_percentile_abs_rho = quantile(abs(r), 0.95),
             Empirical_P = (1 + sum(abs(r) >= abs(obs))) / (length(r) + 1))
}))
print(emp)
fwrite(emp, file.path(table_dir, "random_null_empirical_P.tsv"), sep = "\t")

p_null <- ggplot(null_dt, aes(Rho)) +
  geom_histogram(bins = 40, fill = "grey75", color = "white") +
  geom_vline(data = obs_dt, aes(xintercept = Observed), color = "#B2182B", linewidth = 1) +
  facet_wrap(~ Model, ncol = 1) +
  labs(title = "Observed estrogen-chromatin correlation vs random gene-set families",
       subtitle = paste0(n_random, " random family pairs, size-matched; red line = observed"),
       x = "Spearman rho", y = "Random pairs") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot")
save_plot(p_null, "random_null_distribution", 7, 7)

# ----------------------------
# 10. Exploratory pathway-by-pathway matrix (adjusted residuals, BH FDR)
# ----------------------------
vars_full <- models[["Cluster_cellularity_TILs_adjusted"]]
R <- apply(adjust(t(scores_clean[c(chr_c, est_c), , drop = FALSE]), vars_full), 2, rank)
colnames(R) <- c(chr_c, est_c)
rho_m <- cor(R)[chr_c, est_c, drop = FALSE]

k_cov  <- ncol(model.matrix(~ ., data = cov_df[, vars_full, drop = FALSE])) - 1
df_res <- n - 2 - k_cov                        # approximate df for a partial correlation
tval   <- rho_m * sqrt(df_res / (1 - rho_m^2))
pval   <- 2 * pt(-abs(tval), df_res)

pair_dt <- data.table(Chromatin_set = rep(chr_c, times = length(est_c)),
                      Estrogen_set  = rep(est_c, each  = length(chr_c)),
                      Rho = as.vector(rho_m), P_value = as.vector(pval))
pair_dt[, FDR := p.adjust(P_value, method = "BH")]
pair_dt[, `:=`(Chromatin_label = label_map(Chromatin_set)[match(Chromatin_set, c(chr_c))],
               Estrogen_label  = label_map(Estrogen_set)[match(Estrogen_set, c(est_c))])]
setorder(pair_dt, FDR, P_value)
fwrite(pair_dt, file.path(table_dir, "pairwise_adjusted_correlations_exploratory.tsv"), sep = "\t")
cat("Pairs tested:", nrow(pair_dt), "| FDR < 0.05:", sum(pair_dt$FDR < 0.05), "\n")

lab_chr <- setNames(label_map(chr_c), chr_c)
lab_est <- setNames(label_map(est_c), est_c)
pair_dt[, `:=`(Chromatin_label = lab_chr[Chromatin_set], Estrogen_label = lab_est[Estrogen_set])]
co <- if (length(chr_c) > 2) hclust(dist(rho_m))$order else seq_along(chr_c)
eo <- if (length(est_c) > 2) hclust(dist(t(rho_m)))$order else seq_along(est_c)
pair_dt[, Chromatin_label := factor(Chromatin_label, levels = rev(lab_chr[chr_c[co]]))]
pair_dt[, Estrogen_label  := factor(Estrogen_label,  levels = lab_est[est_c[eo]])]

p_hm <- ggplot(pair_dt, aes(Estrogen_label, Chromatin_label, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman rho") +
  labs(title = "Exploratory: chromatin vs estrogen pathways",
       subtitle = "Residuals after adjusting for cluster, cellularity and TILs; shared genes removed; * FDR < 0.05",
       x = "Estrogen pathways", y = "Chromatin pathways") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), plot.title.position = "plot",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))
save_plot(p_hm, "pairwise_adjusted_heatmap", 13, max(8, 0.2 * length(chr_c) + 4))

# ----------------------------
# 11. README and provenance (fixed file names)
# ----------------------------
writeLines(
  c("pILC estrogen vs chromatin pathway correlation",
    "Seed: 42",
    "Gene sets: MSigDB GO (local C5 GMT), Hallmark (local GMT), Reactome and KEGG Medicus (msigdbr); family assignment by name rules (geneset_selection_audit.tsv).",
    paste0("Set size ", min_set_size, "-", max_set_size, " expressed genes; Jaccard >= ", jaccard_cutoff, " collapsed within family."),
    "Expression: saved VST matrix (expression-filtered genes, made before clustering), ssGSEA (GSVA).",
    "Primary: genes shared by both families removed; one PC1 score per family; Spearman on residuals; permutation P.",
    "Models: unadjusted; cluster; cluster + tumor cellularity + TILs (missing values filled with the cohort median).",
    "Permutation is inside clusters when the cluster is adjusted for.",
    "Null: size-matched random gene-set families processed identically.",
    "Per-cluster rho is descriptive (n about 10 per cluster). The pathway matrix is exploratory (BH FDR).",
    "Clusters were derived from the same expression data; results are not independent validation.",
    "RNA pathway scores show coordinated expression, not causal or mechanistic links."),
  file.path(out_dir, "README_chromatin_estrogen_correlation.txt"))
writeLines(c(paste0("msigdbr version: ", as.character(packageVersion("msigdbr"))),
             paste0("GSVA version: ", as.character(packageVersion("GSVA"))),
             paste0("C5 GMT: ", c5_gmt_file),
             paste0("C5 GMT md5: ", unname(tools::md5sum(c5_gmt_file))),
             paste0("Hallmark GMT: ", h_gmt_file),
             paste0("Hallmark GMT md5: ", unname(tools::md5sum(h_gmt_file)))),
           file.path(out_dir, "database_provenance.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done. Output: ", out_dir)