# ==============================================================================
# Relabel an existing R-NMF membership file (k = 3) to the old GenePattern numbers
#
# Why: the groups from R-NMF and GenePattern are identical (ARI = 1), but NMF
# numbers its clusters arbitrarily. Earlier analyses used the GenePattern
# numbers, so the R labels are converted to them. No NMF or VST rerun is needed.
#
# Input : - R-NMF membership CSV (sample IDs in column 1 + a cluster column)
#         - GenePattern consensus.k.3.gct (tab-separated membership list):
#             line 1: #1.2 | line 2: 31 <TAB> 1 |
#             line 3: Name, Description, membership.ordered | then one line/sample
# Output: Cluster_membership_k3_GPlabels.csv, in the same folder as the input file
#         (the original Cluster_membership_k3.csv is kept unchanged), with exactly two columns
#         and the original sample names (dot format), GenePattern cluster numbers:
#             "sample","cluster"
#             "TP03.117",1
#         The .gct is not changed.
# ==============================================================================

library(mclust)   # adjustedRandIndex()

# ---- 0. Settings (edit here only) --------------------------------------------
r_file  <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/NMF/R_NMFConsensus_50runs/Cluster_membership_k3.csv"
gp_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/Archive/NMF/GenePattern_NMFConsensus/vst_transformed_data.consensus.k.3.gct"
cluster_col <- "cluster"      # name of the label column in r_file (check the head() print)

# Helper: make sample names comparable between files
#   "TP13.1150_S40" -> "TP13-1150"   (remove _S## suffix, dots -> hyphens)
clean_id <- function(x) gsub("\\.", "-", sub("_S[0-9]+$", "", x))

# ---- 1. Read the existing R-NMF membership file (CSV) ------------------------
mem <- read.csv(r_file, check.names = FALSE, stringsAsFactors = FALSE)
print(head(mem))                                   # check the column names here
stopifnot(cluster_col %in% colnames(mem))          # stop if the label column is not found

r_lab <- setNames(as.integer(mem[[cluster_col]]),
                  clean_id(mem[[1]]))              # first column = sample names

# ---- 2. Read the old GenePattern clusters (.gct, tab-separated, NOT csv) -----
gp <- read.table(gp_file, skip = 2, header = TRUE, sep = "\t",
                 check.names = FALSE, quote = "", comment.char = "")
stopifnot("membership.ordered" %in% colnames(gp))

gp_lab <- setNames(as.integer(gp$membership.ordered), clean_id(gp$Name))
cat("GenePattern cluster sizes:\n"); print(table(gp_lab))   # expect 10, 10, 11

# ---- 3. Check that both files contain the same samples -----------------------
stopifnot(setequal(names(r_lab), names(gp_lab)))
gp_lab <- gp_lab[names(r_lab)]                     # same sample order as the R file

# ---- 4. Mapping from the data (R label -> GenePattern label) -----------------
# For each R cluster, take the GenePattern cluster that shares most samples.
overlap <- table(R = r_lab, GP = gp_lab)
print(overlap)                                     # one non-zero cell per row expected

map <- apply(overlap, 1, function(r) as.integer(colnames(overlap)[which.max(r)]))
stopifnot(!anyDuplicated(map))                     # must be one-to-one
cat("Mapping (R label -> GenePattern label):\n"); print(map)

new_lab <- unname(map[as.character(r_lab)])

# ---- 5. Safety checks: the groups must be identical to GenePattern -----------
stopifnot(sum(new_lab != gp_lab) == 0)             # no sample may change group
cat("ARI:", round(adjustedRandIndex(gp_lab, new_lab), 3), "\n")   # expect 1
print(table(new_lab))                                              # expect 10, 10, 11

# ---- 6. Save (same folder, same file name, no date) --------------------------
# Output has exactly two columns, "sample" and "cluster". The sample names are
# written exactly as they were in the input file (mem[[1]], e.g. TP03.117),
# and the row order is unchanged.
out <- data.frame(sample = mem[[1]], cluster = new_lab, stringsAsFactors = FALSE)
# The original file (r_file) is NOT changed. The relabelled file is written next to it.
out_file <- file.path(dirname(r_file), "Cluster_membership_k3_GPlabels.csv")
write.csv(out, out_file, row.names = FALSE)        # quoted text, plain numbers
cat("Saved:", out_file, "\n")
print(head(out))                                   # check the structure
stopifnot(identical(colnames(out), c("sample", "cluster")))
cat("Original kept unchanged:", r_file, "\n")
