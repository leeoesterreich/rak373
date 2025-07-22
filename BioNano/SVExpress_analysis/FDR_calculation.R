```r
# ------------------------------------------------------------------------
# To check the number of rows and columns, along with column headers across each file
# ------------------------------------------------------------------------
 
sv_categories <- list(
  "0-100Kb_downstream"     = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Coordinate_Proximity_based/Downstream_100Kb/lm_test_results_significant_genes.txt",
  "0-10Kb_downstream"      = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Coordinate_Proximity_based/Downstream_10Kb/lm_test_results_significant_genes.txt",
  "gene_body"              = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/OverlapGenes_based/lm_test_results_significant_genes.txt",
  "0-10Kb_upstream"        = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Coordinate_Proximity_based/Upstream_10Kb/lm_test_results_significant_genes.txt",
  "0-100Kb_upstream"       = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Coordinate_Proximity_based/Upstream_100Kb/lm_test_results_significant_genes.txt",
  "1Mb_up_or_downstream"   = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Coordinate_Proximity_based/Up_Down_1Mb/lm_test_results_significant_genes.txt",
  "Distance_weighted_1Mb"   ="/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Distance_weighted/lm_test_results_significant_genes.txt"
)

for (category in names(sv_categories)) {
  cat("====== Category:", category, "======\n")
  
  file_path <- sv_categories[[category]]
  
  if (file.exists(file_path)) {
    df <- tryCatch(read.delim(file_path, header = TRUE, stringsAsFactors = FALSE), 
                   error = function(e) {
                     cat("Error reading file:", e$message, "\n")
                     return(NULL)
                   })
    
    if (!is.null(df)) {
      cat("Number of Rows:", nrow(df), "\n")
      cat("Number of Columns:", ncol(df), "\n")
      cat("Column Names:\n")
      print(colnames(df))
      cat("First Few Rows:\n")
      print(head(df, 5))
    }
  } else {
    cat("File does not exist:", file_path, "\n")
  }
  cat("\n")
}

# ------------------------------
# Data integrity checks
# ------------------------------
for (category in names(sv_categories)) {
  cat("====== Checking Category:", category, "======\n")
  file_path <- sv_categories[[category]]
  
  if (file.exists(file_path)) {
    df <- tryCatch(read.delim(file_path, header = TRUE, stringsAsFactors = FALSE), 
                   error = function(e) {
                     cat("Error reading file:", e$message, "\n")
                     return(NULL)
                   })
    
    if (!is.null(df)) {
      # Check for NA values
      na_counts <- colSums(is.na(df))
      total_na <- sum(na_counts)
      
      cat("Total NA values:", total_na, "\n")
      if (total_na > 0) {
        cat("Columns with NAs:\n")
        print(na_counts[na_counts > 0])
      }

      # Check if numeric columns contain non-numeric entries
      invalid_format_summary <- sapply(df, function(col) {
        if (is.character(col)) return(TRUE)  # definitely non-numeric
        if (!is.numeric(col)) {
          suppressWarnings(!all(!is.na(as.numeric(as.character(col)))))
        } else {
          return(FALSE)
        }
      })
      
      cat("Columns with possible invalid numeric formats:\n")
      print(names(df)[invalid_format_summary])

      if (total_na > 0 || any(invalid_format_summary)) {
        cat("Example rows with issues:\n")
        problem_rows <- df[apply(is.na(df) | df == "", 1, any), ]
        print(head(problem_rows, 3))
      } else {
        cat("No NA values or format issues detected.\n")
      }
    }
  } else {
    cat("File not found:", file_path, "\n")
  }
  cat("\n")
}

# ------------------------------------------------------------
# To compute FDR-adjusted p-values per category and add two new columns to each dataset
# ------------------------------------------------------------
# Create a list to store updated data frames with FDR
sv_data_with_fdr <- list()

for (category in names(sv_categories)) {
  cat("====== Computing FDR for:", category, "======\n")
  file_path <- sv_categories[[category]]
  
  if (file.exists(file_path)) {
    df <- tryCatch(read.delim(file_path, header = TRUE, stringsAsFactors = FALSE), 
                   error = function(e) {
                     cat("Error reading file:", e$message, "\n")
                     return(NULL)
                   })
    
    if (!is.null(df)) {
      # Compute FDR for both p-value columns
      df$no.correction.FDR <- p.adjust(df$no.correction..p.value., method = "fdr")
      df$CNA.correction.FDR <- p.adjust(df$CNA.correction..p.value., method = "fdr")
      
      sv_data_with_fdr[[category]] <- df
      
      # Preview
      print(head(df[, c("no.correction..p.value.", "no.correction.FDR",
                        "CNA.correction..p.value.", "CNA.correction.FDR")]))
    }
  } else {
    cat("File not found:", file_path, "\n")
  }
  cat("\n")
}
# Define output directory
output_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each category and write to CSV
for (category in names(sv_data_with_fdr)) {
  df <- sv_data_with_fdr[[category]]
  
  # Create a safe file name
  file_name <- paste0(gsub("[^a-zA-Z0-9]", "_", category), "_with_FDR.csv")
  file_path <- file.path(output_dir, file_name)
  
  # Write CSV
  write.csv(df, file = file_path, row.names = TRUE)
  
  cat("Exported:", file_path, "\n")
}

# ------------------------------------------------------------
# To filter each model separately, 10% FDR threshold (i.e., FDR < 0.1)
# ------------------------------------------------------------
# Create two new lists to store filtered results
filtered_no_correction <- list()
filtered_CNA_correction <- list()

# Threshold
fdr_threshold <- 0.1

for (category in names(sv_data_with_fdr)) {
  df <- sv_data_with_fdr[[category]]
  
  # Filter separately
  no_corr_filtered <- subset(df, no.correction.FDR < fdr_threshold)
  cna_corr_filtered <- subset(df, CNA.correction.FDR < fdr_threshold)
  
  # Store
  filtered_no_correction[[category]] <- no_corr_filtered
  filtered_CNA_correction[[category]] <- cna_corr_filtered
  
  # Print summary
  cat("=== Category:", category, "===\n")
  cat("No Correction < 0.1:", nrow(no_corr_filtered), "genes\n")
  cat("CNA Correction < 0.1:", nrow(cna_corr_filtered), "genes\n\n")
}

# Output folder
filtered_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Filtered_FDR_0.1"
if (!dir.exists(filtered_dir)) dir.create(filtered_dir, recursive = TRUE)

# Save both filtered sets
for (category in names(filtered_no_correction)) {
  # No correction
  write.csv(filtered_no_correction[[category]],
            file = file.path(filtered_dir, paste0(gsub("[^a-zA-Z0-9]", "_", category), "_no_correction_FDR_lt_0.1.csv")),
            row.names = TRUE)
  
  # CNA correction
  write.csv(filtered_CNA_correction[[category]],
            file = file.path(filtered_dir, paste0(gsub("[^a-zA-Z0-9]", "_", category), "_CNA_correction_FDR_lt_0.1.csv")),
            row.names = TRUE)
}
##########################################################
# To counts number of genes based on both:
# FDR < 0.1 vs ≥ 0.1; p-value < 0.05 vs ≥ 0.05, for both no.correction and CNA.correction columns.

sv_categories <- list(
  "0-100Kb_downstream"     = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/0_100Kb_downstream_with_FDR.csv",
  "gene_body"              = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/gene_body_with_FDR.csv",
  "0-100Kb_upstream"       = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/0_100Kb_upstream_with_FDR.csv",
  "1Mb_up_or_downstream"   = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/1Mb_up_or_downstream_with_FDR.csv",
  "distance_weighted"      = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Distance_weighted_1Mb_with_FDR.csv"
)
# Set p-value threshold
pval_thresh <- 0.05

# Print FDR and p-value stats per category
for (category in names(sv_categories)) {
  file_path <- sv_categories[[category]]

  cat("===== Stats for:", category, "=====\n")

  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)

    # FDR counts
    no_corr_fdr_lt_0.1 <- sum(df$no.correction.FDR < 0.1, na.rm = TRUE)
    no_corr_fdr_ge_0.1 <- sum(df$no.correction.FDR >= 0.1, na.rm = TRUE)

    cna_corr_fdr_lt_0.1 <- sum(df$CNA.correction.FDR < 0.1, na.rm = TRUE)
    cna_corr_fdr_ge_0.1 <- sum(df$CNA.correction.FDR >= 0.1, na.rm = TRUE)

    # p-value counts
    no_corr_p_lt_thresh <- sum(df$no.correction..p.value. < pval_thresh, na.rm = TRUE)
    no_corr_p_ge_thresh <- sum(df$no.correction..p.value. >= pval_thresh, na.rm = TRUE)

    cna_corr_p_lt_thresh <- sum(df$CNA.correction..p.value. < pval_thresh, na.rm = TRUE)
    cna_corr_p_ge_thresh <- sum(df$CNA.correction..p.value. >= pval_thresh, na.rm = TRUE)

    # Output
    cat("FDR < 0.1 (no_correction):     ", no_corr_fdr_lt_0.1, "\n")
    cat("FDR ≥ 0.1 (no_correction):     ", no_corr_fdr_ge_0.1, "\n")
    cat("FDR < 0.1 (CNA_correction):    ", cna_corr_fdr_lt_0.1, "\n")
    cat("FDR ≥ 0.1 (CNA_correction):    ", cna_corr_fdr_ge_0.1, "\n")

    cat("p < 0.05 (no_correction):      ", no_corr_p_lt_thresh, "\n")
    cat("p ≥ 0.05 (no_correction):      ", no_corr_p_ge_thresh, "\n")
    cat("p < 0.05 (CNA_correction):     ", cna_corr_p_lt_thresh, "\n")
    cat("p ≥ 0.05 (CNA_correction):     ", cna_corr_p_ge_thresh, "\n")
  } else {
    cat("File not found:", file_path, "\n")
  }

  cat("\n")
}
#####################################################################################################################################################################
# Identify breast cancer–associated genes (from TCGA/ICGC OncoVar sets) that intersect with SV-associated genes classified by different genomic proximity categories.
# ----------------------------------
# Based on FDR
# ----------------------------------
# Load required libraries
library(data.table)
library(dplyr)
library(readr)

# ----------------------------------
# Load and combine BRCA gene sets
# ----------------------------------
tcga_df <- fread("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Reference/TCGA.BRCA.onco.genes.OncoVar.tsv") %>%
  mutate(Source = "TCGA")

icgc_df <- fread("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Reference/ICGC.BRCA.onco.genes.OncoVar.tsv") %>%
  mutate(Source = "ICGC")

# Combine BRCA genes from both sources
brca_genes <- unique(c(tcga_df$Gene_symbol, icgc_df$Gene_symbol))

# ----------------------------------
# Define SV category file paths
# ----------------------------------
sv_categories <- list(
  "0-100Kb_downstream"     = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Filtered_FDR_0.1/0_100Kb_downstream_CNA_correction_FDR_lt_0.1.csv",
  "gene_body"              = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Filtered_FDR_0.1/gene_body_CNA_correction_FDR_lt_0.1.csv",
  "0-100Kb_upstream"       = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Filtered_FDR_0.1/0_100Kb_upstream_CNA_correction_FDR_lt_0.1.csv",
  "1Mb_up_or_downstream"   = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Filtered_FDR_0.1/1Mb_up_or_downstream_CNA_correction_FDR_lt_0.1.csv",
  "distance_weighted"      = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/Results_with_FDR/Filtered_FDR_0.1/Distance_weighted_1Mb_CNA_correction_FDR_lt_0.1.csv"
)

# Initialize result list
all_results <- list()

for (category in names(sv_categories)) {
  file_path <- sv_categories[[category]]

  if (file.exists(file_path)) {
    # Read gene names from row names
    df <- read.csv(file_path, row.names = 1, stringsAsFactors = FALSE)

    # Get row names (gene names)
    matched_genes <- rownames(df)[rownames(df) %in% brca_genes]

    if (length(matched_genes) > 0) {
      intersected_df <- df[matched_genes, , drop = FALSE]
      intersected_df$Gene <- matched_genes
      intersected_df$Proximity_Category <- category
      all_results[[category]] <- intersected_df
    }
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# Combine all matches
final_df <- bind_rows(all_results)

# Save to file
output_file <- "BRCA_intersection_with_FDR_significant_genes_CNA_correction.txt"
write.table(final_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Print summary
cat("Intersected genes saved to:", output_file, "\n")
cat("Total intersecting BRCA genes:", length(unique(final_df$Gene)), "\n")

########################################################################################################################
# There are 13 BRCA unique genes that are significantly associated with expression (FDR<0.1)
# To Filter SV_df_unique.csv by BRCA Genes in OverlapGenes or NearestNonOverlapGene
# === Load required packages ===
library(readr)
library(dplyr)
library(stringr)

# === File paths ===
sv_path <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/SV_df_unique.csv"
output_path <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVExpress/SV_df_filtered_by_BRCA_significant_gene_list.csv"

# === Your gene list ===
gene_list <- c(
  "ZFP36L1", "APOBEC3B", "ERBB2", "ESR1", "KRAS", "NCOR1", "NF1", "PSIP1", "RB1",
  "AKT1", "BACH2", "ETV6", "KDM6A"
)
gene_list <- unique(gene_list)

# === Load SV file ===
sv_df <- read_csv(sv_path, show_col_types = FALSE)

# === Helper function to check if any gene in OverlapGenes matches gene_list ===
check_overlap_genes <- function(overlap_str) {
  genes <- str_split(as.character(overlap_str), ";")[[1]] %>% str_trim()
  any(genes %in% gene_list)
}

# === Filter rows ===
filtered_sv <- sv_df %>%
  rowwise() %>%
  filter(
    check_overlap_genes(OverlapGenes) ||
    str_trim(as.character(NearestNonOverlapGene)) %in% gene_list
  ) %>%
  ungroup()

# === Save result ===
write_csv(filtered_sv, output_path)
cat("Filtered rows saved to:\n", output_path, "\n")

# === Identify which genes from gene_list were NOT mapped ===
mapped_genes <- c(
  unlist(str_split(na.omit(sv_df$OverlapGenes), ";")),
  na.omit(sv_df$NearestNonOverlapGene)
) %>%
  str_trim() %>%
  unique()

unmapped_genes <- setdiff(gene_list, mapped_genes)

cat("Genes NOT found in OverlapGenes or NearestNonOverlapGene:\n")
print(unmapped_genes)

########################################################################################################################
