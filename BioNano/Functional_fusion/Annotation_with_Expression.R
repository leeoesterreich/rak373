# This R script is used to identify functional fusions from the list of fusions detected by OGM (Bionano).

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)

##### Prepare Fusion df from SV data #####
SV_df <- read.csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/SV_df_unqiue_pILC.csv", stringsAsFactors = FALSE)
fusions = data.frame(fusions = unique(SV_df$PutativeGeneFusion))
fusions = fusions$fusions[-c(1)]

fusion_df = data.frame(fusions = fusions, samples = SV_df[match(fusions, SV_df$PutativeGeneFusion), "sample"], 
                       type = SV_df[match(fusions, SV_df$PutativeGeneFusion), "Type"],
                       geneA = sapply(strsplit(fusions, split = "-"), FUN = function(x) x[1]),
                       chrA = SV_df[match(fusions, SV_df$PutativeGeneFusion), "RefcontigID1"],
                       StartA = SV_df[match(fusions, SV_df$PutativeGeneFusion), "RefStartPos"],
                       EndA = SV_df[match(fusions, SV_df$PutativeGeneFusion), "RefStartPos"],
                       geneB = sapply(strsplit(fusions, split = "-"), FUN = function(x) x[2]),
                       chrB = SV_df[match(fusions, SV_df$PutativeGeneFusion), "RefcontigID2"],
                       StartB = SV_df[match(fusions, SV_df$PutativeGeneFusion), "RefEndPos"],
                       EndB = SV_df[match(fusions, SV_df$PutativeGeneFusion), "RefEndPos"])



# split complex fusions
`%notin%` <- Negate(`%in%`)

for (fus in grep(x = fusion_df$fusions, pattern = ";", value = TRUE)) {
  temp <- subset(fusion_df, fusions %in% fus)
  fusions_ <- unlist(strsplit(temp$fusions, split = ";"))

  temp_ <- as.data.frame(matrix(data = NA, nrow = length(fusions_), ncol = ncol(temp)))
  colnames(temp_) <- colnames(temp)
  temp_[1:length(fusions_), ] <- temp
  temp_$fusions <- fusions_

  fusion_df <- subset(fusion_df, fusions %notin% fus)
  fusion_df <- rbind(fusion_df, temp_)
}

fusion_df$geneA = sapply(strsplit(fusion_df$fusions, split = "-"), FUN = function(x) x[1])
fusion_df$geneB = sapply(strsplit(fusion_df$fusions, split = "-"), FUN = function(x) x[2])

fusion_df$fusions_inv = paste0(fusion_df$geneB, "-", fusion_df$geneA)
table(fusion_df$samples) # Counts number of fusions per sample

# Check missing samples 
table(SV_df$sample)  # all samples in original data
setdiff(unique(SV_df$sample), unique(fusion_df$samples))

##### Annotates gene expression levels per fusion partner and classifies up/down/neutral expression #####
expr_data <- read_csv("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/pILC_log2cpm_matrix.csv")

# Keep Ensembl as reference, use gene names as rownames
expr_data_clean <- expr_data |>
  filter(!is.na(external_gene_name) & external_gene_name != "") |>  # optional
  distinct(external_gene_name, .keep_all = TRUE) |>                 # drop duplicates if any
  column_to_rownames(var = "external_gene_name") |>
  select(-ensembl_gene_id)

# If above command gives error like : Error: unable to find an inherited method for function â€˜selectâ€™ for signature â€˜x = "data.frame"â€™
# then replace this block with the below command and the rest will be same
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expr_data_clean <- expr_data |>                                           
  dplyr::filter(!is.na(external_gene_name) & external_gene_name != "") |>
  dplyr::distinct(external_gene_name, .keep_all = TRUE) |>
  dplyr::select(-ensembl_gene_id) |>
  tibble::column_to_rownames(var = "external_gene_name")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
# Just to check the matrix of gene expression
expr_data_clean["TSPAN6", ] #Expression of "TSPAN6" across samples

# Z-Score Normalization
expr_zscore <- t(scale(t(as.matrix(expr_data_clean))))

# Check for Overexpression of a Gene 
gene <- "CDH3" # Example CDH3
gene_expr <- expr_data_clean[gene, ] |> unlist() |> as.numeric()
names(gene_expr) <- colnames(expr_data_clean)
summary(gene_expr)
sd(gene_expr)
upregulated <- gene_expr > (median(gene_expr) + 0.4 * sd(gene_expr))
downregulated <- gene_expr < (median(gene_expr) - 0.4 * sd(gene_expr))

# Print sample names
names(which(upregulated))
names(which(downregulated))

# Annotate expr_data with fusion_df
fusion_df <- fusion_df |>
  rowwise() |>
  mutate(
    geneA_exp = if (geneA %in% rownames(expr_data_clean) && samples %in% colnames(expr_data_clean)) {
      expr_data_clean[geneA, samples]
    } else {
      NA_real_
    },
    geneB_exp = if (geneB %in% rownames(expr_data_clean) && samples %in% colnames(expr_data_clean)) {
      expr_data_clean[geneB, samples]
    } else {
      NA_real_
    }
  ) |>
  ungroup()

get_expr_status <- function(gene, sample) {
  if (!(gene %in% rownames(expr_data_clean)) || !(sample %in% colnames(expr_data_clean))) {
    return(NA_character_)
  }
  expr_vec <- expr_data_clean[gene, ] |> unlist() |> as.numeric()
  expr_val <- expr_data_clean[gene, sample]
  names(expr_vec) <- colnames(expr_data_clean)
  
  up_thresh <- median(expr_vec) + 0.4 * sd(expr_vec)
  down_thresh <- median(expr_vec) - 0.4 * sd(expr_vec)
  
  if (expr_val > up_thresh) return("up")
  if (expr_val < down_thresh) return("dn")
  return("normal")
}

fusion_df <- fusion_df |>
  rowwise() |>
  mutate(
    geneA_status = get_expr_status(geneA, samples),
    geneB_status = get_expr_status(geneB, samples)
  ) |>
  ungroup()

write.csv(
  fusion_df,
  file = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/fusion_expression_annotation.csv",
  row.names = FALSE
)

##### Combining expression status with functional annotation from OncoVar #####
oncoVar <- read.delim("/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Reference/ICGC.BRCA.all.genes.OncoVar.tsv", sep = "\t", stringsAsFactors = FALSE)
rownames(oncoVar) <- oncoVar$Gene_symbol

# Annotate Each Fusion Gene (A & B) with OncoVar Status
fusion_df$geneA_OncoTS <- ifelse(
  grepl("Oncogene", oncoVar[fusion_df$geneA, "OncoKB"], ignore.case = TRUE), "Oncogene",
  ifelse(grepl("TSG", oncoVar[fusion_df$geneA, "OncoKB"], ignore.case = TRUE), "TSG", NA)
)

fusion_df$geneB_OncoTS <- ifelse(
  grepl("Oncogene", oncoVar[fusion_df$geneB, "OncoKB"], ignore.case = TRUE), "Oncogene",
  ifelse(grepl("TSG", oncoVar[fusion_df$geneB, "OncoKB"], ignore.case = TRUE), "TSG", NA)
)

# Summarize Fusion Expression Status
fusion_df$exp <- ifelse(
  fusion_df$geneA_status == "up" | fusion_df$geneB_status == "up", "up-regulation",
  ifelse(fusion_df$geneA_status == "dn" | fusion_df$geneB_status == "dn", "down-regulation", "-")
)

write.csv(
  fusion_df,
  file = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/fusion_expression_OncoVar_annotated.csv",
  row.names = FALSE
)

# Filter Expression-Relevant Fusions
functional_fusions <- fusion_df %>% 
  filter(exp != "-")
  
write.csv(
  functional_fusions,
  file = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/pILC/functional_fusions_expression_only.csv",
  row.names = FALSE
)

##### chromosomal enrichment analysis #####

# Total expression-relevant fusions involving each chromosome
chr_focus <- c("chr1", "chr8", "chr11", "chr17")

fusion_chr_counts <- table(
  c(
    functional_fusions$chrA[functional_fusions$chrA %in% chr_focus],
    functional_fusions$chrB[functional_fusions$chrB %in% chr_focus]
  )
)
fusion_chr_counts

# Compare to Background (All fusions)
all_chr_counts <- table(
  c(
    fusion_df$chrA[fusion_df$chrA %in% chr_focus],
    fusion_df$chrB[fusion_df$chrB %in% chr_focus]
  )
)
all_chr_counts

# Statistical Enrichment
chr <- "chr1" # Example for chr17
exp_chr <- fusion_chr_counts[chr]
nonexp_chr <- all_chr_counts[chr] - exp_chr
exp_other <- sum(fusion_chr_counts) - exp_chr
nonexp_other <- sum(all_chr_counts) - exp_chr - nonexp_chr

contingency_matrix <- matrix(
  c(exp_chr, nonexp_chr, exp_other, nonexp_other),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Expression = c("Relevant", "NotRelevant"),
    Location = c(paste0(chr), "Other")
  )
)

print(contingency_matrix)
fisher.test(contingency_matrix)

# Loop over all chromosomes of interest:
results <- lapply(c("chr1", "chr8", "chr11", "chr17"), function(chr) {
  exp_chr <- fusion_chr_counts[chr]
  nonexp_chr <- all_chr_counts[chr] - exp_chr
  exp_other <- sum(fusion_chr_counts) - exp_chr
  nonexp_other <- sum(all_chr_counts) - exp_chr - nonexp_chr
  
  mat <- matrix(c(exp_chr, nonexp_chr, exp_other, nonexp_other), nrow = 2)
  ft <- fisher.test(mat)
  
  data.frame(
    chromosome = chr,
    p_value = ft$p.value,
    odds_ratio = ft$estimate,
    lower_CI = ft$conf.int[1],
    upper_CI = ft$conf.int[2]
  )
}) |> bind_rows()
print(results)

# Save a Summary Table and create Plot
 
chr_fusion_enrichment_plot <- ggplot(results, aes(x = chromosome, y = odds_ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(
    aes(ymin = lower_CI, ymax = upper_CI),
    width = 0.2
  ) +
  labs(
    title = "Odds Ratio of Expression-Relevant Fusions by Chromosome",
    y = "Odds Ratio (Fisher Test)",
    x = "Chromosome"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = function(x) {
      b <- pretty(x)
      b[b != 0]  # drop 0 tick & label only
    }
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# Save the summary table
write.csv(
  results,
  file = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Functional_Fusion/chr_fusion_enrichment_summary.csv",
  row.names = FALSE
)

# Save the plot
ggsave(
  filename = "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Functional_Fusion/chr_fusion_enrichment_plot.png",
  plot = chr_fusion_enrichment_plot,
  width = 8,
  height = 5,
  dpi = 300
)

