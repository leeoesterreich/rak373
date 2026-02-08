############################################
## Identify Repeated Fusions That Are Functional
## --------------------------------------------
## This script identifies fusions that occur more than once in the 
## fusion count dataset and checks if they are also present in the
## annotated functional fusion file.
############################################

# Load libraries
library(dplyr)

# ====== Path ======
fusion_count_path <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/SVs/fusion_counts_IPM.csv"
functional_fusion_path <- "/ix1/alee/LO_LAB/Personal/Rahul/BioNano/Functional_Fusion/IPM/functional_fusions_expression_only_IPM.csv"

# ====== Load fusion count file ======
fusion_counts <- read.csv(fusion_count_path, stringsAsFactors = FALSE)

# Ensure columns exist
stopifnot("Count" %in% colnames(fusion_counts))
stopifnot("Fusion" %in% colnames(fusion_counts))

# ====== Filter fusions with Count > 1 ======
repeated_fusions <- fusion_counts %>%
  filter(Count > 1)

cat("Fusions with Count > 1:\n")
print(repeated_fusions)

# ====== Load functional fusion file ======
functional_fusions <- read.csv(functional_fusion_path, stringsAsFactors = FALSE)

# Check for column
stopifnot("fusions" %in% colnames(functional_fusions))

# ====== Check overlaps ======
present_in_functional <- repeated_fusions %>%
  filter(Fusion %in% functional_fusions$fusions)

not_present <- setdiff(repeated_fusions$Fusion, functional_fusions$fusions)

# ====== Output results ======

cat("\nRepeated fusions that ARE present in functional_fusions:\n")
print(present_in_functional)

cat("\nRepeated fusions that are NOT present in functional_fusions:\n")
print(not_present)

# ====== Save results ======

write.csv(repeated_fusions, "repeated_fusions.csv", row.names = FALSE)
write.csv(present_in_functional, "repeated_fusions_present_in_functional.csv", row.names = FALSE)
write.csv(data.frame(Fusion = not_present), "repeated_fusions_NOT_in_functional.csv", row.names = FALSE)

cat("\nResults saved to CSV files in current working directory.\n")
