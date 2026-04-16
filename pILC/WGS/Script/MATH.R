##############################################
## Visualization of MATH Scores (ITH Analysis)
##############################################

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

#------------------------------------------------
# Input file
#------------------------------------------------
math_file <- "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WGS/PCGR/MAFtool/MATH_scores.tsv"

# Read data
math_df <- read_tsv(math_file)

# Clean sample names (optional)
math_df <- math_df %>%
  mutate(Sample = gsub("\\.taf\\.vcf\\.gz", "", Tumor_Sample_Barcode))

#------------------------------------------------
# Calculate summary statistics
#------------------------------------------------
median_math <- median(math_df$MATH_Score)
mean_math <- mean(math_df$MATH_Score)
sd_math <- sd(math_df$MATH_Score)

cat("Median MATH:", median_math, "\nMean:", mean_math, "\nSD:", sd_math, "\n")

# Label high/low ITH based on median split
math_df <- math_df %>%
  mutate(ITH_Group = ifelse(MATH_Score >= median_math, "High ITH", "Low ITH"))

#------------------------------------------------
# Boxplot with jitter (summary view)
#------------------------------------------------
p1 <- ggplot(math_df, aes(x = "", y = MATH_Score, fill = ITH_Group)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.1, size = 3, aes(color = ITH_Group)) +
  geom_hline(yintercept = median_math, linetype = "dashed", color = "red") +
  annotate("text", x = 1.2, y = median_math + 0.5,
           label = paste0("Median = ", round(median_math, 2)),
           color = "red", size = 4) +
  scale_fill_manual(values = c("High ITH" = "#E64B35", "Low ITH" = "#4DBBD5")) +
  scale_color_manual(values = c("High ITH" = "#E64B35", "Low ITH" = "#4DBBD5")) +
  theme_bw(base_size = 14) +
  labs(title = "Intra-tumor Heterogeneity (MATH Scores)",
       y = "MATH Score", x = "") +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

#------------------------------------------------
# Barplot ordered by MATH score
#------------------------------------------------
p2 <- ggplot(math_df %>% arrange(MATH_Score),
             aes(x = reorder(Sample, MATH_Score), y = MATH_Score, fill = ITH_Group)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = round(MATH_Score, 1)), vjust = -0.5, size = 3.5) +
  geom_hline(yintercept = median_math, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("High ITH" = "#E64B35", "Low ITH" = "#4DBBD5")) +
  theme_bw(base_size = 14) +
  labs(title = "MATH Scores Across pILC Samples",
       x = "Sample", y = "MATH Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "top")

#------------------------------------------------
# Save plots (PDF + PNG)
#------------------------------------------------
ggsave("MATH_Score_Boxplot.pdf", p1, width = 5, height = 5)
ggsave("MATH_Score_Barplot.pdf", p2, width = 7, height = 5)

ggsave("MATH_Score_Boxplot.png", p1, width = 5, height = 5, dpi = 300)
ggsave("MATH_Score_Barplot.png", p2, width = 7, height = 5, dpi = 300)

cat("Plots saved as PDF and PNG (300 dpi)\n")

library(patchwork)
combined <- p1 + p2 + plot_layout(ncol = 2)
ggsave("MATH_Score_Combined.pdf", combined, width = 10, height = 5)
ggsave("MATH_Score_combined.png", combined, width = 10, height = 5, dpi = 300)
