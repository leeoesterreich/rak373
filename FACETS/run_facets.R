# After generating the `tumor_normal.snp_pileup.gz` file using `snp-pileup`, run FACETS for purity, ploidy, and CNV analysis.

.libPaths("/ihome/alee/rak373/R/x86_64-pc-linux-gnu-library/4.3")
library(pctGCdata)
library(facets)

set.seed(1234) # Set seed for reproducibility

# Use your actual file path
datafile <- "/ix1/alee/LO_LAB/Personal/Rahul/PDO_CNV/Facets/sample.snp_pileup.gz"

# Check if the file exists
if (!file.exists(datafile)) stop("File not found!")

rcmat <- readSnpMatrix(datafile) # Read the SNP matrix

XX <- preProcSample(rcmat) # Run the FACETS preprocessing step

str(XX) # View the structure of the result

oo <- procSample(XX, cval = 150) # Perform segmentation

oo$dipLogR # Check the Estimated Diploid Log Ratio

fit <- emcncf(oo) # Run Copy Number Calling

head(fit$cncf) # Examine Segment Table
fit$purity
fit$ploidy

plotSample(oo, fit) # Visualize CNV Profile

pdf("sample_facets_plot.pdf", width = 10, height = 6)
plotSample(oo, fit)
dev.off()

png("sample_facets_plot.png", width = 1000, height = 600)
plotSample(oo, fit)
dev.off()

write.csv(fit$cncf, "sample_facets_segments.csv", row.names = FALSE) 

