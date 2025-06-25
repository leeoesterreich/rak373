CNV Analysis Using FACETS This pipeline uses FACETS in Normal-Tumor mode for copy number variation (CNV) analysis.

Before running FACETS, compile the snp-pileup binary from the FACETS source:

    g++ -std=c++11 -I/ihome/alee/rak373/htslib \
    -L/ihome/alee/rak373/htslib -lhts \
    -Wl,-rpath=/ihome/alee/rak373/htslib \
    -o snp-pileup snp-pileup.cpp
    
./snp-pileup --help  # Test installation
Steps

Tool mskcc/facets-suite (https://github.com/mskcc/facets-suite/tree/master)

VCF File

 Download the SNP VCF and its index file from the NCBI SNP database:
    wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz 
    wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz.tbi 
Input Preparation From sorted BAM files, generate tumor_normal.snp_pileup.gz using the snp-pileup-wrapper.R script provided in the FACETS suite.

Run FACETS by using the instructions provided over https://github.com/mskcc/facets/blob/master/vignettes/FACETS.pdf Alternatively, we can run FACETS using the run-facets-wrapper.R script provided in the FACETS suite, with a slight modification to the script. In run-facets-wrapper.R, use facets::readSnpMatrix instead of facetsSuite version called read_snp_matrix
