Standalone OncodriveCLUSTL 1.2.0 for pILC WES
================================================

Files
-----
1. install_oncodriveclustl_1.2.0.sh
   Creates an independent Python 3.10 environment, installs
   OncodriveCLUSTL 1.2.0, configures BGDATA_LOCAL, and preloads hg38.

2. OncodriveCLUSTL_Input.Rmd
   Creates the SNV-only mutation input and a dedicated CLUSTL coding-
   region file from the existing OncodriveFML GRCh38 CDS file.

3. OncodriveCLUSTL.sh
   SLURM batch script for the full coding OncodriveCLUSTL analysis.

Recommended order
-----------------

A. Install once on a login node:

   chmod +x install_oncodriveclustl_1.2.0.sh
   ./install_oncodriveclustl_1.2.0.sh

B. Run all chunks in OncodriveCLUSTL_Input.Rmd to create:

   .../OncodriveCLUSTL/input/pILC_WES_oncodriveclustl_SNV.tsv.gz
   .../OncodriveCLUSTL/reference/cds_hg38_oncodriveclustl.tsv.gz

C. Submit the analysis:

   chmod +x OncodriveCLUSTL.sh
   sbatch OncodriveCLUSTL.sh

Main parameters
---------------
Genome: hg38
Variants: coding SNVs only
Minimum mutations per element: 2
Minimum mutations per cluster: 2
Smoothing window: 11 bp
Cluster window: 11 bp
Context: trinucleotide
Simulations: 1,000
Simulation mode: region_restricted
Simulation window: 31 bp
Signature calculation: region_normalized
Coding segments: concatenated

Important distinction from OncodriveFML
---------------------------------------
OncodriveCLUSTL does not use CADD, an OncodriveFML configuration file,
tabix, or the FML indel settings. It uses SNVs, genomic element
coordinates, and the hg38 reference sequence.
