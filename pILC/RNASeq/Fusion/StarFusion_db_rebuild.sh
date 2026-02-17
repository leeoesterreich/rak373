#!/bin/bash
# Before rebuild indices, make sure to download (wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz)
# and tar -zxvf GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz
# rebuild_ctat_genome_lib_indices for StarFusion (https://github.com/STAR-Fusion/STAR-Fusion/wiki/rebuild_ctat_genome_lib_indices)

#SBATCH --job-name=rebuild_indices
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=rebuild_indices_%j.log

# Change to the ctat_genome_lib_build_dir
cd /ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/StarFusion_db/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir

# Rebuild indices
/software/rhel9/manual/install/star-fusion/1.15.1/python3.11/lib/STAR-Fusion/ctat-genome-lib-builder/util/rebuild_indices.pl .
