#!/bin/bash
#SBATCH --job-name=CADD_hg38
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveFML/slurm_log/CADD_hg38_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON/OncodriveFML/slurm_log/CADD_hg38_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -Eeuo pipefail

CADD_DIR="/ix1/alee/LO_LAB/Personal/Rahul/tools/OncodriveFML/resources/CADD/v1.7/GRCh38"

CADD_URL="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
INDEX_URL="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz.tbi"
MD5_URL="https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/MD5SUMs"

mkdir -p "$CADD_DIR"
cd "$CADD_DIR"

echo "Download started: $(date)"
echo "Destination: $CADD_DIR"

df -h "$CADD_DIR"

# -c allows an interrupted download to resume.
wget -c "$CADD_URL" \
    -O whole_genome_SNVs.tsv.gz

wget -c "$INDEX_URL" \
    -O whole_genome_SNVs.tsv.gz.tbi

wget -c "$MD5_URL" \
    -O MD5SUMs

echo
echo "Downloaded files:"
ls -lh \
    whole_genome_SNVs.tsv.gz \
    whole_genome_SNVs.tsv.gz.tbi \
    MD5SUMs

echo
echo "Verifying MD5 checksums..."

grep -E \
    'whole_genome_SNVs\.tsv\.gz($| )|whole_genome_SNVs\.tsv\.gz\.tbi($| )' \
    MD5SUMs \
    > MD5SUMs_SNV_only

md5sum -c MD5SUMs_SNV_only

echo
echo "Download and verification completed: $(date)"
