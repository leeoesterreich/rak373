#This script is used to download human CollecTRI table.

cd /ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/decoupler/

# Download the human CollecTRI table
wget https://static.omnipathdb.org/tables/collectri-human.tsv.gz
# (if wget is not available: curl -L -O https://static.omnipathdb.org/tables/collectri-human.tsv.gz)

# Check that the download is complete (listed size is about 1.3 MB)
ls -lh collectri-human.tsv.gz
gunzip -t collectri-human.tsv.gz && echo "gzip OK"

# Look at the column names and the first rows
zcat collectri-human.tsv.gz | head -5
zcat collectri-human.tsv.gz | wc -l