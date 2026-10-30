#!/bin/bash
#SBATCH --job-name=rnafusion_refs_dl
#SBATCH --clusters=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/refs_dl_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/refs_dl_%j.err

set -eu
module load awscli/2.28.6

DST=/ix1/alee/LO_LAB/General/References/rnafusion_4.1.3
SRC=s3://nf-core-awsmegatests/rnafusion/references/

# More parallel transfers, set in a temporary config so ~/.aws is not touched
export AWS_CONFIG_FILE=$(mktemp)
printf "[default]\ns3 =\n  max_concurrent_requests = 16\n" > ${AWS_CONFIG_FILE}

mkdir -p ${DST}
aws --no-sign-request s3 sync ${SRC} ${DST}/

echo "files: $(find ${DST} -type f | wc -l)  (the S3 listing has 338 objects under GRCh38)"
echo "size:  $(du -sh ${DST} | cut -f1)  (the S3 listing says 130.8 GiB)"
