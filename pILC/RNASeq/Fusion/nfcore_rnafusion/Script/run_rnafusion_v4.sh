#!/bin/bash
#SBATCH --job-name=rnafusion_v4
#SBATCH --clusters=htc
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/nfcore_rnafusion/logs/fusion_v4_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/nfcore_rnafusion/logs/fusion_v4_%j.err
#SBATCH --cpus-per-task=4            # head job coordinates and submits jobs
#SBATCH --mem=8G
#SBATCH --time=72:00:00              # if the run outlives this, resubmit the same script (-resume continues)

# Usage: sbatch run_rnafusion_v4.sh                       (defaults below)
#        sbatch run_rnafusion_v4.sh <SAMPLESHEET> <OUTDIR> (optional override)

BASE=/ix1/alee/LO_LAB/Personal/Rahul/pILC/RNASeq/RNAfusion/nfcore_rnafusion
VAST4=/vast/alee/Rahul/PDO/RNAFusion_4.1.3
SHEET=${1:-${BASE}/samplesheet_fusion.csv}
OUTDIR=${2:-${VAST4}/results}
REFS=/ix1/alee/LO_LAB/General/References/rnafusion_4.1.3    # --genomes_base (contains GRCh38/)
CONFIG=${BASE}/rnafusion_v4.config
CRC_CONFIG=/software/rhel9/manual/install/nf-core/pipelines/config/htc.config
VERSION=4.1.3

module load nextflow/25.10.3 singularity/4.3.2
export NXF_VER=25.10.4                                      # 4.1.3 needs Nextflow >= 25.10.4
export NXF_HOME=${VAST4}/.nextflow
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_SINGULARITY_CACHEDIR=${VAST4}/singularity_cache
export SINGULARITY_CACHEDIR=${VAST4}/singularity_oci_cache
export SINGULARITY_TMPDIR=${VAST4}/singularity_tmp
mkdir -p ${VAST4} ${NXF_SINGULARITY_CACHEDIR} ${SINGULARITY_CACHEDIR} ${SINGULARITY_TMPDIR} ${BASE}/logs
cd ${BASE}

[[ -f ${CONFIG} ]] || { echo "ERROR: config not found: ${CONFIG}" >&2; exit 1; }
[[ -f ${SHEET} ]]  || { echo "ERROR: samplesheet not found: ${SHEET}" >&2; exit 1; }
[[ -f ${REFS}/GRCh38/fusion_report_db/mitelman.db ]] || { echo "ERROR: references incomplete under ${REFS}" >&2; exit 1; }

echo "Samplesheet: ${SHEET} | Outdir: ${OUTDIR} | Version: ${VERSION} | Nextflow: ${NXF_VER} | Refs: ${REFS}"

nextflow run nf-core/rnafusion -r ${VERSION} -profile singularity \
    -c ${CRC_CONFIG} -c ${CONFIG} -work-dir ${VAST4}/work \
    --input ${SHEET} --outdir ${OUTDIR} \
    --genomes_base ${REFS} \
    --tools all --no_cosmic --publish_dir_mode link \
    -resume
