#!/bin/bash
#SBATCH --job-name=clear_vast
#SBATCH --clusters=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/clear_vast_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs/clear_vast_%j.err
# Deletes work/ and results/ on /vast after batch has been copied AND verified.
# Keeps singularity_cache, singularity_oci_cache and .nextflow.
V=/vast/alee/Rahul/PDO/RNAFusion_4.1.3
LOGS=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/logs
DEST=/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion/20240315_Organoid_RNAseq/results

grep -q "VERIFIED: copy is complete" $(ls -t ${LOGS}/copy_res_*.out | head -n 1) 2>/dev/null \
    || { echo "ERROR: no VERIFIED line in the batch copy log. Nothing deleted." >&2; exit 1; }
N=$(ls -d ${DEST}/*/ 2>/dev/null | wc -l)
[[ ${N} -ge 19 ]] || { echo "ERROR: only ${N} folders in the copy (expected 19). Nothing deleted." >&2; exit 1; }

for d in work results; do
    if [[ -d ${V}/${d} ]]; then
        echo "Deleting ${V}/${d}"
        find ${V}/${d} -mindepth 1 -maxdepth 1 -print0 | xargs -0 -P 8 -n 1 rm -rf
        rm -rf ${V}/${d}
    fi
done
echo "Done. Remaining in ${V}:"; ls -la ${V}
