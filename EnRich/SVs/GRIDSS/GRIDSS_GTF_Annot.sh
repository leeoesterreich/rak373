#!/bin/bash
#SBATCH --job-name=GRIDSS_GTF_Annot
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH --output=slurm_logs/gridss_gtf_annot_%j.out
#SBATCH --error=slurm_logs/gridss_gtf_annot_%j.err

set -euo pipefail

# -----------------------------
# Modules
# -----------------------------
module purge
module load bedtools/2.31.1
module load bcftools/1.22
module load htslib/1.22

# -----------------------------
# Paths
# -----------------------------
BASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/GRIDSS_012425"
OUTBASE="${BASE}/GRIDSS_ANNOT_GTF"
mkdir -p "$OUTBASE" slurm_logs

REFDIR="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human"
GENE_BED="${REFDIR}/gencode.v19.genes.bed"

SAMPLE_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet_012425.csv"

# -----------------------------
# Checks
# -----------------------------
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
[[ -s "$GENE_BED" ]] || { echo "ERROR: gene BED not found: $GENE_BED"; exit 1; }
[[ -s "$SAMPLE_SHEET" ]] || { echo "ERROR: sample sheet not found: $SAMPLE_SHEET"; exit 1; }

# -----------------------------
# Read samples from CSV (first column)
# -----------------------------
mapfile -t SAMPLES < <(cut -d',' -f1 "$SAMPLE_SHEET" | tail -n +2)

echo "[INFO] Samples from sample sheet:"
printf '  %s\n' "${SAMPLES[@]}"

# -----------------------------
# Per-sample processing (BND only)
# -----------------------------
for sample in "${SAMPLES[@]}"; do
  sample="$(echo "$sample" | xargs)"
  [[ -n "$sample" ]] || continue

  echo "------------------------------------------------------------"
  echo "[INFO] Sample: $sample"

  sample_dir="${BASE}/${sample}"
  vcf="${sample_dir}/${sample}.gridss.filtered.vcf.gz"

  if [[ ! -s "$vcf" ]]; then
    echo "[WARN] Missing filtered GRIDSS VCF: $vcf  (skipping)"
    continue
  fi

  outdir="${OUTBASE}/${sample}"
  mkdir -p "$outdir"

  bed="${outdir}/${sample}.BND.breakends.bed"
  withgenes="${outdir}/${sample}.BND.breakends.withGenes.tsv"
  pairs="${outdir}/${sample}.BND.pairs.tsv"

  # Ensure index
  if [[ ! -s "${vcf}.tbi" ]]; then
    echo "[INFO] Indexing VCF: $vcf"
    tabix -f -p vcf "$vcf" || true
  fi

  # ------------------------------------------------------------
  # Extract breakends:
  #  - Side A: CHROM:POS from the VCF record
  #  - Side B: chr:pos parsed from ALT breakend notation
  #    ALT examples:
  #      N]chr6:152316374]
  #      ]chr6:152316374]N
  #      N[chr6:152316374[
  #      [chr6:152316374[N
  # ------------------------------------------------------------
  bcftools query -i 'INFO/SVTYPE="BND"' \
    -f '%CHROM\t%POS\t%ID\t%ALT\n' "$vcf" \
  | awk 'BEGIN{OFS="\t"}
  function emit(chr, pos, id, side,   start, end){
    if(chr=="" || pos=="" || pos=="." ) return
    start = pos-1; end = pos
    if(start<0) start=0
    print chr, start, end, id, side
  }
  {
    chrA=$1; posA=$2; id=$3; alt=$4
    emit(chrA, posA, id, "A")

    # Extract "something:number" from ALT by stripping brackets and bases
    tmp=alt
    gsub(/\[/, " ", tmp); gsub(/\]/, " ", tmp)
    # now tmp contains tokens; find token with colon and numeric pos
    n=split(tmp, a, /[ \t]/)
    for(i=1;i<=n;i++){
      if(a[i] ~ /:/){
        split(a[i], b, ":")
        if(b[2] ~ /^[0-9]+$/){
          chrB=b[1]; posB=b[2]
          emit(chrB, posB, id, "B")
          break
        }
      }
    }
  }' > "$bed"

  if [[ ! -s "$bed" ]]; then
    echo "[WARN] No BND breakends extracted for $sample (empty bed)."
    continue
  fi

  # Intersect with genes
  bedtools intersect -wa -wb \
    -a "$bed" \
    -b "$GENE_BED" \
  | awk 'BEGIN{OFS="\t"}{print $4,$5,$1,$3,$9,$10}' \
  > "$withgenes"

  # Build fusion-style pair table
  awk -F'\t' 'BEGIN{OFS="\t"}
    {
      id=$1; side=$2; chr=$3; pos=$4; gene=$5
      if(side=="A"){
        chrA[id]=chr; posA[id]=pos
        geneA[id]=geneA[id] (geneA[id]?";":"") gene
      } else if(side=="B"){
        chrB[id]=chr; posB[id]=pos
        geneB[id]=geneB[id] (geneB[id]?";":"") gene
      }
    }
    END{
      print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SV_ID","Sample"
      for(id in chrA){
        print chrA[id],posA[id],(geneA[id]?geneA[id]:"NA"), \
              (chrB[id]?chrB[id]:"NA"),(posB[id]?posB[id]:"NA"),(geneB[id]?geneB[id]:"NA"), \
              id,"'"$sample"'"
      }
    }' "$withgenes" > "$pairs"

  echo "[INFO] Wrote: $pairs"
done

# -----------------------------
# Cohort-level table
# -----------------------------
Summary="${OUTBASE}/ALL_SAMPLES.BND.pairs.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
FNR==1{next}  # skip header in each per-sample pairs file
{
  key = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7
  s = $8

  if(!((key SUBSEP s) in seen)){
    seen[key SUBSEP s] = 1
    nsamp[key]++
    slist[key] = (slist[key] ? slist[key] "," s : s)
  }

  chrA[key]=$1; posA[key]=$2; geneA[key]=$3
  chrB[key]=$4; posB[key]=$5; geneB[key]=$6
  sv[key]=$7
}
END{
  print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SV_ID","Unique_Samples","Sample_List"
  for(k in nsamp){
    print chrA[k],posA[k],geneA[k],chrB[k],posB[k],geneB[k],sv[k],nsamp[k],slist[k]
  }
}' "${OUTBASE}"/*/*.BND.pairs.tsv > "$Summary" 2>/dev/null || true

echo "============================================================"
echo "[DONE]"
echo "Summary table: $Summary"
