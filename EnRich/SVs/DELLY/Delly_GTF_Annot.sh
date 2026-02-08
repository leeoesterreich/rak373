#!/bin/bash
#SBATCH --job-name=Delly_GTF_Annot
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --output=slurm_log/delly_gtf_annot_%j.out
#SBATCH --error=slurm_log/delly_gtf_annot_%j.err

set -euo pipefail

# -----------------------------
# Paths
# -----------------------------
BASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Delly_010526"
OUTBASE="${BASE}/Delly_ANNOT_GTF"
mkdir -p "$OUTBASE"

REFDIR="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human"
GENE_BED="${REFDIR}/gencode.v19.genes.bed"

SAMPLE_SHEET="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/sample_sheet_010526.csv"

# -----------------------------
# Modules
# -----------------------------
module purge
module load bedtools/2.31.1
module load bcftools/1.22
module load htslib/1.22

# -----------------------------
# Checks
# -----------------------------
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }
[[ -s "$GENE_BED" ]] || { echo "ERROR: gene BED not found: $GENE_BED"; exit 1; }
[[ -s "$SAMPLE_SHEET" ]] || { echo "ERROR: sample sheet not found: $SAMPLE_SHEET"; exit 1; }

# -----------------------------
# Read samples from CSV (col1, skip header)
# -----------------------------
mapfile -t SAMPLES < <(cut -d',' -f1 "$SAMPLE_SHEET" | tail -n +2)

echo "[INFO] Samples from sample sheet:"
printf '  %s\n' "${SAMPLES[@]}"

# -----------------------------
# Per-sample processing (BND only)
# -----------------------------
for sample in "${SAMPLES[@]}"; do
  echo "------------------------------------------------------------"
  echo "[INFO] Sample: $sample"

  vcf="${BASE}/${sample}/${sample}_delly_filtered_merged.vcf.gz"
  if [[ ! -s "$vcf" ]]; then
    echo "[WARN] Missing VCF: $vcf  (skipping)"
    continue
  fi

  outdir="${OUTBASE}/${sample}"
  mkdir -p "$outdir"

  bed="${outdir}/${sample}.BND.breakends.bed"
  withgenes="${outdir}/${sample}.BND.breakends.withGenes.tsv"
  pairs="${outdir}/${sample}.BND.pairs.tsv"
  support="${outdir}/${sample}.BND.support.tsv"

  # ------------------------------------------------------------
  # Support map: SV_ID -> SR, PE (from filtered merged VCF)
  # ------------------------------------------------------------
  bcftools query -i 'INFO/SVTYPE=="BND"' \
    -f '%ID\t%INFO/SR\t%INFO/PE\n' \
    "$vcf" \
  | awk 'BEGIN{OFS="\t"}{
      sr=($2=="."||$2=="")?0:$2;
      pe=($3=="."||$3=="")?0:$3;
      print $1,sr,pe
    }' \
  | sort -u > "$support"

  # ------------------------------------------------------------
  # Extract chr/pos for both breakends (A and B) as 1-based positions
  # ------------------------------------------------------------
  bcftools query -i 'INFO/SVTYPE=="BND"' \
    -f '%CHROM\t%POS\t%ID\tA\n%INFO/CHR2\t%INFO/POS2\t%ID\tB\n' \
    "$vcf" \
  | awk 'BEGIN{OFS="\t"}
    {
      chr=$1; pos=$2; id=$3; side=$4
      if(chr=="" || pos=="" || pos=="." ) next
      start=pos-1; end=pos
      if(start<0) start=0
      print chr,start,end,id,side
    }' > "$bed"

  # ------------------------------------------------------------
  # Intersect breakends with genes
  # NOTE: This keeps only overlaps (same as your original).
  # If you want to retain non-overlaps too, switch to: bedtools intersect -loj
  # ------------------------------------------------------------
  bedtools intersect -wa -wb \
    -a "$bed" \
    -b "$GENE_BED" \
  | awk 'BEGIN{OFS="\t"}{print $4,$5,$1,$3,$9,$10}' \
  > "$withgenes"

  # ------------------------------------------------------------
  # Build fusion-style pair table + add SR/PE/TOTAL_SUPPORT
  # ------------------------------------------------------------
  awk -F'\t' 'BEGIN{OFS="\t"}
    # File 1: support map
    FNR==NR{
      id=$1; sr[id]=$2; pe[id]=$3
      next
    }
    # File 2: withgenes
    {
      id=$1; side=$2
      if(side=="A"){
        chrA[id]=$3; posA[id]=$4
        geneA[id]=geneA[id] (geneA[id]?";":"") $5
      } else if(side=="B"){
        chrB[id]=$3; posB[id]=$4
        geneB[id]=geneB[id] (geneB[id]?";":"") $5
      }
    }
    END{
      print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SV_ID","Sample","SR","PE","TOTAL_SUPPORT"
      for(id in chrA){
        s=(id in sr ? sr[id] : 0)
        p=(id in pe ? pe[id] : 0)
        gA=(geneA[id]?geneA[id]:"NA")
        gB=(geneB[id]?geneB[id]:"NA")
        print chrA[id],posA[id],gA,chrB[id],posB[id],gB,id,"'"$sample"'",s,p,(s+p)
      }
    }' "$support" "$withgenes" > "$pairs"

  echo "[INFO] Wrote: $pairs"
done

# -----------------------------
# Cohort-level table (aggregated by breakpoint+genes+SV_ID)
# Adds:
#   - Max_TOTAL_SUPPORT across samples
#   - Support_BySample listing SR/PE/TOT per sample
# -----------------------------
Summary="${OUTBASE}/ALL_SAMPLES.BND.pairs.tsv"

awk -F'\t' 'BEGIN{OFS="\t"}
FNR==1{next}  # skip header in each per-sample pairs file
{
  # Per-sample columns:
  # 1 ChrA, 2 PositionA, 3 Gene_ChrA, 4 ChrB, 5 PositionB, 6 Gene_ChrB, 7 SV_ID, 8 Sample, 9 SR, 10 PE, 11 TOTAL_SUPPORT
  key = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7

  s=$8
  sr=($9==""?0:$9)
  pe=($10==""?0:$10)
  tot=($11==""?0:$11)

  # count unique samples per key
  if(!((key SUBSEP s) in seen)){
    seen[key SUBSEP s]=1
    nsamp[key]++
    slist[key]=(slist[key]?slist[key] "," s:s)
  }

  # store fields once
  chrA[key]=$1; posA[key]=$2; geneA[key]=$3
  chrB[key]=$4; posB[key]=$5; geneB[key]=$6
  sv[key]=$7

  # max total support across samples
  if(!(key in maxTot) || tot>maxTot[key]) maxTot[key]=tot

  # per-sample support string
  if(!((key SUBSEP s) in supSeen)){
    supSeen[key SUBSEP s]=1
    supList[key]=(supList[key]?supList[key] ",":"") s ":SR=" sr ";PE=" pe ";TOT=" tot
  }
}
END{
  print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SV_ID","Unique_Samples","Sample_List","Max_TOTAL_SUPPORT","Support_BySample"
  for(k in nsamp){
    print chrA[k],posA[k],geneA[k],chrB[k],posB[k],geneB[k],sv[k],nsamp[k],slist[k],maxTot[k],supList[k]
  }
}' "${OUTBASE}"/*/*.BND.pairs.tsv > "$Summary"

echo "============================================================"
echo "[DONE]"
echo "Summary table: $Summary"
