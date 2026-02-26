#!/bin/bash
#SBATCH --job-name=Delly_BND_GTF_Annot
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
BASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/Delly/Delly_010526"
OUTBASE="${BASE}/Delly_ANNOT_GTF"
mkdir -p "$OUTBASE" slurm_log

REFDIR="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human"
GENE_BED="${REFDIR}/gencode.v19.genes.bed"

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

# -----------------------------
# Discover sample VCFs automatically
# -----------------------------
mapfile -t VCFS < <(find "$BASE" -type f -name '*_delly_filtered_merged.vcf.gz' | sort)

if [[ ${#VCFS[@]} -eq 0 ]]; then
  echo "ERROR: No *_delly_filtered_merged.vcf.gz found under: $BASE"
  exit 1
fi

echo "[INFO] Found ${#VCFS[@]} DELLY merged VCFs"

# -----------------------------
# Helper: normalize chrom names to "chr*"
# -----------------------------
norm_chr_stream() {
  awk 'BEGIN{OFS="\t"}{
    c=$1
    if(c !~ /^chr/) c="chr" c
    $1=c
    print
  }'
}

# -----------------------------
# Per-VCF processing (BND only)
# Writes ONLY: ${sample}.BND.pairs.tsv
# -----------------------------
for vcf in "${VCFS[@]}"; do
  # sample is filename without suffix
  # e.g. TP20_M388_delly_filtered_merged.vcf.gz -> TP20_M388
  sample="$(basename "$vcf" _delly_filtered_merged.vcf.gz)"

  echo "------------------------------------------------------------"
  echo "[INFO] Sample: $sample"
  echo "[INFO] VCF: $vcf"

  outdir="${OUTBASE}/${sample}"
  mkdir -p "$outdir"

  pairs="${outdir}/${sample}.BND.pairs.tsv"

  awk -F'\t' -v OFS='\t' -v SAMPLE="$sample" '
    # 1) support stream
    NR==FNR{
      id=$1
      sr[id]=($2=="."||$2==""?0:$2)
      pe[id]=($3=="."||$3==""?0:$3)
      next
    }

    # 2) breakends+genes stream (bedtools -loj)
    {
      id=$4
      side=$5
      chr=$1
      pos=$3

      gname=$9
      if(gname=="." || gname=="") gname="NA"

      if(side=="A"){
        chrA[id]=chr; posA[id]=pos
        if(gname!="NA") geneA[id]=geneA[id] (geneA[id]?";":"") gname
      } else if(side=="B"){
        chrB[id]=chr; posB[id]=pos
        if(gname!="NA") geneB[id]=geneB[id] (geneB[id]?";":"") gname
      }
    }

    END{
      print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SV_ID","Sample","SR","PE","TOTAL_SUPPORT"

      for(id in chrA){
        s=(id in sr ? sr[id] : 0)
        p=(id in pe ? pe[id] : 0)
        gA=(geneA[id] && geneA[id]!="") ? geneA[id] : "NA"
        gB=(geneB[id] && geneB[id]!="") ? geneB[id] : "NA"
        cB=(id in chrB) ? chrB[id] : "NA"
        pB=(id in posB) ? posB[id] : "NA"
        print chrA[id],posA[id],gA,cB,pB,gB,id,SAMPLE,s,p,(s+p)
      }

      for(id in chrB){
        if(id in chrA) continue
        s=(id in sr ? sr[id] : 0)
        p=(id in pe ? pe[id] : 0)
        gB=(geneB[id] && geneB[id]!="") ? geneB[id] : "NA"
        print "NA","NA","NA",chrB[id],posB[id],gB,id,SAMPLE,s,p,(s+p)
      }
    }
  ' \
  <(
    bcftools query -i 'INFO/SVTYPE=="BND"' -f '%ID\t%INFO/SR\t%INFO/PE\n' "$vcf" \
    | awk 'BEGIN{OFS="\t"}{
        sr=($2=="."||$2=="")?0:$2
        pe=($3=="."||$3=="")?0:$3
        print $1,sr,pe
      }' | sort -u
  ) \
  <(
    bcftools query -i 'INFO/SVTYPE=="BND"' \
      -f '%CHROM\t%POS\t%ID\tA\n%INFO/CHR2\t%INFO/POS2\t%ID\tB\n' "$vcf" \
    | awk 'BEGIN{OFS="\t"}
      {
        chr=$1; pos=$2; id=$3; side=$4
        if(chr=="" || pos=="" || pos==".") next
        start=pos-1; end=pos
        if(start<0) start=0
        print chr,start,end,id,side
      }' \
    | norm_chr_stream \
    | bedtools intersect -loj -a - -b "$GENE_BED"
  ) > "$pairs"

  echo "[INFO] Wrote: $pairs"
done

# -----------------------------
# Cohort-level summary
# -----------------------------
Summary="${OUTBASE}/ALL_SAMPLES.BND.pairs.tsv"

mapfile -t PAIRS_FILES < <(find "$OUTBASE" -type f -name "*.BND.pairs.tsv" ! -name "ALL_SAMPLES.*" | sort)

if [[ ${#PAIRS_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No per-sample *.BND.pairs.tsv found under $OUTBASE"
  exit 1
fi

awk -F'\t' 'BEGIN{OFS="\t"}
FNR==1{next}
{
  key = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7

  s=$8
  sr=($9==""?0:$9)
  pe=($10==""?0:$10)
  tot=($11==""?0:$11)

  if(!((key SUBSEP s) in seen)){
    seen[key SUBSEP s]=1
    nsamp[key]++
    slist[key]=(slist[key]?slist[key] "," s:s)
  }

  chrA[key]=$1; posA[key]=$2; geneA[key]=$3
  chrB[key]=$4; posB[key]=$5; geneB[key]=$6
  sv[key]=$7

  if(!(key in maxTot) || tot>maxTot[key]) maxTot[key]=tot

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
}' "${PAIRS_FILES[@]}" > "$Summary"

echo "[DONE] Summary table: $Summary"
