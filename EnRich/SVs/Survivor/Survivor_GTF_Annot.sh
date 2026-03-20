#!/bin/bash
#SBATCH --job-name=Survivor_SV_GTF_Annot
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --output=slurm_log/survivor_gtf_annot_%j.out
#SBATCH --error=slurm_log/survivor_gtf_annot_%j.err

set -euo pipefail

# -----------------------------
# Paths (UPDATED for your Batch2_021326)
# -----------------------------
BASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20250124/Survivor_012425/merged_bnd/"
OUTBASE="/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/EnRich/20250124/Survivor_012425/merged_bnd_GTF_ANNOT"
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
command -v bgzip   >/dev/null 2>&1 || { echo "ERROR: bgzip not found"; exit 1; }
command -v tabix   >/dev/null 2>&1 || { echo "ERROR: tabix not found"; exit 1; }
[[ -s "$GENE_BED" ]] || { echo "ERROR: gene BED not found: $GENE_BED"; exit 1; }
[[ -d "$BASE" ]] || { echo "ERROR: BASE dir not found: $BASE"; exit 1; }

# -----------------------------
# Discover VCFs
# -----------------------------
mapfile -t VCFS < <(find "$BASE" -maxdepth 1 -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) | sort)
if [[ ${#VCFS[@]} -eq 0 ]]; then
  echo "ERROR: No VCFs found under: $BASE"
  exit 1
fi

echo "[INFO] Found ${#VCFS[@]} VCFs"
printf '[INFO]  - %s\n' "${VCFS[@]##*/}"

TMPVCF="${OUTBASE}/tmp_vcfs"
mkdir -p "$TMPVCF"

ensure_bgzip() {
  local in="$1"
  if [[ "$in" == *.vcf.gz ]]; then
    echo "$in"
    return 0
  fi
  local bn out
  bn="$(basename "$in")"
  out="${TMPVCF}/${bn}.gz"
  [[ -s "$out" ]] || { echo "[INFO] bgzip: $in -> $out" >&2; bgzip -c "$in" > "$out"; }
  [[ -s "${out}.tbi" ]] || tabix -f -p vcf "$out"
  echo "$out"
}

# -----------------------------
# Per-sample annotation
# -----------------------------
for vcf0 in "${VCFS[@]}"; do
  vcf="$(ensure_bgzip "$vcf0")"

  b="$(basename "$vcf0")"
  # TP19_M876_merged.BND.vcf -> TP19_M876
  sample="${b%%_merged*}"
  sample="${sample%%.vcf*}"
  [[ -n "$sample" ]] || sample="SAMPLE"

  echo "------------------------------------------------------------"
  echo "[INFO] Sample: $sample"
  echo "[INFO] VCF: $vcf0"

  outdir="${OUTBASE}/${sample}"
  mkdir -p "$outdir"

  out_tsv="${outdir}/${sample}.SV.pairs.tsv"
  tmp_out="${out_tsv}.tmp"

  # Do not let one bad file kill the whole job
  set +e
  (
    awk -F'\t' -v OFS='\t' -v SAMPLE="$sample" '
      NR==FNR{
        id=$1
        svtype[id]=$2
        svlen[id]=$3
        supp[id]=$4
        vec[id]=$5
        next
      }
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
        print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SV_ID","Sample","SVTYPE","SVLEN","SUPP","SUPP_VEC"
        for(id in chrA){
          gA=(geneA[id] && geneA[id]!="") ? geneA[id] : "NA"
          gB=(geneB[id] && geneB[id]!="") ? geneB[id] : "NA"
          cB=(id in chrB) ? chrB[id] : "NA"
          pB=(id in posB) ? posB[id] : "NA"
          t =(id in svtype) ? svtype[id] : "NA"
          l =(id in svlen)  ? svlen[id]  : "NA"
          sp=(id in supp)   ? supp[id]   : "NA"
          v =(id in vec)    ? vec[id]    : "NA"
          print chrA[id],posA[id],gA,cB,pB,gB,id,SAMPLE,t,l,sp,v
        }
      }
    ' \
    <(
      bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC\t%INFO/CHR2\t%INFO/END\n' "$vcf" \
      | awk 'BEGIN{OFS="\t"}
        {
          chrom=$1; pos=$2; id=$3; t=$4; l=$5; s=$6; v=$7; c2=$8; e2=$9
          if(id=="" || id=="."){
            if(c2=="" || c2==".") c2="NA"
            if(e2=="" || e2==".") e2="NA"
            if(t==""  || t==".")  t="NA"
            id=chrom ":" pos ":" c2 ":" e2 ":" t
          }
          t=(t==""||t==".")?"NA":t
          l=(l==""||l==".")?"NA":l
          s=(s==""||s==".")?"NA":s
          v=(v==""||v==".")?"NA":v
          print id,t,l,s,v
        }' | sort -u
    ) \
    <(
      bcftools query -f '%CHROM\t%POS\t%ID\t%ALT\t%INFO/CHR2\t%INFO/END\t%INFO/SVTYPE\n' "$vcf" \
      | awk 'BEGIN{OFS="\t"}
        function is_int(x){ return (x ~ /^[0-9]+$/) }
        function norm_chr(c){ if(c!="" && c!="." && c !~ /^chr/) c="chr" c; return c }
        function is_std_chr(c){ c=norm_chr(c); return (c ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/) }

        # Extract partner chr:pos from breakend ALT (handles N]chr:pos], ]chr:pos]N, N[chr:pos[, etc.)
        function alt_partner(alt,   hit){
          hit=alt
          gsub(/[\[\]]/, "", hit)
          gsub(/[ACGTN]/, "", hit)
          if (match(hit, /(chr)?([0-9]+|X|Y):[0-9]+/)) {
            hit=substr(hit, RSTART, RLENGTH)
            return hit
          }
          return ""
        }

        {
          chrA=$1; posA=$2; id=$3; alt=$4; c2=$5; e2=$6; t=$7

          if(id=="" || id=="."){
            if(c2=="" || c2==".") c2="NA"
            if(e2=="" || e2==".") e2="NA"
            if(t==""  || t==".")  t="NA"
            id=chrA ":" posA ":" c2 ":" e2 ":" t
          }

          if(!is_int(posA) || !is_std_chr(chrA)) next

          chrB=""; posB=""
          # Prefer CHR2/END if END is numeric
          if(c2!="" && c2!="." && e2!="" && e2!="." && is_int(e2)){
            chrB=c2; posB=e2
          } else {
            p = alt_partner(alt)
            if(p!=""){
              split(p, a, ":")
              chrB=a[1]; posB=a[2]
            }
          }

          if(chrB=="" || posB=="" || !is_int(posB) || !is_std_chr(chrB)) next

          chrA=norm_chr(chrA)
          chrB=norm_chr(chrB)

          sA=posA-1; eA=posA; if(sA<0) sA=0
          sB=posB-1; eB=posB; if(sB<0) sB=0

          print chrA, sA, eA, id, "A"
          print chrB, sB, eB, id, "B"
        }' \
      | bedtools intersect -loj -a - -b "$GENE_BED"
    ) > "$tmp_out"
  )
  rc=$?
  set -e

  if [[ $rc -ne 0 ]]; then
    echo "[ERROR] Failed sample: $sample (VCF: $vcf0). Check slurm_log/*.err"
    rm -f "$tmp_out"
    continue
  fi

  mv -f "$tmp_out" "$out_tsv"
  echo "[INFO] Wrote: $out_tsv (rows excluding header: $(awk 'END{print (NR>0?NR-1:0)}' "$out_tsv"))"
done

# -----------------------------
# Cohort-level summary
# -----------------------------
SUMMARY="${OUTBASE}/ALL_SAMPLES.SV.pairs.tsv"
mapfile -t TSVS < <(find "$OUTBASE" -type f -name "*.SV.pairs.tsv" ! -name "ALL_SAMPLES.*" | sort)

if [[ ${#TSVS[@]} -eq 0 ]]; then
  echo "ERROR: No per-sample *.SV.pairs.tsv found under $OUTBASE"
  exit 1
fi

data_lines=$(awk 'FNR==1{next} {c++} END{print c+0}' "${TSVS[@]}")
if [[ "$data_lines" -eq 0 ]]; then
  echo "WARNING: No SV rows found in any per-sample TSV (header-only). Writing cohort header only."
  echo -e "ChrA\tPositionA\tGene_ChrA\tChrB\tPositionB\tGene_ChrB\tSVTYPE\tUnique_Samples\tSample_List\tMax_SUPP\tSUPP_VEC_bySample" > "$SUMMARY"
  echo "[DONE] Cohort summary (header only): $SUMMARY"
  exit 0
fi

awk -F'\t' 'BEGIN{OFS="\t"}
FNR==1{next}
{
  key = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $9
  samp=$8
  supp=$11
  vec=$12

  if(!((key SUBSEP samp) in seen)){
    seen[key SUBSEP samp]=1
    nsamp[key]++
    slist[key]=(slist[key]?slist[key] "," samp:samp)
  }

  chrA[key]=$1; posA[key]=$2; geneA[key]=$3
  chrB[key]=$4; posB[key]=$5; geneB[key]=$6
  type[key]=$9

  if(supp ~ /^[0-9]+$/){
    if(!(key in maxSupp) || supp>maxSupp[key]) maxSupp[key]=supp
  } else {
    if(!(key in maxSupp)) maxSupp[key]="NA"
  }

  vecset[key]= (vecset[key]? vecset[key] ";" : "") samp ":" vec
}
END{
  print "ChrA","PositionA","Gene_ChrA","ChrB","PositionB","Gene_ChrB","SVTYPE","Unique_Samples","Sample_List","Max_SUPP","SUPP_VEC_bySample"
  for(k in nsamp){
    print chrA[k],posA[k],geneA[k],chrB[k],posB[k],geneB[k],type[k],nsamp[k],slist[k],(k in maxSupp?maxSupp[k]:"NA"),vecset[k]
  }
}' "${TSVS[@]}" > "$SUMMARY"

echo "[DONE] Cohort summary: $SUMMARY"