#!/usr/bin/env bash

# Required to run BreakID

set -euo pipefail

# ============================================================
# Prepare hg19 .nib files for BreakID (ALL contigs)
# - Downloads UCSC hg19 FASTA contigs (chr*.fa)
# - Builds ONLY what is needed from blatSrc35: jkweb.a + faToNib
# - Converts ALL chr*.fa (main + random + hap + chrUn) -> .nib
#
# Output:  $NIB_DIR/<contig>.nib  (for every <contig>.fa in $FA_DIR)
# ============================================================

BASE="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human"
FA_DIR="${BASE}/hg19_fa"
NIB_DIR="${BASE}/hg19_nib_files"
SRC_DIR="${BASE}/hg19_kent_src"          # will contain blatSrc/
FA2NIB_BIN="${HOME}/bin/faToNib"

UCSC_RSYNC="rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
BLAT_ZIP_URL="https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip"

# ----------------------------
# 0) Setup
# ----------------------------
mkdir -p "$FA_DIR" "$NIB_DIR" "$SRC_DIR" "$(dirname "$FA2NIB_BIN")"

module purge >/dev/null 2>&1 || true
module load gcc >/dev/null 2>&1 || true

# ----------------------------
# 1) Download hg19 chr*.fa (only if missing)
# ----------------------------
if [[ ! -s "$FA_DIR/chr1.fa" ]]; then
  echo "[INFO] Downloading hg19 chr*.fa(.gz) to $FA_DIR ..."
  cd "$FA_DIR"
  rsync -avzP \
    --include="chr*.fa.gz" \
    --include="chr*.fa" \
    --exclude="*" \
    "$UCSC_RSYNC" .

  # unzip if needed
  gunzip -f *.fa.gz 2>/dev/null || true
fi

[[ -s "$FA_DIR/chr1.fa" ]] || { echo "[ERROR] Missing $FA_DIR/chr1.fa after download" >&2; exit 1; }

# ----------------------------
# 2) Build faToNib (only what is needed)
# ----------------------------
if [[ ! -x "$FA2NIB_BIN" ]]; then
  echo "[INFO] Building faToNib from blatSrc35 into $FA2NIB_BIN ..."

  cd "$SRC_DIR"
  if [[ ! -d "$SRC_DIR/blatSrc" ]]; then
    wget -O blatSrc35.zip "$BLAT_ZIP_URL"
    unzip -q blatSrc35.zip
  fi
  [[ -d "$SRC_DIR/blatSrc" ]] || { echo "[ERROR] blatSrc directory not found under $SRC_DIR" >&2; exit 1; }

  # Build jkweb.a (library needed by faToNib)
  cd "$SRC_DIR/blatSrc/lib"
  make -j 8

  # Ensure jkweb.a is available at lib/jkweb.a (faToNib link expects this)
  cd "$SRC_DIR/blatSrc"
  if [[ -f "lib/x86_64/jkweb.a" ]]; then
    ln -sf x86_64/jkweb.a lib/jkweb.a
  fi
  [[ -f "lib/jkweb.a" ]] || { echo "[ERROR] jkweb.a not found at $SRC_DIR/blatSrc/lib/jkweb.a" >&2; exit 1; }

  # Build faToNib only (do NOT build whole kent tree)
  cd "$SRC_DIR/blatSrc/utils/faToNib"
  make clean >/dev/null 2>&1 || true
  make -j 8

  # Install faToNib
  [[ -x "$SRC_DIR/blatSrc/utils/faToNib/faToNib" ]] || { echo "[ERROR] faToNib build failed" >&2; exit 1; }
  cp -f "$SRC_DIR/blatSrc/utils/faToNib/faToNib" "$FA2NIB_BIN"
fi

[[ -x "$FA2NIB_BIN" ]] || { echo "[ERROR] faToNib not executable at $FA2NIB_BIN" >&2; exit 1; }

# ----------------------------
# 3) Convert FASTA -> NIB (ALL contigs)
# ----------------------------
echo "[INFO] Converting ALL chr*.fa -> chr*.nib into $NIB_DIR ..."
shopt -s nullglob
fa_list=( "$FA_DIR"/*.fa )
[[ ${#fa_list[@]} -gt 0 ]] || { echo "[ERROR] No .fa files found in $FA_DIR" >&2; exit 1; }

n_total=0
n_new=0

for fa in "${fa_list[@]}"; do
  base=$(basename "$fa" .fa)
  out="${NIB_DIR}/${base}.nib"
  ((n_total+=1))

  # skip if already created
  if [[ -s "$out" ]]; then
    continue
  fi

  echo "[INFO] Converting $base"
  "$FA2NIB_BIN" "$fa" "$out"
  ((n_new+=1))
done

# ----------------------------
# 4) Verify
# ----------------------------
n=$(ls -1 "$NIB_DIR"/*.nib 2>/dev/null | wc -l | tr -d ' ')
echo "[INFO] FASTA contigs seen: $n_total"
echo "[INFO] New NIBs created:  $n_new"
echo "[INFO] Total NIBs now:    $n"

if [[ "$n" -lt 25 ]]; then
  echo "[ERROR] NIB generation looks incomplete (found $n .nib files in $NIB_DIR)" >&2
  exit 1
fi

echo "[INFO] DONE. Example outputs:"
ls -1 "$NIB_DIR"/*.nib | head
