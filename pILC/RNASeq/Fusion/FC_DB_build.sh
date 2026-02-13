#!/usr/bin/env bash
#SBATCH --job-name=FusionCatcher_db
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1-00:00
#SBATCH --output=slurm_log/FC_DB_download_%j.out
#SBATCH --error=slurm_log/FC_DB_download_%j.err

set -euo pipefail

DBROOT="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/FusionCatcher_db"
BUILD="human_v102"
URL_BASE="https://sourceforge.net/projects/fusioncatcher/files/data"
PARTS=(aa ab ac ad)

# ----------------------------
# Prep
# ----------------------------
mkdir -p slurm_log
mkdir -p "$DBROOT"
cd "$DBROOT"

echo "[$(date)] Starting FusionCatcher DB download"
echo "DBROOT=$DBROOT"
echo "BUILD=$BUILD"
echo "Node: $(hostname)"
echo

# Helpful: quick connectivity check (won't fail the job if blocked)
echo "[$(date)] Connectivity check..."
( curl -Is "${URL_BASE}/${BUILD}.tar.gz.aa/download" | head -n 5 ) || true
echo

# ----------------------------
# Download parts (resumable)
# ----------------------------
for p in "${PARTS[@]}"; do
  out="${BUILD}.tar.gz.${p}"
  url="${URL_BASE}/${BUILD}.tar.gz.${p}/download"

  echo "[$(date)] Downloading $out"
  # -c resume, --tries/--waitretry robust, -O enforce expected filename
  wget -c --tries=20 --waitretry=10 --timeout=60 --continue -O "$out" "$url"

  # Verify the downloaded part is not an HTML redirect/error page
  # (SourceForge sometimes returns HTML if blocked)
  if file "$out" | grep -qi "HTML"; then
    echo "ERROR: $out looks like HTML (download blocked or redirected improperly)."
    echo "First 5 lines:"
    head -n 5 "$out" || true
    exit 2
  fi

  echo "[$(date)] Done: $out  ($(du -h "$out" | awk '{print $1}'))"
  echo
done

# ----------------------------
# Extract
# ----------------------------
echo "[$(date)] Extracting to ${DBROOT}/${BUILD}/ ..."
# If folder exists and is non-empty, avoid accidental overwrites
if [[ -d "$BUILD" ]] && [[ "$(ls -A "$BUILD" 2>/dev/null | wc -l)" -gt 0 ]]; then
  echo "ERROR: Directory '$BUILD' already exists and is not empty. Remove/rename it first."
  exit 3
fi

cat "${BUILD}.tar.gz."* | tar xz

# ----------------------------
# Symlink "current"
# ----------------------------
ln -sfn "$BUILD" current

# ----------------------------
# Basic sanity checks
# ----------------------------
echo "[$(date)] Sanity checks:"
echo "Listing top of current/:"
ls -lah current | head -n 50

echo
echo "Disk usage:"
du -sh "$BUILD" || true

echo
echo "Sample subdirs:"
find current -maxdepth 2 -type d | head -n 40

echo
echo "[$(date)] DONE. DB ready at: $DBROOT/current"
