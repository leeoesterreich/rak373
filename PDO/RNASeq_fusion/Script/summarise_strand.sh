#!/bin/bash
# Usage: Script/summarise_strand.sh <folder with <sample>/lib_format_counts.json>
# Writes <folder>/strand_summary.tsv. Call: reverse/forward if fraction > 0.8, else unstranded_or_mixed.
BASE=${BASE:-/ix1/alee/LO_LAB/Personal/Rahul/PDO/RNASeq_fusion}
IN_DIR=${1:-${BASE}/strand_check}
OUT=${IN_DIR}/strand_summary.tsv

python3 - "${IN_DIR}" "${OUT}" <<'PYEOF'
import glob, json, os, sys
in_dir, out = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(os.path.join(in_dir, "*", "lib_format_counts.json")))
if not files:
    sys.exit(f"ERROR: no lib_format_counts.json under {in_dir}")
cols = ["sample","expected_format","forward_ISF","reverse_ISR","unstranded_IU",
        "single_fwd_SF","single_rev_SR","forward_fraction","reverse_fraction",
        "compatible_fragment_ratio","call"]
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for f in files:
        d = json.load(open(f))
        isf, isr, iu = d.get("ISF",0), d.get("ISR",0), d.get("IU",0)
        tot = isf + isr + iu
        ff = isf/tot if tot else 0.0
        rf = isr/tot if tot else 0.0
        call = "reverse" if rf > 0.8 else "forward" if ff > 0.8 else "unstranded_or_mixed"
        row = [os.path.basename(os.path.dirname(f)), d.get("expected_format","NA"),
               isf, isr, iu, d.get("SF",0), d.get("SR",0),
               f"{ff:.3f}", f"{rf:.3f}",
               f"{d.get('compatible_fragment_ratio',float('nan')):.3f}", call]
        fh.write("\t".join(map(str,row)) + "\n")
PYEOF

echo "Written: ${OUT}"
column -t -s $'\t' "${OUT}"
