#!/usr/bin/env python3
import os
import sys
import deeptools.countReadsPerBin as crpb

chunk_list = sys.argv[1]          # file with BAMs (one per line)
out_tsv    = sys.argv[2]          # output tsv
mode       = sys.argv[3]          # "short" or "long"
threads    = int(sys.argv[4])     # e.g. 10
blacklist  = sys.argv[5]          # bed file (unused here)
binbed     = sys.argv[6]          # fixed 100kb bins BED path

min_len = 90  if mode == "short" else 151
max_len = 150 if mode == "short" else 220

with open(chunk_list) as f:
    bamfiles = [x.strip() for x in f if x.strip()]

# sample names from BAM filenames
samplenames = [os.path.basename(b).replace(".filter.coordsort.bam", "") for b in bamfiles]

cr = crpb.CountReadsPerBin(
    bamfiles,
    bedFile=[binbed],             # pass path, not file handle
    minMappingQuality=30,
    ignoreDuplicates=False,
    numberOfProcessors=threads,
    verbose=False,
    samFlag_exclude=256,
    extendReads=True,
    minFragmentLength=min_len,
    maxFragmentLength=max_len,
    out_file_for_raw_data=out_tsv,
    blackListFileName=None,       # keep all bins; mask later after merge
)

cr.run()

# add header row
tmp = out_tsv + ".tmp"
with open(out_tsv, "r") as fin, open(tmp, "w") as fout:
    fout.write("\t".join(["chrom", "start", "end"] + samplenames) + "\n")
    for line in fin:
        fout.write(line)

os.replace(tmp, out_tsv)
print("Wrote (with header):", out_tsv)
