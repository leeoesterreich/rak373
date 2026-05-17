#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import py2bit
import statsmodels.api as sm

# ------------------------------------------------------------
# Calculate GC content for each filtered bin + LOESS correction
# ------------------------------------------------------------

INDIR = "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/Fragmentation_analysis/Regional_Fragment_analysis/02_crpb"

SHORT_IN = os.path.join(
    INDIR,
    "short_rawdata_countreadsperbin.blacklistFiltered.tsv"
)

LONG_IN = os.path.join(
    INDIR,
    "long_rawdata_countreadsperbin.blacklistFiltered.tsv"
)

TWOBIT = "/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Human/Human_ref_hg38_p7/hg38.2bit"

OUTDIR = INDIR


def gc_fraction_from_2bit(tb, chrom, start, end):
    """
    Calculate GC fraction for one genomic interval using py2bit.
    Handles chr/non-chr chromosome naming differences.
    """

    chrom = str(chrom)
    start = int(start)
    end = int(end)

    chroms = tb.chroms()

    chrom_use = chrom

    if chrom_use not in chroms:
        if chrom.startswith("chr"):
            alt_chrom = chrom.replace("chr", "", 1)
        else:
            alt_chrom = "chr" + chrom

        if alt_chrom in chroms:
            chrom_use = alt_chrom
        else:
            return np.nan

    chrom_size = chroms[chrom_use]

    start = max(0, start)
    end = min(end, chrom_size)

    if end <= start:
        return np.nan

    try:
        seq = tb.sequence(chrom_use, start, end)
    except Exception:
        return np.nan

    if seq is None:
        return np.nan

    seq = seq.upper()

    valid = (
        seq.count("A") +
        seq.count("C") +
        seq.count("G") +
        seq.count("T")
    )

    if valid == 0:
        return np.nan

    gc_count = seq.count("G") + seq.count("C")

    return gc_count / valid


def gc_loess_correct(tb, in_tsv, out_prefix, frac=0.75):
    print("\nProcessing:", in_tsv, flush=True)

    df = pd.read_csv(in_tsv, sep="\t")

    df = df.rename(
        columns={
            df.columns[0]: "chrom",
            df.columns[1]: "start",
            df.columns[2]: "end"
        }
    )

    df["start"] = pd.to_numeric(df["start"], errors="raise").astype(int)
    df["end"] = pd.to_numeric(df["end"], errors="raise").astype(int)

    samplenames = list(df.columns[3:])

    print("  bins:", df.shape[0], flush=True)
    print("  samples:", len(samplenames), flush=True)

    # ---- GC content ----
    gc = [
        gc_fraction_from_2bit(tb, c, s, e)
        for c, s, e in zip(df["chrom"], df["start"], df["end"])
    ]

    df["GC"] = gc

    print("  GC NaNs:", int(df["GC"].isna().sum()), "out of", df.shape[0], flush=True)

    df = df.dropna(subset=["GC"]).copy()

    print("  bins after GC filter:", df.shape[0], flush=True)

    # Mean readcount across samples
    df["mean_readcount_over_all_samples"] = df[samplenames].mean(axis=1)

    # Sort deterministically by GC and genomic coordinate
    df = df.sort_values(["GC", "chrom", "start", "end"]).reset_index(drop=True)

    # ---- LOESS per sample ----
    corrected = {}

    for i, sample in enumerate(samplenames, 1):
        y = df[sample].values
        x = df["GC"].values

        loess = sm.nonparametric.lowess(
            endog=y,
            exog=x,
            frac=frac,
            is_sorted=True,
            return_sorted=False
        )

        corrected[sample] = y - loess + np.median(y)

        if i % 25 == 0 or i == len(samplenames):
            print("    LOESS done:", i, "/", len(samplenames), flush=True)

    # ---- Save outputs ----
    df_corr = pd.DataFrame(corrected)

    df_corr.insert(0, "end", df["end"].values)
    df_corr.insert(0, "start", df["start"].values)
    df_corr.insert(0, "chrom", df["chrom"].values)

    # Full table with GC
    out1 = out_prefix + "_readcounts_corr_and_GC.csv"

    df_out = df[[
        "chrom",
        "start",
        "end",
        "GC",
        "mean_readcount_over_all_samples"
    ]].join(df_corr[samplenames])

    df_out.to_csv(out1, index=False)

    # Matrix form: samples x bins
    df_corr["bin"] = (
        df["chrom"].astype(str) + "_" +
        df["start"].astype(str) + "_" +
        df["end"].astype(str)
    )

    mat = df_corr.set_index("bin")[samplenames].transpose()

    out2 = out_prefix + "_corr_readcounts.csv"
    mat.to_csv(out2)

    print("  wrote:", flush=True)
    print("   ", out1, flush=True)
    print("   ", out2, flush=True)

    return mat


def main():
    print("Input short file:", SHORT_IN, flush=True)
    print("Input long file :", LONG_IN, flush=True)
    print("Reference 2bit  :", TWOBIT, flush=True)

    # Quick checks
    df_test = pd.read_csv(SHORT_IN, sep="\t", nrows=5)
    print("\nFirst 5 bins from short file:", flush=True)
    print(df_test.iloc[:, 0:3], flush=True)

    tb = py2bit.open(TWOBIT)

    try:
        chrom_keys = list(tb.chroms().keys())[:20]
        print("\nFirst 20 chroms in 2bit:", chrom_keys, flush=True)
        print("chr1 in 2bit?", "chr1" in tb.chroms(), flush=True)
        print("1 in 2bit?", "1" in tb.chroms(), flush=True)

        short_mat = gc_loess_correct(
            tb,
            SHORT_IN,
            os.path.join(OUTDIR, "short")
        )

        long_mat = gc_loess_correct(
            tb,
            LONG_IN,
            os.path.join(OUTDIR, "long")
        )

        print("\nshort matrix:", short_mat.shape, flush=True)
        print("long matrix :", long_mat.shape, flush=True)

        assert short_mat.index.equals(long_mat.index), "Short and long sample names are not aligned"
        assert short_mat.columns.equals(long_mat.columns), "Short and long genomic bins are not aligned"

        print("\nShort/long matrices aligned", flush=True)

    finally:
        tb.close()


if __name__ == "__main__":
    main()
