#!/usr/bin/env python3
"""
Merge Arriba, STAR-Fusion and FusionCatcher calls (nf-core/rnafusion 4.1.3 output)
into per-sample consensus fusions, attach FusionInspector read support, flag likely
artifacts and count recurrence across samples.

Usage (one results folder per batch; samples are found from results/arriba):
  python build_consensus_v4.py --results <results_batch1> <results_batch2> <results_batch3> \
         --outdir <dir>
Optional: --samplesheet <sheet1.csv> <sheet2.csv> ... (same order) to take the sample list
from samplesheets instead.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

# ----------------------------------------------------------------------------
# Settings you may want to change
# ----------------------------------------------------------------------------
TOL = 10                    # breakpoints within +/- TOL bp at BOTH ends = same event
READTHROUGH_MAX_DIST = 100_000   # same chromosome and closer than this = readthrough/adjacent
# Gene families / duplicated regions that often give mapping artifacts
PARALOG_REGEX = re.compile(
    r"^(NBPF|NOTCH2NL|SEC22B|GOLGA6|GOLGA8|NPIP|TBC1D3|ARL17|KANSL1|BOLA2|"
    r"SMG1|ANKRD20A|FAM90A|USP17L|TP53TG3|CT45A|SPATA31)")
# FusionCatcher description tags that point to normal/benign or blacklisted events
FC_NORMAL_TAGS = {"banned", "healthy", "non_tumor_cells", "1000genomes", "gtex"}
# FusionCatcher description tags meaning "listed in a cancer database" (label only, not a filter)
FC_KNOWN_TAGS = {"known", "tcga", "18cancers", "mitelman", "ccle3", "tcga3"}
TOOL_PRIORITY = ["starfusion", "arriba", "fusioncatcher"]   # which tool names the event

COMMON_COLS = ["sample", "tool", "gene5", "gene3", "chr5", "pos5", "chr3", "pos3",
               "support", "tags"]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def parse_bp(text):
    """'chr19:52584197:+' or '19:52584197' -> ('19', 52584197)."""
    parts = str(text).split(":")
    return parts[0].replace("chr", ""), int(parts[1])


def to_num(series):
    return pd.to_numeric(series, errors="coerce").fillna(0)


def read_table(path):
    """Read a tab-separated file; return empty frame if file missing or empty."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t", dtype=str)
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def build_frame(rows):
    df = pd.DataFrame(rows, columns=COMMON_COLS)
    for c in ("gene5", "gene3", "tags"):          # empty cells arrive as NaN; keep them as text
        df[c] = df[c].fillna("").astype(str)
    return df


# ----------------------------------------------------------------------------
# One reader per tool: all return rows in the same common layout
# ----------------------------------------------------------------------------
def read_arriba(path, sample):
    df, rows = read_table(path), []
    for _, r in df.iterrows():
        c5, p5 = parse_bp(r["breakpoint1"])
        c3, p3 = parse_bp(r["breakpoint2"])
        support = sum(to_num(pd.Series([r.get(k)])).iloc[0]
                      for k in ("split_reads1", "split_reads2", "discordant_mates"))
        rows.append([sample, "arriba", r["gene1"], r["gene2"], c5, p5, c3, p3,
                     support, "confidence=" + str(r.get("confidence"))])
    return build_frame(rows)


def read_starfusion(path, sample):
    df, rows = read_table(path), []
    for _, r in df.iterrows():
        c5, p5 = parse_bp(r["LeftBreakpoint"])
        c3, p3 = parse_bp(r["RightBreakpoint"])
        support = float(r["JunctionReadCount"]) + float(r["SpanningFragCount"])
        rows.append([sample, "starfusion", r["LeftGene"].split("^")[0],
                     r["RightGene"].split("^")[0], c5, p5, c3, p3, support, ""])
    return build_frame(rows)


def read_fusioncatcher(path, sample):
    df, rows = read_table(path), []
    for _, r in df.iterrows():
        c5, p5 = parse_bp(r["Fusion_point_for_gene_1(5end_fusion_partner)"])
        c3, p3 = parse_bp(r["Fusion_point_for_gene_2(3end_fusion_partner)"])
        support = (float(r["Spanning_pairs"]) + float(r["Spanning_unique_reads"]))
        rows.append([sample, "fusioncatcher", r["Gene_1_symbol(5end_fusion_partner)"],
                     r["Gene_2_symbol(3end_fusion_partner)"], c5, p5, c3, p3, support,
                     r["Fusion_description"]])
    return build_frame(rows)


def read_fusioninspector(path, sample):
    """FusionInspector read support: used only as extra evidence."""
    df, rows = read_table(path), []
    for _, r in df.iterrows():
        c5, p5 = parse_bp(r["LeftBreakpoint"])
        c3, p3 = parse_bp(r["RightBreakpoint"])
        rows.append([sample, c5, p5, c3, p3, float(r["JunctionReadCount"]),
                     float(r["SpanningFragCount"])])
    return pd.DataFrame(rows, columns=["sample", "chr5", "pos5", "chr3", "pos3",
                                       "fi_junction", "fi_spanning"])


# ----------------------------------------------------------------------------
# Greedy clustering by breakpoint (works on a list of (chr5,pos5,chr3,pos3))
# ----------------------------------------------------------------------------
def assign_clusters(coords):
    reps, ids, n_clusters = {}, [], 0    # reps: (chr5, chr3) -> [(pos5, pos3, id)]
    for c5, p5, c3, p3 in coords:
        found = None
        for q5, q3, cid in reps.get((c5, c3), []):
            if abs(p5 - q5) <= TOL and abs(p3 - q3) <= TOL:
                found = cid
                break
        if found is None:               # new cluster
            found = n_clusters
            n_clusters += 1
            reps.setdefault((c5, c3), []).append((p5, p3, found))
        ids.append(found)
    return ids


# ----------------------------------------------------------------------------
# Per-sample consensus
# ----------------------------------------------------------------------------
def consensus_for_sample(calls, fi):
    calls = calls.reset_index(drop=True)
    calls["event"] = assign_clusters(zip(calls.chr5, calls.pos5, calls.chr3, calls.pos3))
    out = []
    for _, grp in calls.groupby("event"):
        grp = grp.copy()
        grp["prio"] = grp["tool"].map(TOOL_PRIORITY.index)
        best = grp.sort_values("prio").iloc[0]            # call that names the event
        tools = sorted(grp["tool"].unique())
        rec = {"sample": best["sample"], "gene5": best["gene5"], "gene3": best["gene3"],
               "chr5": best["chr5"], "pos5": best["pos5"],
               "chr3": best["chr3"], "pos3": best["pos3"],
               "tools": ",".join(tools), "n_tools": len(tools),
               "all_gene_names": "|".join(sorted({f"{a}--{b}" for a, b in
                                                  zip(grp.gene5, grp.gene3)}))}
        for t in TOOL_PRIORITY:                           # support per tool
            sub = grp[grp.tool == t]
            rec[f"support_{t}"] = sub["support"].max() if len(sub) else 0
        # FusionInspector: best support among rows matching both breakpoints
        m = fi[(fi.chr5 == best["chr5"]) & (fi.chr3 == best["chr3"]) &
               ((fi.pos5 - best["pos5"]).abs() <= TOL) &
               ((fi.pos3 - best["pos3"]).abs() <= TOL)]
        rec["fi_junction"] = m["fi_junction"].max() if len(m) else float("nan")
        rec["fi_spanning"] = m["fi_spanning"].max() if len(m) else float("nan")
        # Artifact flags (information only; nothing is removed)
        rec["flag_readthrough"] = bool(best["chr5"] == best["chr3"] and
                                       abs(best["pos5"] - best["pos3"]) < READTHROUGH_MAX_DIST)
        rec["flag_paralog"] = bool(PARALOG_REGEX.match(str(best["gene5"])) or
                                   PARALOG_REGEX.match(str(best["gene3"])) or
                                   any(PARALOG_REGEX.match(str(g)) for g in
                                       "|".join(grp.gene5).split("|") + "|".join(grp.gene3).split("|")))
        fc_tags = set(",".join(grp[grp.tool == "fusioncatcher"]["tags"]).split(","))
        rec["flag_fc_normal"] = bool(fc_tags & FC_NORMAL_TAGS)
        # Known-fusion label: database tags FusionCatcher reported, e.g. "known;mitelman;tcga".
        # Empty = FusionCatcher gave no database tag (Arriba/STAR-Fusion do not write them)
        rec["known_fusion_db"] = ";".join(sorted({t.strip() for t in fc_tags} & FC_KNOWN_TAGS))
        rec["flag_unannotated"] = bool(str(best["gene5"]).startswith("ENSG") or
                                       str(best["gene3"]).startswith("ENSG"))
        out.append(rec)
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", nargs="+", required=True,
                    help="results folder of each batch")
    ap.add_argument("--samplesheet", nargs="+", default=None,
                    help="optional: samplesheet of each batch, same order as --results")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min_tools", type=int, default=2)
    a = ap.parse_args()

    if a.samplesheet and len(a.results) != len(a.samplesheet):
        sys.exit("ERROR: give one samplesheet per results folder, in the same order")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # (results folder, sample, batch label) for every sample of every batch
    jobs = []
    for i, res in enumerate(a.results):
        R = Path(res)
        label = R.parent.name
        if a.samplesheet:
            samples_b = pd.read_csv(a.samplesheet[i], dtype=str)["sample"].drop_duplicates().tolist()
        else:                            # sample names from the Arriba result files
            suffix = ".arriba.fusions.tsv"
            samples_b = sorted(p.name[:-len(suffix)] for p in (R / "arriba").glob("*" + suffix))
        if not samples_b:
            sys.exit(f"ERROR: no samples found for {R}")
        print(f"{label}: {len(samples_b)} samples")
        for s in samples_b:
            jobs.append((R, s, label))
    names = [j[1] for j in jobs]
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        sys.exit(f"ERROR: sample names used in more than one batch: {dups}")
    samples = names

    all_calls, per_sample = [], []
    for R, s, label in jobs:
        parts = [
            read_arriba(R / "arriba" / f"{s}.arriba.fusions.tsv", s),
            read_starfusion(R / "starfusion" / f"{s}.starfusion.fusion_predictions.tsv", s),
            read_fusioncatcher(R / "fusioncatcher" / f"{s}.fusion-genes.txt", s),
        ]
        parts = [p for p in parts if len(p)]
        if not parts:                    # no caller reported anything for this sample
            print(f"WARNING: no calls for {s}", file=sys.stderr)
            continue
        calls = pd.concat(parts, ignore_index=True)
        calls["batch"] = label
        fi = read_fusioninspector(R / "fusioninspector" / s /
                                  f"{s}.FusionInspector.fusions.tsv", s)
        res_s = consensus_for_sample(calls, fi)
        res_s["batch"] = label
        all_calls.append(calls)
        per_sample.append(res_s)

    all_calls = pd.concat(all_calls, ignore_index=True)
    cons = pd.concat(per_sample, ignore_index=True)

    # Recurrence: same breakpoint (+/- TOL) in different samples = same global event
    cons["event_id"] = assign_clusters(zip(cons.chr5, cons.pos5, cons.chr3, cons.pos3))
    n_samples = cons.groupby("event_id")["sample"].nunique()
    cons["n_samples_with_event"] = cons["event_id"].map(n_samples)
    flags = ["flag_readthrough", "flag_paralog", "flag_fc_normal", "flag_unannotated"]
    cons["any_artifact_flag"] = cons[flags].any(axis=1)

    rec = (cons.groupby("event_id")
           .agg(gene5=("gene5", "first"), gene3=("gene3", "first"),
                chr5=("chr5", "first"), pos5=("pos5", "first"),
                chr3=("chr3", "first"), pos3=("pos3", "first"),
                n_samples=("sample", "nunique"),
                samples=("sample", lambda x: ",".join(sorted(set(x)))),
                batches=("batch", lambda x: ",".join(sorted(set(x)))),
                max_tools=("n_tools", "max"),
                any_artifact_flag=("any_artifact_flag", "max"))
           .sort_values("n_samples", ascending=False).reset_index())

    # Shortlist = called by at least --min_tools callers. Artifact flags, FusionInspector reads,
    # known_fusion_db tags and recurrence are kept as columns for the reader to judge, not as filters
    short = cons[cons.n_tools >= a.min_tools]
    short = short.sort_values(["n_tools", "fi_junction"], ascending=False)

    all_calls.to_csv(out / "fusions_all_calls.tsv", sep="\t", index=False)
    cons.to_csv(out / "fusions_consensus_per_sample.tsv", sep="\t", index=False)
    rec.to_csv(out / "fusions_event_recurrence.tsv", sep="\t", index=False)
    short.to_csv(out / "fusions_shortlist.tsv", sep="\t", index=False)

    print(f"samples read           : {all_calls['sample'].nunique()} of {len(samples)}")
    print(f"raw calls              : {len(all_calls)}")
    print(f"per-sample events      : {len(cons)}")
    print(f"events with >=2 tools  : {(cons.n_tools >= 2).sum()}")
    print(f"events with >=3 tools  : {(cons.n_tools >= 3).sum()}")
    print(f"flagged as artifact    : {cons.any_artifact_flag.sum()}")
    print(f"known_fusion_db labelled: {(cons.known_fusion_db != '').sum()}")
    print(f"shortlist rows         : {len(short)}")


if __name__ == "__main__":
    main()