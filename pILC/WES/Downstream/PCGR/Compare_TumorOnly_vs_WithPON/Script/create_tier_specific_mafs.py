#!/usr/bin/env python3

import os
import csv
import glob
import subprocess
from pathlib import Path


# ============================================================
# Configuration
# ============================================================

# Change this if your comparison folder is somewhere else.
COMPARE_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Compare_TumorOnly_vs_WithPON"
)

# PCGR output folders
TUMOR_ONLY_PCGR_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/Tumor_only/pcgr_outputs"
)

WITH_PON_PCGR_DIR = Path(
    "/ix1/alee/LO_LAB/Personal/Rahul/pILC/WES/Downstream/PCGR/With_PON/pcgr_outputs"
)

# Tier VCF folders
TIER1_VCF_DIR = COMPARE_DIR / "Tier_1_common"
TIER2_VCF_DIR = COMPARE_DIR / "Tier_2_With_PON_only"
TIER3_VCF_DIR = COMPARE_DIR / "Tier_3_Tumor_only_only"

# Output folders
OUT_DIR = COMPARE_DIR / "Tier_MAF"
TIER1_MAF_DIR = OUT_DIR / "Tier_1_common_maf"
TIER2_MAF_DIR = OUT_DIR / "Tier_2_With_PON_only_maf"
TIER3_MAF_DIR = OUT_DIR / "Tier_3_Tumor_only_only_maf"
COMBINED_DIR = OUT_DIR / "Combined_MAF"

for d in [TIER1_MAF_DIR, TIER2_MAF_DIR, TIER3_MAF_DIR, COMBINED_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# ============================================================
# Helper functions
# ============================================================

def normalize_chr(chrom):
    """
    Normalize chromosome names so chr1 and 1 can match.
    """
    chrom = str(chrom).strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return chrom


def run_cmd(cmd):
    result = subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return result.stdout


def extract_vcf_keys(vcf_file):
    """
    Extract variant keys from VCF:
    CHROM, POS, REF, ALT

    Returns a set of tuples:
    (chrom, pos, ref, alt)
    """
    cmd = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\n",
        str(vcf_file)
    ]

    output = run_cmd(cmd)

    keys = set()
    for line in output.splitlines():
        if not line.strip():
            continue

        chrom, pos, ref, alt = line.rstrip("\n").split("\t")

        # If ALT has comma-separated alleles, split them.
        # In your workflow these should usually already be split,
        # but this makes the script more robust.
        for alt_allele in alt.split(","):
            key = (
                normalize_chr(chrom),
                str(pos),
                ref,
                alt_allele
            )
            keys.add(key)

    return keys


def find_pcgr_maf(pcgr_base_dir, sample):
    """
    Find PCGR MAF file for one sample.
    Expected common pattern:
    <sample>.pcgr.grch38.maf
    """
    sample_dir = pcgr_base_dir / sample

    candidates = [
        sample_dir / f"{sample}.pcgr.grch38.maf",
        sample_dir / f"{sample}.pcgr.grch38.pass.maf",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    maf_files = sorted(sample_dir.glob("*.maf"))
    if maf_files:
        return maf_files[0]

    raise FileNotFoundError(f"No MAF file found for sample {sample} in {sample_dir}")


def maf_row_key(row):
    """
    Build variant key from a MAF row.

    Standard MAF columns:
    Chromosome
    Start_Position
    Reference_Allele
    Tumor_Seq_Allele2
    """
    chrom = normalize_chr(row.get("Chromosome", ""))
    pos = str(row.get("Start_Position", ""))
    ref = row.get("Reference_Allele", "")

    # Most standard MAF files use Tumor_Seq_Allele2 for ALT.
    alt = row.get("Tumor_Seq_Allele2", "")

    # Fallbacks in case PCGR/vcf2maf uses slightly different naming.
    if alt == "" and "Tumor_Seq_Allele1" in row:
        alt = row.get("Tumor_Seq_Allele1", "")

    return chrom, pos, ref, alt


def subset_maf_by_vcf_keys(maf_file, vcf_keys, out_maf):
    """
    Subset a PCGR-generated MAF using variant keys extracted from tier VCF.
    Preserves comment lines and header.
    """
    comment_lines = []
    header = None
    matched_rows = []

    with open(maf_file, "r", newline="") as f:
        for line in f:
            if line.startswith("#"):
                comment_lines.append(line)
                continue

            header = line.rstrip("\n").split("\t")
            break

        if header is None:
            raise ValueError(f"No MAF header found in {maf_file}")

        reader = csv.DictReader(f, fieldnames=header, delimiter="\t")

        for row in reader:
            key = maf_row_key(row)
            if key in vcf_keys:
                matched_rows.append(row)

    with open(out_maf, "w", newline="") as out:
        for line in comment_lines:
            out.write(line)

        writer = csv.DictWriter(
            out,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(matched_rows)

    return len(matched_rows), len(vcf_keys)


def get_sample_from_vcf(vcf_file, tier_suffix):
    name = vcf_file.name
    return name.replace(tier_suffix, "")


def process_tier(tier_name, vcf_dir, tier_suffix, maf_source_dir, out_dir):
    """
    Process one tier.
    """
    print(f"\n==============================")
    print(f"Processing {tier_name}")
    print(f"VCF directory: {vcf_dir}")
    print(f"MAF source:    {maf_source_dir}")
    print(f"Output dir:    {out_dir}")
    print(f"==============================")

    summary_rows = []

    vcf_files = sorted(vcf_dir.glob(f"*{tier_suffix}"))

    if not vcf_files:
        print(f"WARNING: No VCF files found in {vcf_dir}")
        return summary_rows

    for vcf_file in vcf_files:
        sample = get_sample_from_vcf(vcf_file, tier_suffix)

        print(f"\nSample: {sample}")
        print(f"Tier VCF: {vcf_file}")

        try:
            maf_file = find_pcgr_maf(maf_source_dir, sample)
            print(f"Source MAF: {maf_file}")

            vcf_keys = extract_vcf_keys(vcf_file)

            out_maf = out_dir / f"{sample}.{tier_name}.maf"

            matched_count, vcf_key_count = subset_maf_by_vcf_keys(
                maf_file=maf_file,
                vcf_keys=vcf_keys,
                out_maf=out_maf
            )

            print(f"VCF variants: {vcf_key_count}")
            print(f"MAF matched rows: {matched_count}")
            print(f"Output MAF: {out_maf}")

            summary_rows.append({
                "Sample": sample,
                "Tier": tier_name,
                "VCF": str(vcf_file),
                "Source_MAF": str(maf_file),
                "Output_MAF": str(out_maf),
                "VCF_variant_count": vcf_key_count,
                "MAF_matched_count": matched_count,
                "Status": "OK"
            })

        except Exception as e:
            print(f"ERROR processing {sample}: {e}")
            summary_rows.append({
                "Sample": sample,
                "Tier": tier_name,
                "VCF": str(vcf_file),
                "Source_MAF": "NA",
                "Output_MAF": "NA",
                "VCF_variant_count": "NA",
                "MAF_matched_count": "NA",
                "Status": f"ERROR: {e}"
            })

    return summary_rows


def combine_mafs(maf_dir, combined_maf):
    """
    Combine per-sample MAFs into one MAF.
    Keeps header from the first file only.
    Skips comment lines from later files.
    """
    maf_files = sorted(maf_dir.glob("*.maf"))

    if not maf_files:
        print(f"WARNING: No MAF files found to combine in {maf_dir}")
        return

    header_written = False

    with open(combined_maf, "w") as out:
        for maf_file in maf_files:
            with open(maf_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        if not header_written:
                            out.write(line)
                        continue

                    # Header line starts with Hugo_Symbol in standard MAF
                    if line.startswith("Hugo_Symbol"):
                        if not header_written:
                            out.write(line)
                            header_written = True
                        continue

                    out.write(line)

    print(f"Combined MAF written: {combined_maf}")


def write_summary(summary_rows, out_csv):
    fieldnames = [
        "Sample",
        "Tier",
        "VCF",
        "Source_MAF",
        "Output_MAF",
        "VCF_variant_count",
        "MAF_matched_count",
        "Status"
    ]

    with open(out_csv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\nSummary written: {out_csv}")


# ============================================================
# Main
# ============================================================

all_summary_rows = []

# Tier 1: common variants
# Use With_PON MAF as source
all_summary_rows.extend(
    process_tier(
        tier_name="tier1.common",
        vcf_dir=TIER1_VCF_DIR,
        tier_suffix=".tier1.common.vcf.gz",
        maf_source_dir=WITH_PON_PCGR_DIR,
        out_dir=TIER1_MAF_DIR
    )
)

# Tier 2: With_PON-only variants
# Use With_PON MAF as source
all_summary_rows.extend(
    process_tier(
        tier_name="tier2.with_pon_only",
        vcf_dir=TIER2_VCF_DIR,
        tier_suffix=".tier2.with_pon_only.vcf.gz",
        maf_source_dir=WITH_PON_PCGR_DIR,
        out_dir=TIER2_MAF_DIR
    )
)

# Tier 3: Tumor-only-only variants
# Use Tumor-only MAF as source
all_summary_rows.extend(
    process_tier(
        tier_name="tier3.tumor_only_only",
        vcf_dir=TIER3_VCF_DIR,
        tier_suffix=".tier3.tumor_only_only.vcf.gz",
        maf_source_dir=TUMOR_ONLY_PCGR_DIR,
        out_dir=TIER3_MAF_DIR
    )
)

# Write summary
summary_csv = OUT_DIR / "tier_maf_creation_summary.csv"
write_summary(all_summary_rows, summary_csv)

# Combine MAFs across samples
combine_mafs(
    TIER1_MAF_DIR,
    COMBINED_DIR / "Tier_1_common_all_samples.maf"
)

combine_mafs(
    TIER2_MAF_DIR,
    COMBINED_DIR / "Tier_2_With_PON_only_all_samples.maf"
)

combine_mafs(
    TIER3_MAF_DIR,
    COMBINED_DIR / "Tier_3_Tumor_only_only_all_samples.maf"
)

print("\nAll done.")