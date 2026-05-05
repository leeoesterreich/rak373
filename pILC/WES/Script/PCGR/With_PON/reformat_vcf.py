from pysam import VariantFile
import os
import argparse
import subprocess


def run_cmd(cmd):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def add_info_if_missing(header, key, number, vtype, description):
    if key not in header.info:
        header.info.add(key, number=number, type=vtype, description=description)


def get_depth(sample_data):
    dp = sample_data.get("DP", 0)
    if dp is None:
        return 0
    try:
        return int(dp)
    except Exception:
        return 0


def get_alt_af(sample_data, dp):
    """
    Prefer FORMAT/AF if present.
    If AF is unavailable, calculate AF from AD[1] / DP.
    """
    af = sample_data.get("AF", None)

    if af is not None:
        try:
            if isinstance(af, (list, tuple)):
                return round(float(af[0]), 3)
            return round(float(af), 3)
        except Exception:
            pass

    ad = sample_data.get("AD", None)

    try:
        if ad is not None and len(ad) > 1 and dp > 0:
            return round(float(ad[1]) / float(dp), 3)
    except Exception:
        pass

    return 0.0


def reformat_vcf(vcf_file, out, reference, tumor_sample=None):
    tmp_vcf = "tmp_.vcf"

    run_cmd([
        "bcftools", "norm",
        "-f", reference,
        "-m", "-both",
        vcf_file,
        "-o", tmp_vcf
    ])

    with VariantFile(tmp_vcf) as fr:
        header = fr.header

        add_info_if_missing(header, "TDP", 1, "Integer", "Tumor sample depth")
        add_info_if_missing(header, "NDP", 1, "Integer", "Normal sample depth")
        add_info_if_missing(header, "TAF", 1, "Float", "Tumor sample AF")
        add_info_if_missing(header, "NAF", 1, "Float", "Normal sample AF")

        samples = list(header.samples)
        print("Samples in VCF:", samples)

        if len(samples) == 0:
            raise Exception("No sample found in VCF header.")

        if tumor_sample is not None:
            if tumor_sample in samples:
                tumor_idx = samples.index(tumor_sample)
            elif tumor_sample in ["0", "1"]:
                tumor_idx = int(tumor_sample)
            else:
                raise Exception(f"{tumor_sample} is not in samples {samples} recorded in VCF")
        else:
            tumor_idx = 0

        if tumor_idx >= len(samples):
            raise Exception(f"Tumor sample index {tumor_idx} is invalid for samples {samples}")

        # Tumor-only VCF: no normal sample exists
        if len(samples) == 1:
            normal_idx = None
            print("Detected tumor-only VCF. Setting NDP=0 and NAF=0.0")
        else:
            normal_idx = 1 - tumor_idx
            print(f"Detected paired VCF. Tumor={samples[tumor_idx]}, Normal={samples[normal_idx]}")

        with VariantFile(out, "w", header=header) as fw:
            for record in fr:
                tumor_data = record.samples[tumor_idx]

                tdp = get_depth(tumor_data)
                taf = get_alt_af(tumor_data, tdp)

                record.info["TDP"] = tdp
                record.info["TAF"] = taf

                if normal_idx is None:
                    record.info["NDP"] = 0
                    record.info["NAF"] = 0.0
                else:
                    normal_data = record.samples[normal_idx]
                    ndp = get_depth(normal_data)
                    naf = get_alt_af(normal_data, ndp)

                    record.info["NDP"] = ndp
                    record.info["NAF"] = naf

                fw.write(record)

    os.remove(tmp_vcf)

    run_cmd(["bgzip", "-f", out])
    run_cmd(["tabix", "-f", "-p", "vcf", out + ".gz"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat a VCF for PCGR compatibility.")
    parser.add_argument("-vcf_file", required=True, help="Input VCF file")
    parser.add_argument("-out", required=True, help="Output VCF file without .gz")
    parser.add_argument("-reference", required=True, help="Reference FASTA file")
    parser.add_argument("-tumor_sample", required=False, help="Tumor sample name or index")

    args = parser.parse_args()

    reformat_vcf(
        vcf_file=args.vcf_file,
        out=args.out,
        reference=args.reference,
        tumor_sample=args.tumor_sample
    )