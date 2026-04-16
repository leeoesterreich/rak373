from pysam import VariantFile
import os
import argparse

def reformat_vcf(vcf_file, out, reference, tumor_sample=None):
    os.system(f'bcftools norm -f {reference} -m -both {vcf_file} -o tmp_.vcf')

    with VariantFile('tmp_.vcf') as fr:
        header = fr.header
        header.info.add('TDP', number=1, type='Integer', description='Tumor sample depth')
        header.info.add('NDP', number=1, type='Integer', description='Normal sample depth')
        header.info.add('TAF', number=1, type='Float', description='Tumor sample AF')
        header.info.add('NAF', number=1, type='Float', description='Normal sample AF')
        samples = list(header.samples)

        if tumor_sample is not None:
            if tumor_sample not in [0, 1, '0', '1']:
                if tumor_sample in samples:
                    tumor_idx = samples.index(tumor_sample)
                    normal_idx = 1 - tumor_idx
                else:
                    raise Exception(f'{tumor_sample} is not in samples {samples} recorded in vcf')
            else:
                tumor_idx = int(tumor_sample)
                normal_idx = 1 - tumor_idx
        else:
            tumor_idx = guess_tumor_idx(vcf_file)
            normal_idx = 1 - tumor_idx

        with VariantFile(out, 'w', header=header) as fw:
            for record in fr:
                tdp = record.samples[tumor_idx].get('DP', 0)
                ndp = record.samples[normal_idx].get('DP', 0)

                record.info['TDP'] = tdp
                record.info['NDP'] = ndp

                try:
                    taf = round(record.samples[tumor_idx]['AD'][1] / tdp, 3) if tdp else 0.0
                except (IndexError, TypeError):
                    taf = 0.0

                try:
                    naf = round(record.samples[normal_idx]['AD'][1] / ndp, 3) if ndp else 0.0
                except (IndexError, TypeError):
                    naf = 0.0

                record.info['TAF'] = taf
                record.info['NAF'] = naf

                fw.write(record)

    os.remove('tmp_.vcf')
    os.system(f'bgzip {out}')
    os.system(f'tabix {out}.gz')

def guess_tumor_idx(vcf_file):
    tumor_is_first = 0
    tumor_is_second = 0

    with VariantFile(vcf_file) as fr:
        samples = list(fr.header.samples)
        formats = list(fr.header.formats)
        if 'AF' not in formats:
            raise Exception('No AF in format info to detect tumor sample')

        for record in fr:
            if record.samples[0]['AF'][0] > record.samples[1]['AF'][0]:
                tumor_is_first += 1
            else:
                tumor_is_second += 1

    tumor_idx = tumor_is_second >= tumor_is_first
    print(f'Inferred tumor sample is {samples[tumor_idx]}')
    return tumor_idx

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Reformat a VCF for PCGR compatibility.")
    parser.add_argument("-vcf_file", required=True, help="Input VCF file (decompressed)")
    parser.add_argument("-out", required=True, help="Output VCF file (uncompressed, will be bgzipped)")
    parser.add_argument("-reference", required=True, help="Reference FASTA file")
    parser.add_argument("-tumor_sample", required=False, help="Tumor sample name or index (0/1)")

    args = parser.parse_args()

    reformat_vcf(
        vcf_file=args.vcf_file,
        out=args.out,
        reference=args.reference,
        tumor_sample=args.tumor_sample
    )
