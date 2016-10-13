#!/usr/bin/env python
import argparse

from cnv_pipeline.pipeline import run_cnv

if __name__ == '__main__':
    parser = argparse.ArgumentParser("CNV PIPELINE")
    parser.add_argument('-v', '--vcf', help='VCF file for sample pair [REQUIRED]')
    parser.add_argument('-s', '--sample_dir', help='Sample-specific intermediate output dir [OPTIONAL]', required=True)
    parser.add_argument('-t', '--tumor', help='Tumor BAM [REQUIRED]', required=True)
    parser.add_argument('-n', '--normal', help='Normal BAM [REQUIRED]', required=True)
    parser.add_argument('-a', '--adtex_dir', help='ADTEx output dir', default=None)
    parser.add_argument('-ct', '--col_tumor', help='VCF column index for tumor, 1-based', type=int, default=10)
    parser.add_argument('-cn', '--col_normal', help='VCF column index for tumor, 1-based', type=int, default=11)
    parser.add_argument("--ploidy", help="Most common ploidy in the tumour sample [2]", type=int, default=2)
    parser.add_argument("--minReadDepth", help="The threshold for minimum read depth for each exon [10]",
                        type=int, default=10)

    args = parser.parse_args()
    run_cnv(vcf_path=args.vcf, sample_dir=args.sample_dir, adtex_dir=args.adtex_dir,
            tumor_bam=args.tumor, normal_bam=args.normal,
            col_tumor=args.col_tumor, col_normal=args.col_normal,
            ploidy=args.ploidy, min_read_depth=args.minReadDepth)
