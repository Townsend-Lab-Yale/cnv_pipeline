import os
import sys

from baf_from_vcf import baf_from_vcf
from get_loh_intervals_adtex import finalize_loh


def run_cnv(vcf_path, adtex_dir=None,
            baf_path=None, feather_path=None, col_normal=10, col_tumor=11,
            format_dp_index=5, mq_cutoff=30, chroms=None):

    if baf_path is None:
        baf_path = "{}_baf.txt".format(os.path.splitext(vcf_path)[0])
    if feather_path is None:
        feather_path = "{}_saas.feather".format(os.path.splitext(vcf_path)[0])
    if chroms is None:
        chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'

    baf_from_vcf(vcf_path, baf_path, feather_path=feather_path, col_tumor=col_tumor, col_normal=col_normal,
                 format_dp_index=format_dp_index, mq_cutoff=mq_cutoff,
                 chroms=chroms)

    build_coverage()

    run_adtex()

    finalize_loh(adtex_dir)


def build_coverage():
    pass


def run_adtex():
    pass






if __name__ == '__main__':
    vcf_path = sys.argv[1]
    run_cnv(vcf_path)
