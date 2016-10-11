import os
import sys

from baf_from_vcf import baf_from_vcf
from build_coverage_files import build_genome_file, build_coverage_files
from get_loh_intervals_adtex import finalize_loh


def run_cnv(vcf_path, sample_dir=None, adtex_dir=None, tumor_bam=None, normal_bam=None,
            baf_path=None, feather_path=None, tumor_cov_path=None, normal_cov_path=None,
            col_normal=10, col_tumor=11,
            format_dp_index=5, mq_cutoff=30, chroms=None,
            target_path=None):
    """Run pipeline.

    At minimum, requires vcf_path and sample_dir (for storing output).

    Args:
        sample_dir (str): Must be sample specific to prevent file overwrite.
    """

    if baf_path is None:
        baf_path = os.path.join(sample_dir, "baf.txt")
    if feather_path is None:
        feather_path = os.path.join(sample_dir, "saas.feather")
    if tumor_cov_path is None:
        tumor_cov_path = os.path.join(sample_dir, "tumor_cov.bed")
    if normal_cov_path is None:
        normal_cov_path = os.path.join(sample_dir, "normal_cov.bed")
    if chroms is None:
        chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'
    if target_path is None:
        os.environ.get('CODING_REGIONS')
    genome_path = os.path.join(sample_dir, "genome.txt")

    baf_from_vcf(vcf_path, baf_path, feather_path=feather_path,
                 col_tumor=col_tumor, col_normal=col_normal,
                 format_dp_index=format_dp_index, mq_cutoff=mq_cutoff,
                 chroms=chroms)

    build_genome_file(sample_bam=normal_bam, genome_path=genome_path)

    build_coverage_files(tumor_bam=tumor_bam, normal_bam=normal_bam, genome_path=genome_path,
                         tumor_cov_path=tumor_cov_path, normal_cov_path=normal_cov_path,
                         target_bed_path=target_path)

    run_adtex()

    finalize_loh(adtex_dir)


def run_adtex():
    pass


if __name__ == '__main__':
    vcf_path = sys.argv[1]
    run_cnv(vcf_path)
