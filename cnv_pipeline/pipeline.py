import os
import shlex
import subprocess
import argparse

import configparser

from cnv_pipeline.baf_from_vcf import baf_from_vcf
from cnv_pipeline.build_coverage_files import build_genome_file, build_coverage_files

from cnv_pipeline.get_loh_intervals_adtex import finalize_loh


this_dir = os.path.dirname(os.path.realpath(__file__))
config = configparser.ConfigParser()


def run_cnv(vcf_path, sample_dir=None, adtex_dir=None, tumor_bam=None, normal_bam=None,
            baf_path=None, feather_path=None, tumor_cov_path=None, normal_cov_path=None,
            col_normal=10, col_tumor=11,
            format_dp_index=5, mq_cutoff=30, chroms=None,
            target_path=None, ploidy=2, min_read_depth=10):
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
    if adtex_dir is None:
        adtex_dir = os.path.join(sample_dir, 'adtex_output')
    if chroms is None:
        chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'
    if target_path is None:
        config.read(os.path.join(this_dir, 'config.ini'))
        target_path = config.get('paths', 'CODING_REGIONS')
    genome_path = os.path.join(sample_dir, "genome.txt")
    if not os.path.exists(sample_dir):
        os.mkdir(sample_dir)

    baf_from_vcf(vcf_path, baf_path, feather_path=feather_path,
                 col_tumor=col_tumor, col_normal=col_normal,
                 format_dp_index=format_dp_index, mq_cutoff=mq_cutoff,
                 chroms=chroms)

    build_genome_file(sample_bam=normal_bam, genome_path=genome_path)

    build_coverage_files(tumor_bam=tumor_bam, normal_bam=normal_bam, genome_path=genome_path,
                         tumor_cov_path=tumor_cov_path, normal_cov_path=normal_cov_path,
                         target_bed_path=target_path)

    run_adtex(normal_cov_path=normal_cov_path,
              tumor_cov_path=tumor_cov_path,
              adtex_dir=adtex_dir,
              baf_path=baf_path,
              target_path=target_path,
              ploidy=ploidy, min_read_depth=min_read_depth)

    finalize_loh(adtex_dir)


def run_adtex(normal_cov_path=None, tumor_cov_path=None, adtex_dir=None, baf_path=None, target_path=None,
              stdout_path=None, ploidy=2, min_read_depth=10):
    """ Example call from bash:
        python2 ADTEx_sgg.py --DOC \
        -n ${pdir}/cov_normal.bed \
        -t ${pdir}/cov_tumor.bed \
        -o ${pdir}/adtex_try1 --baf ${pdir}/baf.txt \
        --bed $coding_bed --estimatePloidy --plot \
        > ${pdir}/run_info.txt 2>&1
    """
    if stdout_path is None:
        stdout_path = os.path.join(adtex_dir, 'run_info.txt')

    config.read(os.path.join(this_dir, 'config.ini'))

    python2_path = config.get('paths', 'PYTHON2')
    adtex_script = config.get('paths', 'ADTEX')

    cmd = ("{python2} {adtex_script} --DOC -n {normal_cov_path} -t {tumor_cov_path} "
           "-o {adtex_dir} --baf {baf_path} --bed {target_path} --estimatePloidy --plot "
           "--ploidy {ploidy} --min_read_depth {mrd}")
    cmd = cmd.format(python2=python2_path,
                     adtex_script=adtex_script,
                     normal_cov_path=normal_cov_path,
                     tumor_cov_path=tumor_cov_path,
                     adtex_dir=adtex_dir,
                     baf_path=baf_path,
                     target_path=target_path,
                     ploidy=ploidy, mrd=min_read_depth)
    print("Running ADTEx with command:\n{cmd}")
    args = shlex.split(cmd)
    with open(stdout_path, 'w') as outfile:
        proc = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=outfile, stderr=outfile)
        proc.communicate()
    print("ADTEx run complete")


def main():
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
