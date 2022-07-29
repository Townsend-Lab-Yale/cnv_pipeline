import os
import sys
import shlex
import subprocess
import argparse
import contextlib

from .baf_from_vcf import baf_from_vcf
from .build_coverage_files import build_genome_file, build_coverage_files
from .config import load_config, this_dir
from .get_loh_intervals_adtex import finalize_loh
from .trim_vcf import trim_vcf


def run_cnv(vcf_path, sample_dir=None, adtex_dir=None, tumor_bam=None, normal_bam=None,
            baf_path=None, parquet_path=None, tumor_cov_path=None, normal_cov_path=None,
            tumor_id=None, normal_id=None, saas_only=False, adtex_stdout='-',
            mq_cutoff=30, chroms=None, vcf_out=None,
            target_path=None, ploidy=None, min_read_depth=10, sample_id='Tumor',
            ratio_min=0.4, ratio_max=0.6, min_tumor=20, min_normal=10, min_gq=90):
    """Run pipeline.

    At minimum, requires vcf_path and sample_dir (for storing output).

    Note:
        sample_dir (str): Must be sample specific to prevent file overwrite.
    """

    if baf_path is None:
        baf_path = os.path.join(sample_dir, "baf.txt")
    if parquet_path is None:
        parquet_path = os.path.join(sample_dir, "saas.parquet")
    if tumor_cov_path is None:
        tumor_cov_path = os.path.join(sample_dir, "tumor_cov.bed")
    if normal_cov_path is None:
        normal_cov_path = os.path.join(sample_dir, "normal_cov.bed")
    if adtex_dir is None:
        adtex_dir = os.path.join(sample_dir, 'adtex_output')
    if chroms is None:
        chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'
    if target_path is None:
        config = load_config()
        target_path = config.get('paths', 'CODING_REGIONS')
    if vcf_out is None:
        vcf_out = os.path.join(sample_dir, "snps_trimmed.vcf")
    genome_path = os.path.join(sample_dir, "genome.txt")
    if not os.path.exists(sample_dir):
        os.mkdir(sample_dir)
    if normal_id is None:
        normal_id = tumor_id + 'N'

    trim_vcf(vcf_in=vcf_path, tumor_id=tumor_id, normal_id=normal_id,
             ratio_min=ratio_min, ratio_max=ratio_max, min_depth_n=min_normal,
             min_depth_t=min_tumor, min_gq_n=min_gq, vcf_out=vcf_out)

    baf_from_vcf(vcf_out, baf_path, parquet_path=parquet_path,
                 tumor_id=tumor_id, normal_id=normal_id,
                 mq_cutoff=mq_cutoff, chroms=chroms)

    run_saasCNV(sample_id=sample_id, sample_dir=sample_dir,
                baf_path=parquet_path, stdout_path='-')

    if not saas_only:
        build_genome_file(sample_bam=normal_bam, genome_path=genome_path)

        build_coverage_files(tumor_bam=tumor_bam, normal_bam=normal_bam, genome_path=genome_path,
                             tumor_cov_path=tumor_cov_path, normal_cov_path=normal_cov_path,
                             target_bed_path=target_path)

        run_adtex(normal_cov_path=normal_cov_path,
                  tumor_cov_path=tumor_cov_path,
                  adtex_dir=adtex_dir,
                  baf_path=baf_path,
                  target_path=target_path,
                  ploidy=ploidy, min_read_depth=min_read_depth,
                  stdout_path=adtex_stdout)

        finalize_loh(adtex_dir)


def run_saasCNV(sample_id=None, sample_dir=None, baf_path=None, stdout_path='-'):
    """Example call from bash:
    Rscript run_saas.R {s_id} {sample_dir} {baf_path} 50 30 FALSE 0.05 0.05
    """
    script_path = os.path.join(this_dir, 'run_saas.R')
    config = load_config()
    rscript_path = config.get('paths', 'RSCRIPT', fallback='Rscript')
    sample_path = os.path.realpath(sample_dir)
    baf_path = os.path.realpath(baf_path)
    cmd = (f"{rscript_path} {script_path} {sample_id} {sample_path} {baf_path}"
           " 50 30 FALSE 0.05 0.05")
    print("Running saasCNV with command:\n  {}".format(cmd))
    args = shlex.split(cmd)
    with smart_open(stdout_path) as outfile:
        proc = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=outfile, stderr=outfile)
        proc.communicate()
    print("saasCNV run complete")


def run_adtex(normal_cov_path=None, tumor_cov_path=None, adtex_dir=None, baf_path=None, target_path=None,
              stdout_path=None, ploidy=None, min_read_depth=10):
    """ Example call from bash:
        python2 ADTEx_sgg.py --DOC \
        -n ${pdir}/cov_normal.bed \
        -t ${pdir}/cov_tumor.bed \
        -o ${pdir}/adtex_try1 --baf ${pdir}/baf.txt \
        --bed $coding_bed --estimatePloidy --plot \
        > ${pdir}/run_info.txt 2>&1
    """
    if stdout_path is None:
        stdout_path = '-'  # will write to stdout
    ploidy_str = '--ploidy {}'.format(ploidy) if ploidy is not None else ''

    config = load_config()
    python_path = config.get('paths', 'PYTHON', fallback=sys.executable)
    adtex_script = config.get('paths', 'ADTEX', fallback='adtex.py')
    cmd = ("{python_path} {adtex_script} --DOC -n {normal_cov_path} -t {tumor_cov_path} "
           "-o {adtex_dir} --baf {baf_path} --bed {target_path} --estimatePloidy --plot "
           "{ploidy_str} --minReadDepth {mrd}")
    cmd = cmd.format(python_path=python_path,
                     adtex_script=adtex_script,
                     normal_cov_path=normal_cov_path,
                     tumor_cov_path=tumor_cov_path,
                     adtex_dir=adtex_dir,
                     baf_path=baf_path,
                     target_path=target_path,
                     ploidy_str=ploidy_str, mrd=min_read_depth)
    print("Running ADTEx with command:\n  {}".format(cmd))
    args = shlex.split(cmd)
    with smart_open(stdout_path) as outfile:
        proc = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=outfile, stderr=outfile)
        proc.communicate()
    print("ADTEx run complete")


@contextlib.contextmanager
def smart_open(filename=None):
    """For opening file -OR- writing to stdout. via https://stackoverflow.com/a/17603000"""
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def main():
    if __name__ == "__main__" and __package__ is None:
        __package__ = "cnv_pipeline"

    parser = argparse.ArgumentParser("CNV PIPELINE")
    parser.add_argument('-v', '--vcf', help='VCF file for sample pair', required=True)
    parser.add_argument('-s', '--sample_dir', help='Sample-specific output dir', required=True)
    parser.add_argument('-t', '--tumor', help='Tumor BAM', required=False)
    parser.add_argument('-n', '--normal', help='Normal BAM', required=False)
    parser.add_argument('-tid', '--tumor_id', help='Tumor name, for vcf extraction', required=True)
    parser.add_argument('-nid', '--normal_id', help='Normal name, for vcf extraction', required=True)

    parser.add_argument('--saas_only', help='Only run saasCNV, not ADTEx.', action='store_true', default=False)

    parser.add_argument('-rmin', '--ratio_min', help='Min alt ratio in normal sample [0.4]', type=float, default=0.4)
    parser.add_argument('-rmax', '--ratio_max', help='Max alt ratio in normal sample [0.6]', type=float, default=0.6)
    parser.add_argument('-tmin', '--min_tumor', help='Min reads from tumor sample [20]', type=int, default=20)
    parser.add_argument('-nmin', '--min_normal', help='Min reads from normal sample [10]', type=int, default=10)
    parser.add_argument('-gq', '--min_gq', help='Genotype quality cutoff for normal sample [90]', type=int, default=90)
    parser.add_argument('-a', '--adtex_dir', help='ADTEx output dir', default=None)
    parser.add_argument("--ploidy", help="Most common ploidy in the tumour sample", type=int, default=None)
    parser.add_argument("--minReadDepth", help="The ADTEx threshold for minimum read depth for each exon [10]",
                        type=int, default=10)
    parser.add_argument('-ao', '--adtex_stdout', help='ADTEx stdout path if overriding STDOUT', default='-')

    args = parser.parse_args()
    vcf_dict = dict(ratio_min=args.ratio_min, ratio_max=args.ratio_max,
                    min_tumor=args.min_tumor, min_normal=args.min_normal,
                    min_gq=args.min_gq)
    run_cnv(vcf_path=args.vcf, sample_dir=args.sample_dir, adtex_dir=args.adtex_dir,
            tumor_bam=args.tumor, normal_bam=args.normal,
            tumor_id=args.tumor_id, normal_id=args.normal_id,
            saas_only=args.saas_only,
            ploidy=args.ploidy, min_read_depth=args.minReadDepth,
            adtex_stdout=args.adtex_stdout, sample_id=args.tumor_id, **vcf_dict)
