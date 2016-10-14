import subprocess
import shlex

from cnv_pipeline.config import load_config


def build_genome_file(sample_bam=None, genome_path=None):
    """Build genome file for use with bedtools coverage --sorted mode."""

    config = load_config()
    sam_path = config.get('paths', 'SAMTOOLS', 'samtools')

    """samtools view -H {bam} | grep -P "@SQ\tSN:" | sed 's/@SQ\tSN://' | sed 's/\tLN:/\t/' > {genome_path}"""
    cmd1 = "{sam_path} view -H {bam}".format(sam_path=sam_path, bam=sample_bam)
    cmd2 = 'grep -P "@SQ\tSN:"'
    cmd3 = "sed 's/@SQ\tSN://'"
    cmd4 = "sed 's/\tLN:/\t/'"

    with open(genome_path, 'w') as genome_file:
        p1 = subprocess.Popen(shlex.split(cmd1), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(shlex.split(cmd2), stdout=subprocess.PIPE, stdin=p1.stdout)
        p3 = subprocess.Popen(shlex.split(cmd3), stdout=subprocess.PIPE, stdin=p2.stdout)
        p4 = subprocess.Popen(shlex.split(cmd4), stdout=genome_file, stdin=p3.stdout)
    p4.communicate()


def build_coverage_files(tumor_bam=None, normal_bam=None, genome_path=None,
                         tumor_cov_path=None, normal_cov_path=None,
                         target_bed_path=None):
    """Build normal and tumor coverage files for whole exome.

    bedtools coverage -g $g -d -sorted -a $CODING_REGIONS -b normal.bam > cov_normal.bed;
    bedtools coverage -g $g -d -sorted -a $CODING_REGIONS -b tumor.bam > cov_tumor.bed;
    """

    config = load_config()
    bedtools_path = config.get('paths', 'BEDTOOLS', 'bedtools')

    cmd_template = "{bedtools_path} coverage -g {g} -d -sorted -a {target} -b {bam};"

    for (bam, out_path) in [(tumor_bam, tumor_cov_path), (normal_bam, normal_cov_path)]:
        print("Generating coverage for {}".format(bam))
        with open(out_path, 'w') as out:
            cmd = cmd_template.format(bedtools_path=bedtools_path, g=genome_path,
                                      target=target_bed_path, bam=bam)
            print("...bedtools command: {}".format(cmd))
            proc = subprocess.Popen(shlex.split(cmd), stdin=subprocess.DEVNULL, stdout=out)
            proc.communicate()
