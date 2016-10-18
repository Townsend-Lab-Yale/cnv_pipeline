import os
import subprocess
import shlex

from .config import load_config


def trim_vcf(vcf_in=None, tumor_id=None, normal_id=None,
             ratio_min=0.4, ratio_max=0.6, min_depth_n=10, min_depth_t=20,
             min_gq_n=90, vcf_out=None, sample_dir=None):
    """Create new filtered vcf file for tumor and normal sample."""
    if normal_id is None:
        normal_id = tumor_id + 'N'
    if vcf_out is None:
        vcf_out = os.path.join(sample_dir, 'snps_trimmed.vcf')
    config = load_config()
    gatk = config.get('gatk', 'cmd', fallback='gatk')
    fasta = config.get('paths', 'FASTA')
    cmd = (
        """{gatk} -T SelectVariants -R {fasta} -V {vcf} -sn {normal} -sn {tumor} """
        """-selectType SNP --restrictAllelesTo BIALLELIC -select 'vc.getGenotype("{normal}").isHet() """
        """&& vc.getGenotype("{normal}").getGQ() > {gq} """
        """&& 1.0 * vc.getGenotype("{normal}").getAD().1 /  vc.getGenotype("{normal}").getDP() > {ratio_min} """
        """&&  1.0 * vc.getGenotype("{normal}").getAD().1 /  vc.getGenotype("{normal}").getDP() < {ratio_max} """
        """&& vc.getGenotype("{normal}").getDP() > {min_depth_n} """
        """&& vc.getGenotype("{tumor}").getDP() > {min_depth_t}'""")
    cmd = cmd.format(gatk=gatk, vcf=vcf_in, fasta=fasta, normal=normal_id, tumor=tumor_id,
                     gq=min_gq_n, ratio_min=ratio_min, ratio_max=ratio_max,
                     min_depth_n=min_depth_n, min_depth_t=min_depth_t)
    print("VCF trim command: {}".format(cmd))
    with open(vcf_out, 'w') as new_vcf:
        subprocess.check_call(shlex.split(cmd), stdout=new_vcf)
    print("Created vcf: {}".format(vcf_out))
