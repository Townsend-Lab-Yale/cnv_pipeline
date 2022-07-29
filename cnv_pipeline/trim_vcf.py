import os
import subprocess
import shlex

from .config import GATK_ALIAS


def trim_vcf(vcf_in=None, tumor_id=None, normal_id=None, ref_fasta=None,
             ratio_min=0.4, ratio_max=0.6, min_depth_n=10, min_depth_t=20,
             min_gq_n=90, vcf_out=None, sample_dir=None):
    """Create new filtered vcf file for tumor and normal sample."""
    if normal_id is None:
        normal_id = tumor_id + 'N'
    if vcf_out is None:
        vcf_out = os.path.join(sample_dir, 'snps_trimmed.vcf')
    cmd = (
        f"""{GATK_ALIAS} SelectVariants -R {ref_fasta} -V {vcf_in} -sn {normal_id} -sn {tumor_id} """
        f"""--select-type-to-include SNP --restrict-alleles-to BIALLELIC -select 'vc.getGenotype("{normal_id}").isHet() """
        f"""&& vc.getGenotype("{normal_id}").getDP() > {min_depth_n - 1} """
        f"""&& vc.getGenotype("{tumor_id}").getDP() > {min_depth_t - 1} """
        f"""&& vc.getGenotype("{normal_id}").getGQ() > {min_gq_n} """
        f"""&& 1.0 * vc.getGenotype("{normal_id}").getAD().1 /  vc.getGenotype("{normal_id}").getDP() > {ratio_min} """
        f"""&& 1.0 * vc.getGenotype("{normal_id}").getAD().1 /  vc.getGenotype("{normal_id}").getDP() < {ratio_max}' """ 
        f"""--output {vcf_out}""")
    print("VCF trim command: {}".format(cmd))
    subprocess.check_call(shlex.split(cmd))
    print("Created vcf: {}".format(vcf_out))
