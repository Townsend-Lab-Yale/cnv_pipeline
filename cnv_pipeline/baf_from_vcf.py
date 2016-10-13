import os
import subprocess

import feather
# import pandas as pd

script_dir = os.path.dirname(os.path.realpath(__file__))


def baf_from_vcf(vcf_path, baf_path, feather_path=None, col_tumor=11, col_normal=10,
                 format_dp_index=5, mq_cutoff=30,
                 chroms='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'):
    """
    Args:
        patient_id (str): used for saving saasCNV-style snp data to feather.
        vcf_path (str): full or relative vcf path
        col_{tumor,normal} (int): index of {tumor,normal} column in vcf. 1-based.
        format_dp_index (int): index of DP in vcf FORMAT column. 1-based.
        MQ_cutoff (int): cutoff for MQ, in INFO column.

    Intermediate files:
        <baf_path>.feather: created by vcf2table.R
    """
    if feather_path is None:
        feather_path = baf_path + '.feather'

    # RUN Rscript
    print("Running vcf to baf conversion R script.")
    rscript_path = os.path.join(script_dir, 'vcf2table.R')
    # R: Parse vcf, write data to feather
    subprocess.check_call(['Rscript', rscript_path, vcf_path, feather_path,
                           str(col_tumor), str(col_normal), str(format_dp_index), chroms, str(mq_cutoff)])
    # Import results
    df = feather.read_dataframe(feather_path)

    # Finalize dataframe
    print("Finalizing dataframe.")
    df.rename(columns={'CHROM': 'chrom', 'POS': 'SNP_loc'}, inplace=True)
    df['control_doc'] = df['Normal.REF.DP'] + df['Normal.ALT.DP']
    df['tumor_doc'] = df['Tumor.REF.DP'] + df['Tumor.ALT.DP']
    df['control_BAF'] = df['Normal.ALT.DP'] / df['control_doc']
    df['tumor_BAF'] = df['Tumor.ALT.DP'] / df['tumor_doc']
    # mbaf_fun = lambda x: x if x >= 0.5 else 1-x
    # df.control_BAF = df.control_BAF.apply(mbaf_fun)
    # df.tumor_BAF = df.tumor_BAF.apply(mbaf_fun)
    df = df[['chrom', 'SNP_loc', 'control_BAF', 'tumor_BAF', 'control_doc', 'tumor_doc']].copy()
    # df.chrom = df.chrom.apply(lambda x: 'chr' + str(x) if x != 'MT' else 'chrM')
    df.dropna(inplace=True)  # Drop null entries
    df.to_csv(baf_path, sep='\t', index=False)
    print("Saved BAF data to file: {}".format(baf_path))


# if __name__ == '__main__':
#     import sys
#     patient_id = sys.argv[1]
#     vcf_path = sys.argv[2]
#     baf_path = sys.argv[3]
#     baf_from_vcf(patient_id, vcf_path, baf_path)
