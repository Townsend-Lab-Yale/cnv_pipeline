import os
import sys
import subprocess

import feather
import pandas as pd

script_dir = os.path.dirname(os.path.realpath(__file__))

def baf_from_vcf(patient_id, vcf_path, baf_path):
    """
    Args:
        patient_id (str): used for saving saasCNV-style snp data to feather.
    ."""
    initial_dir = os.getcwd()
    dest_dir = os.path.dirname(vcf_path)
    os.chdir(dest_dir)

    table_basename = '{}_baf'.format(patient_id)
    feather_path = os.path.join(dest_dir, table_basename + '.feather')

    # RUN Rscript
    print("Running R script.")
    rscript_path = os.path.join(script_dir, 'vcf2table.R')
    # R: Parse vcf, write data to feather
    subprocess.check_call(['Rscript', rscript_path, vcf_path, table_basename])

    # Finalize dataframe
    print("Finalizing dataframe.")
    df = feather.read_dataframe(feather_path)
    df.rename(columns={'CHROM': 'chrom', 'POS': 'SNP_loc'}, inplace=True)
    df['control_doc'] = df['Normal.REF.DP'] + df['Normal.ALT.DP']
    df['tumor_doc'] = df['Tumor.REF.DP'] + df['Tumor.ALT.DP']
    df['control_BAF'] = df['Normal.ALT.DP'] / df['control_doc']
    df['tumor_BAF'] = df['Tumor.ALT.DP'] / df['tumor_doc']
    # mbaf_fun = lambda x: x if x >= 0.5 else 1-x
    # df.control_BAF = df.control_BAF.apply(mbaf_fun)
    # df.tumor_BAF = df.tumor_BAF.apply(mbaf_fun)
    df = df[['chrom', 'SNP_loc', 'control_BAF', 'tumor_BAF', 'control_doc', 'tumor_doc']]
    # df.chrom = df.chrom.apply(lambda x: 'chr' + str(x) if x != 'MT' else 'chrM')
    df.dropna(inplace=True)  # Drop null entries
    df.to_csv(baf_path, sep='\t', index=False)
    print("Saved BAF data to file: {}".format(baf_path))

if __name__ == '__main__':
    patient_id = sys.argv[1]
    vcf_path = sys.argv[2]
    baf_path = sys.argv[3]
    baf_from_vcf(patient_id, vcf_path, baf_path)
