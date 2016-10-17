import os
import re
import subprocess

import pandas as pd
import feather
# import pandas as pd

script_dir = os.path.dirname(os.path.realpath(__file__))


def get_vcf_properties(vcf_path, tumor_id=None, normal_id=None):
    """Locate tumor and normal columns. Identify AD index within FORMAT."""
    skip = 0
    with open(vcf_path, 'r') as file:
        for line in file:
            if line.startswith('##'):
                skip += 1
            else:
                break
    df = pd.read_csv(vcf_path, sep='\t', skiprows=skip, nrows=100)
    col_tumor = list(df.columns).index(tumor_id) + 1  # 1-based
    col_normal = list(df.columns).index(normal_id) + 1  # 1-based
    # get 'AD' location
    # index_ad = df.FORMAT.str.split(':').apply(lambda l: l.index('AD')).unique()
    # if len(index_ad) > 1:
    #     raise Exception('Multiple AD locations within format string')
    # index_ad = index_ad[0] + 1  # 1-based
    return col_tumor, col_normal, skip


def get_mq(info):
    v = info.split(';')
    vals = [i for i in v if re.match('^MQ=[0-9\.]+$', i)]
    if len(vals) != 1:
        return pd.np.nan
    mq = float(vals[0].split('=')[1])
    return mq


def parse_format(f, col_gt=0, col_ad=1):
    """0/1:61,217:278:99:5815,0,1176"""
    vals = f.split(':')
    n_ref, n_alt = vals[col_ad].split(',')
    return vals[col_gt], int(n_ref), int(n_alt)


def baf_from_vcf(vcf_path, baf_path, feather_path=None, tumor_id=None, normal_id=None,
                 mq_cutoff=30, chroms='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'):
    """Args:
        patient_id (str): used for saving saasCNV-style snp data to feather.
        vcf_path (str): full or relative vcf path
        col_{tumor,normal} (int): index of {tumor,normal} column in vcf. 1-based.
        format_ad_index (int): index of AD in vcf FORMAT column. 1-based.
        MQ_cutoff (int): cutoff for MQ, in INFO column.

    Intermediate files:
        <baf_path>.feather: created by vcf2table.R
    """
    chrom_list = chroms.split(',')
    if feather_path is None:
        feather_path = baf_path + '.feather'
    v = get_vcf_properties(vcf_path=vcf_path, tumor_id=tumor_id,
                           normal_id=normal_id)
    col_tumor, col_normal, skip = v
    df = pd.read_csv(vcf_path, sep='\t', skiprows=skip)
    df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
    df['MQ'] = df.INFO.apply(lambda info: get_mq(info))
    for col, names in [(col_tumor, ['Tumor.GT', 'Tumor.REF.DP', 'Tumor.ALT.DP']),
                       (col_normal, ['Normal.GT', 'Normal.REF.DP', 'Normal.ALT.DP'])]:
        vals = df.ix[:, col].apply(lambda f: parse_format(f))
        df[names[0]], df[names[1]], df[names[2]] = zip(*vals)
        df[names[1]] = df[names[1]].astype(pd.np.int64)
        df[names[2]] = df[names[2]].astype(pd.np.int64)
    df = df.loc[(df.MQ > mq_cutoff) & (df.CHROM.isin(chrom_list)),
                ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "MQ",
                 "Normal.GT", "Normal.REF.DP", "Normal.ALT.DP", "Tumor.GT",
                 "Tumor.REF.DP", "Tumor.ALT.DP"]].copy()
    feather.write_dataframe(df, feather_path)  # for saasCNV

    # Finalize baf file
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


def baf_from_vcf_R(vcf_path, baf_path, feather_path=None, tumor_id=None, normal_id=None,
                 mq_cutoff=30,
                 chroms='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'):
    """
    Args:
        patient_id (str): used for saving saasCNV-style snp data to feather.
        vcf_path (str): full or relative vcf path
        col_{tumor,normal} (int): index of {tumor,normal} column in vcf. 1-based.
        format_ad_index (int): index of AD in vcf FORMAT column. 1-based.
        MQ_cutoff (int): cutoff for MQ, in INFO column.

    Intermediate files:
        <baf_path>.feather: created by vcf2table.R
    """
    if feather_path is None:
        feather_path = baf_path + '.feather'

    if os.path.exists(baf_path) and os.path.exists(feather_path):
        print("BAF files exist. Skipping vcf conversion.")
        return

    col_tumor, col_normal, index_ad = get_vcf_properties(vcf_path, tumor_id=tumor_id, normal_id=normal_id)

    if not os.path.exists(feather_path):
        # RUN Rscript
        print("Running vcf to baf conversion R script.")
        rscript_path = os.path.join(script_dir, 'vcf2table.R')
        # R: Parse vcf, write data to feather
        subprocess.check_call(['Rscript', rscript_path, vcf_path, feather_path,
                               str(col_tumor), str(col_normal), str(index_ad), chroms, str(mq_cutoff)])
    else:
        print("Loading pre-existing baf feather file ({}).".format(feather_path))
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
