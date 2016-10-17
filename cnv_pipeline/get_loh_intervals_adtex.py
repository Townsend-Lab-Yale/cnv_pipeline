"""Parse ADTEx results, identify LOH-SNP CNV intervals, trim and plot."""

import os
import shlex
import subprocess

import pandas as pd

from cnv_pipeline.plot_chr_axis import plot_chr_axis, plot_chr_intervals

chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
          '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']


def finalize_loh(proj_dir):
    zygosity_df, loh_segs = prep_loh_dataframes(proj_dir)
    loh_segs = trim_loh_intervals(zygosity_df, loh_segs, out_dir=proj_dir)  # min_ratio=0.8
    plot_loh(loh_segs, zygosity_df, out_dir=proj_dir)


def prep_loh_dataframes(proj_dir):
    """Builds preliminary LOH dataframes from ADTEx results, ready for trimming.

    loh_snps.bed: contains SNPs marked as loh in proj_dir/zygosity/zygosity.res
    cnv_interval.bed: all CNV segments, extracted from proj_dir/cnv.result

    Args:
        proj_dir (str): adtex output dir. Includes cnv.result and 'zygosity' dir.

    Returns:
        z (pd.DataFrame): raw zygosity dataframe, from ADTEx.
        loh_segs (pd.DataFrame): cnv intervals containing LOH SNPs.
    """
    print("Building preliminary loh dataframes in {}".format(proj_dir))
    res_path = os.path.join(proj_dir, 'zygosity', 'zygosity.res')
    cnv_path = os.path.join(proj_dir, 'cnv.result')
    # TEMPORARY FILES
    snp_out = os.path.join(proj_dir, 'temp_loh_snps.bed')
    cnv_out = os.path.join(proj_dir, 'temp_cnv_interval.bed')
    loh_bed = os.path.join(proj_dir, 'temp_loh_intervals.bed')

    z = pd.read_csv(res_path, sep='\t', dtype={'chrom': str})
    z.chrom = z.chrom.astype('category', categories=chroms, ordered=True)
    # SAVE SNP BED FILE
    z_loh = z[z.zygosity == 'LOH']
    loh_snp_bed = pd.concat([z_loh.chrom, z_loh.SNP_loc - 1, z_loh.SNP_loc], axis=1)
    loh_snp_bed.to_csv(snp_out, index=False, sep='\t', header=False)
    # SAVE CNV_INTERVAL BED FILE
    cnv = pd.read_csv(cnv_path, sep='\t', 
                      dtype={'chr': str})
    cnv.rename(columns={'chr': 'chrom'}, inplace=True)
    interval_bed = pd.concat([cnv.chrom, cnv.CNV_start - 1, cnv.CNV_end], axis=1).drop_duplicates()
    interval_bed.to_csv(cnv_out, index=False, sep='\t', header=False)

    # BEDTOOLS INTERSECTION
    """intersectBed -b loh_snps.bed -a cnv_interval.bed -wa -u  > loh_intervals.bed;"""
    print("Running bedtools intersection")
    cmd = "intersectBed -b {snps} -a {cnv} -wa -u".format(snps=snp_out, cnv=cnv_out)
    with open(loh_bed, 'w') as out:
        proc = subprocess.Popen(shlex.split(cmd), stdin=subprocess.DEVNULL, stdout=out)
        proc.communicate()

    loh_segs = pd.read_csv(loh_bed, sep='\t', header=None, 
                           names=['chrom', 'pos_start', 'pos_end'], dtype={'chrom': str})
    loh_segs.chrom = loh_segs.chrom.astype('category', categories=chroms, ordered=True)
    loh_segs['orig_start'] = loh_segs['pos_start']
    loh_segs['orig_end'] = loh_segs['pos_end']
    
    # Remove temporary files
    for i in [snp_out, cnv_out, loh_bed]:
        os.remove(i)

    return z, loh_segs


def trim_loh_intervals(z, loh_segs, out_dir=None, min_ratio=0.8):
    """Trims/filters LOH segments based on breadth of LOH SNPs.

    Args (via prep_loh_dataframes):
        z (pd.DataFrame): raw zygosity dataframe, from ADTEx.
        loh_segs (pd.DataFrame): cnv intervals containing LOH SNPs. 0-based start coord.

    Returns:
        loh_segs (df): stripped loh segs. start pos is 0-based.
    """
    
    print("Trimming loh intervals.")
    if out_dir is None:
        out_dir = os.getcwd()
    
    drop_path = os.path.join(out_dir, 'loh_intervals_dropped.bed')
    final_path = os.path.join(out_dir, 'loh_intervals_final.bed')

    z.sort_values(['chrom', 'SNP_loc'], inplace=True)  # Ensure z is sorted
    for ind, seg in loh_segs.iterrows():
        try:
            l_a = z[(z.chrom == seg.chrom) & (z.SNP_loc > seg.pos_start) & (z.zygosity == 'LOH')].iloc[0]
        except IndexError:
            break
        l_b = z[(z.chrom == seg.chrom) & (z.SNP_loc <= seg.pos_end) & (z.zygosity == 'LOH')].iloc[-1]
        l_d = l_b.SNP_loc - l_a.SNP_loc
        s_d = seg.pos_end - seg.pos_start - 1
        seg_ratio = l_d / s_d if s_d != 0 else 0
        # delete seg if only one LOH point
        if l_a.name == l_b.name:
            loh_segs.loc[seg.name, ['pos_start', 'pos_end']] = (pd.np.nan, pd.np.nan)
            continue  # !!!
        # update seg if smaller than ratio
        if seg_ratio < min_ratio:
            loh_segs.loc[seg.name, ['pos_start', 'pos_end']] = (l_a.SNP_loc - 1, l_b.SNP_loc)

    # dropped
    loh_segs.loc[loh_segs.isnull().any(axis=1), ['chrom', 'orig_start', 'orig_end']]\
        .to_csv(drop_path, sep='\t', index=False, header=False)

    # Export final intervals
    # loh_segs.pos_start = loh_segs.pos_start + 1  # convert to 1-based coords
    loh_segs.dropna(inplace=True)
    loh_segs.pos_start = loh_segs.pos_start.astype(int)
    loh_segs.pos_end = loh_segs.pos_end.astype(int)
    loh_segs[['chrom', 'pos_start', 'pos_end']].to_csv(final_path, sep='\t', index=False, header=False)
    return loh_segs


def plot_loh(loh, z, out_dir=None, y_var='tumor_BAF'):
    """.
    Args:
        y_var (str): 'tumor_BAF' or 'mirrored_BAF'
    """
    print('Plotting LOH.')
    if out_dir is None:
        out_dir = os.getcwd()

    use_Y = (z.chrom == 'Y').any()
    use_MT = (z.chrom == 'MT').any()

    out_path = os.path.join(out_dir, 'LOH_plot.png')

    loh.pos_start = loh.pos_start + 1  # convert from BED format
    
    # Plot raw data using zygosity calls for colors
    ax = plot_chr_axis('chrom', 'SNP_loc', y_var, data=z[z.zygosity == 'LOH'], figsize=(19,3),
                       color='r', ylim=(0,1), use_Y=use_Y, use_MT=use_MT)
    plot_chr_axis('chrom', 'SNP_loc', y_var, data=z[z.zygosity == 'ASCNA'],
                  color='g', format_axis=False, ax=ax, use_Y=use_Y, use_MT=use_MT)
    plot_chr_axis('chrom', 'SNP_loc', y_var, data=z[z.zygosity == 'HET'],
                  color='b', format_axis=False, ax=ax, use_Y=use_Y, use_MT=use_MT)

    # Draw rectangles for intervals
    plot_chr_intervals(ax, 'chrom', 'pos_start', 'pos_end', data=loh, use_Y=use_Y, use_MT=use_MT)

    ax.set_xlabel('')
    ax.set_ylabel('Tumor BAF')
    hf = ax.get_figure()
    hf.tight_layout()
    hf.savefig(out_path)


if __name__ == '__main__':
    import sys
    proj_dir = sys.argv[1]
    finalize_loh(proj_dir)
