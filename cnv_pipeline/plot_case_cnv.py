import os
os.environ['QT_QPA_PLATFORM']='offscreen'

import cnv_pipeline.plot_chr_axis as pcnv

import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')


g = pcnv.GenomeInfo()

AX_DICT=dict(xticks=g.size_df.mids,
             xticklabels=list(g.size_df.index.values),
             yticks=[],
             xlim=[0,g.genome_size])


def plot_case_cnv(case_id, tumor_ids=None, feather_dict=None, cnv_dict=None, dim='baf'):
    """Plot cnv profile from SAAS-CNV for single case, return fig/axis handles.

    Args:
        tumor_ids (iterable): tumor ids.
        feather_dict [dict]: dictionary mapping tumor ids to cnv data.
        dim (str): 'baf' or 'lrd' for minor b-allele frequency or log2 ratio of depths, respectively
    Returns:
        hf: matplotlib figure handle
        axs: matplotlib axis handles.
    """

    hf, axs = plt.subplots(len(tumor_ids), 1, figsize=(14,10), sharex=True,
    #                        gridspec_kw={'height_ratios':[10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]},
                           subplot_kw=AX_DICT)
    # for ax in axs:
    #     format_axis(ax, g)

    if dim == 'baf':
        dim_str = 'minor allele frequencies'
        ylim = (0, 1)
    elif dim == 'lrd':
        dim_str = 'T:N depth ratio (log2)'
        ylim = (-4, 4)
    else:
        raise Exception("Invalid dim parameter ({}). Must specify 'baf' or 'lrd' as dim parameter.".format(dim))

    for sample_id, ax in zip(tumor_ids, axs):
        temp = pd.read_table(cnv_dict[sample_id], dtype={'chr': str})
        temp['chr'] = temp['chr'].str.lstrip('chr')
        temp_loss = temp[temp.CNV.isin(['loss', 'LOH'])]
        temp_gain = temp[temp.CNV.isin(['gain'])]

        baf = pd.read_feather(feather_dict[sample_id])
        baf['baf'] = baf['Tumor.ALT.DP'] / (baf['Tumor.REF.DP'] + baf['Tumor.ALT.DP'])
        baf['lrd'] = pd.np.log2((baf['Tumor.ALT.DP'] + baf['Tumor.REF.DP']) / (baf['Normal.ALT.DP'] + baf['Normal.REF.DP']))

        pcnv.plot_chr_axis(chrom='CHROM', pos='POS', y=dim, data=baf, markersize=1, alpha=0.2, ylim=ylim, ax=ax)
        if len(temp_loss):
            pcnv.plot_chr_intervals(chrom='chr', posA='posStart', posB='posEnd', data=temp_loss, facecolor='b', ax=ax)
        if len(temp_gain):
            pcnv.plot_chr_intervals(chrom='chr', posA='posStart', posB='posEnd', data=temp_gain, facecolor='r', ax=ax)
        ax.set_ylabel(sample_id)
    ax.set_xlabel('Chromosome')

    hf = ax.get_figure()
    # hf.savefig('scna.pdf')
    axs[0].set_title('SCNA results ({})'.format(dim_str))
    return hf, axs


def plot_case_cnv_samples_table(sample_file='samples.txt'):
    """Plot cnv figures for each case. Assumes default paths for saas-cnv output.

    Args:
        s (pd.DataFrame): rows are {patient_id,case_id}, {tumor_id,sample_id}, ...
    """
    s = pd.read_table(sample_file, dtype={'patient_id': str, 'case_id': str, 'tumor_id': str, 'sample_id': str})
    s.rename(columns={'patient_id': 'case_id', 'tumor_id': 'sample_id'}, inplace=True)

    sample_lists = s.groupby('case_id')['sample_id'].apply(lambda s: sorted(list(s.values)))

    for case_id, sample_list in sample_lists.iteritems():
        # LOAD DATA (CNV INPUT AND OUTPUT)
        feather_dict = {s: '{}/saas.feather'.format(s) for s in sample_list}
        cnv_dict = {s: '{}/saasCNV_results/mid_res/seq.cnv.txt'.format(s) for s in sample_list}
        for dim in ['lrd', 'baf']:
            hf, axs = plot_case_cnv(case_id, tumor_ids=sample_list, feather_dict=feather_dict, cnv_dict=cnv_dict, dim=dim)
            hf.savefig('{}_{}.png'.format(case_id, dim), bbox_inches='tight')
