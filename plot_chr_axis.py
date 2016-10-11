import os

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
import matplotlib.patches as patches


class GenomeInfo():
    """Loads basic genome attributes required for plotting genomic data."""
    def __init__(self, use_Y=False, use_MT=False):
        """Create attributes: size_df, genome_size, start_dict, lines."""
        chr_len = ([('1', 249250621), ('2', 243199373), ('3', 198022430), ('4', 191154276),
                  ('5', 180915260), ('6', 171115067), ('7', 159138663), ('8', 146364022),
                  ('9', 141213431), ('10', 135534747), ('11', 135006516), ('12', 133851895),
                  ('13', 115169878), ('14', 107349540), ('15', 102531392), ('16', 90354753),
                  ('17', 81195210), ('18', 78077248), ('19', 59128983), ('20', 63025520),
                  ('21', 48129895), ('22', 51304566), ('X', 155270560), ('Y', 59373566), ('MT', 16569)])
        p_len_dict = {'1': 125000000, '2': 93300000, '3': 91000000, '4': 50400000, '5': 48400000, 
            '6': 61000000, '7': 59900000, '8': 45600000, '9': 49000000, '10': 40200000, 
            '11': 53700000, '12': 35800000, '13': 17900000, '14': 17600000, '15': 19000000, 
            '16': 36600000, '17': 24000000, '18': 17200000, '19': 26500000, '20': 27500000, 
            '21': 13200000, '22': 14700000, 'X': 60600000, 'Y': 12500000}
        sizes = pd.DataFrame(chr_len, columns=['chrom', 'n_sites'])
        sizes['increment'] = pd.Series([0] + list(sizes.n_sites.cumsum()[:-1]))
        sizes['start'] = sizes['increment'] + 1
        mids = []
        for i in range(len(sizes)):
            try:
                mids.append(0.5 * (sizes.start.iloc[i] + sizes.start.iloc[i+1]))
            except IndexError:
                mids.append(sizes.start.iloc[i] + 0.5 * sizes.n_sites.iloc[i-1])
        sizes['mids'] = mids
        sizes['centro'] = sizes.apply(lambda r: p_len_dict.get(r.chrom, round(r.n_sites/2)) + r.increment, axis=1)        
        sizes.set_index('chrom', inplace=True)
        if not use_Y:
            sizes.drop('Y', inplace=True)
        if not use_MT:
            sizes.drop('MT', inplace=True)
        self.p_len_dict = p_len_dict
        self.size_df = sizes
        self.genome_size = sizes.iloc[-1].increment + sizes.iloc[-1].n_sites
        self.start_dict = dict(sizes.start.iteritems())
        self.lines = [((x, 0), (x, 0.5)) for x in sizes.start]  # chromosome boundaries

    def get_genome_pos(self, chrom, pos):
        return pos + self.start_dict[str(chrom)]


def plot_chr_axis(chrom, pos, y, data=None, ax=None, figsize=None, use_Y=False, use_MT=False, 
                  alpha=0.5, markersize=2, ylim=None, format_axis=True, **dict_kw):
    """Plot with genome x-axis."""
    if data is not None:
        chrom = data[chrom]
        pos = data[pos]
        y = data[y]
    g = GenomeInfo(use_Y=use_Y, use_MT=use_MT)
    d = pd.concat([pd.Series(chrom), pd.Series(pos), pd.Series(y)], axis=1)  # position dataframe
    d.columns = ['chrom', 'pos', 'y']
    d['g_pos'] = d.apply(lambda r: g.get_genome_pos(r.chrom, r.pos), axis=1)
    if figsize is None:
        figsize = (14,3)
    if ax is None:
        hfig, ax = plt.subplots(1, figsize=figsize)
    d.plot(x='g_pos', y='y', kind='line', marker='.',
        linewidth=0, alpha=alpha, markersize=markersize, ax=ax, legend=False, **dict_kw)
    ymin, ymax = ax.get_ylim() if ylim is None else ylim
    if format_axis:
        g.lines = [((x, ymin), (x, ymax)) for x in g.size_df.start]  # chromosome boundaries
        lc = mc.LineCollection(g.lines, linewidths=0.5, colors='gray')  # chr boundaries
        ax.add_collection(lc)
        ax.xaxis.set_ticks_position('none') 
        ax.set_xticks(g.size_df.mids)
        ax.set_xticklabels(list(g.size_df.index.values))
        ax.grid('off', axis='x')
        ax.grid('on', which='minor', axis='y')
    ax.set_ylim([ymin, ymax])
    return ax


def plot_chr_intervals(ax, chrom, posA, posB, data=None, facecolor='r', alpha=0.1, ylim=None, **plot_kw):
    """Plot intervals on pre-existing genome axis."""
    if data is not None:
        chrom = data[chrom]
        posA = data[posA]
        posB = data[posB]
    g = GenomeInfo()
    d = pd.concat([pd.Series(chrom), pd.Series(posA), pd.Series(posB)], axis=1)  # position dataframe
    d.columns = ['chrom', 'posA', 'posB']
    d['g_posA'] = d.apply(lambda r: g.get_genome_pos(r.chrom, r.posA), axis=1)
    d['width'] = d.posB - d.posA + 1
    ymin, ymax = ax.get_ylim()
    for ind, seg in d.iterrows():
        xs = seg.g_posA
        xd = seg.width
        ax.add_patch(patches.Rectangle((xs, ymin), xd, ymax, alpha=alpha, facecolor=facecolor, edgecolor='none'))
