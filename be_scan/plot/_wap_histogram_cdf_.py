"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196-3D-Clustering-v7.py Created on Mon Jun  1 19:03:10 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF

#%% plot_wap_hist_cdf() - plots histogram and ECDF for iterwaps output

def plot_wap_hist_cdf(df_waps, cond='wap', sns_context='talk', sns_palette='deep',
                      show_kde=False, show_pval=True, figtitle='Summed PWES Histogram and ECDF',
                      out_prefix='', save_pdf=True, save_png=True, return_plot=True):

    # Plotting parameters
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']

    df_plot = df_waps.loc[df_waps['iter_num'] != 'my_wap']
    my_wap = float(df_waps.loc[df_waps['iter_num'] == 'my_wap'][cond])
    # set up style parameters for plotting in seaborn/mpl
    sns.set_context(sns_context)
    if sns_palette is not None:
        sns.set_palette(sns_palette)
    # generate histogram of wap scores w/ sns distplot
    fig, (ax1,ax2) = plt.subplots(ncols=2, figsize=(10,4))
    sns.distplot(a=df_plot[cond], bins=50, kde=show_kde, hist_kws={'alpha':0.8}, ax=ax1)
    # plot actual WAP as vertical line
    ax1.axvline(x=my_wap, c='black', ls='--')
    ax1.text(x=0.05, y=0.9, s='WAP = {}'.format(str(round(my_wap))),
             transform=ax1.transAxes)
    # aesthetics
    for edge,spine in ax1.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    ax1.set_xlabel('Summed PWES')
    ax1.set_ylabel('Frequency')
    # generate ecdf and find pval for my_wap; if p = 1, find inequality
    wap_ecdf = ECDF(df_plot['wap'], side='left')
    # find pval for my_wap; if p = 1, find inequality
    if wap_ecdf(my_wap) == 1:
        pval = float(1 / df_plot.shape[0])
        str_p = 'p < ' + str(pval)
    else:
        pval = float(1 - wap_ecdf(my_wap))
        str_p = 'p = ' + str(pval)
    # generate cdf plot with sns lineplot
    sns.lineplot(x=wap_ecdf.x, y=wap_ecdf.y, ax=ax2)
    # aesthetics
    for edge,spine in ax2.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    ax2.set_xlabel('Summed PWES')
    ax2.set_ylabel('')
    fig.suptitle(figtitle + ' (%s)' % str_p)

    if save_pdf:
        plt.savefig(out_prefix + '_WAP_analysis.pdf', format='pdf', bbox_inches='tight')
    if save_png:
        plt.savefig(out_prefix + '_WAP_analysis.png', format='png', bbox_inches='tight')
    if return_plot:
        return fig
    else:
        plt.close()
        return
