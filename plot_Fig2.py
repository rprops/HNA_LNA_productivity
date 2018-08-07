#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:47:32 2018

@author: prubbens
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator

#plt.style.use('ggplot')
from os import listdir
import scipy as sc
from statsmodels.graphics.boxplots import violinplot
import seaborn as sns
sns.set_color_codes()
tips = sns.load_dataset("tips")

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
sns.set_style("ticks")

fs_INL_hna = pd.read_csv('FS_final/Inland_fs_scores_HNA_5seq10.csv', index_col=0, header=0)
fs_MICH_hna = pd.read_csv('FS_final/Michigan_fs_scores_HNA_5seq10.csv', index_col=0, header=0)
fs_MUS_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10.csv', index_col=0, header=0)
fs_INL_lna = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10.csv', index_col=0, header=0)
fs_MICH_lna = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10.csv', index_col=0, header=0)
fs_MUS_lna = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10.csv', index_col=0, header=0)

mfe_lassoRL_INL_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_INL_5seq10_HNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_INL_hna['Lake'] = 'Inland'
mfe_lassoRL_INL_hna['Variable Selection'] = 'Randomized Lasso'
mfe_lassoRL_INL_hna['Target'] = 'HNAcc'
mfe_lassoRL_MICH_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MICH_5seq10_HNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_MICH_hna['Lake'] = 'Michigan'
mfe_lassoRL_MICH_hna['Variable Selection'] = 'Randomized Lasso'
mfe_lassoRL_MICH_hna['Target'] = 'HNAcc'
mfe_lassoRL_MUS_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_MUS_hna['Lake'] = 'Muskegon'
mfe_lassoRL_MUS_hna['Variable Selection'] = 'Randomized Lasso'
mfe_lassoRL_MUS_hna['Target'] = 'HNAcc'
mfe_lassoRL_INL_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_INL_5seq10_LNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_INL_lna['Lake'] = 'Inland'
mfe_lassoRL_INL_lna['Variable Selection'] = 'Randomized Lasso'
mfe_lassoRL_INL_lna['Target'] = 'LNAcc'
mfe_lassoRL_MICH_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MICH_5seq10_LNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_MICH_lna['Lake'] = 'Michigan'
mfe_lassoRL_MICH_lna['Variable Selection'] = 'Randomized Lasso'
mfe_lassoRL_MICH_lna['Target'] = 'LNAcc'
mfe_lassoRL_MUS_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_MUS_lna['Lake'] = 'Muskegon'
mfe_lassoRL_MUS_lna['Variable Selection'] = 'Randomized Lasso'
mfe_lassoRL_MUS_lna['Target'] = 'LNAcc'

def patch_violinplot():
    """Patch seaborn's violinplot in current axis to workaround matplotlib's bug ##5423."""
    from matplotlib.collections import PolyCollection
    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, PolyCollection):
            art.set_edgecolor((0.1, 0.1, 0.1))
            
def calc_otus_ranking(fs_rank,mfe,fs_header,mfe_header): 
    ranking = mfe.loc[:,mfe_header]
    n_otus = np.zeros(len(ranking))
    i=0
    for rank in ranking:
       n_otus[i] = fs_rank[fs_rank.loc[:,fs_header] <= rank].shape[0]
       i+=1
    mfe.loc[:,'Number of taxa'] = n_otus
    return mfe     
            

def plot_R2CV_Lake(df): 
    df.loc[:,'Number of taxa'] = df.loc[:,'Number of taxa'].astype(int)
    df = pd.melt(df,id_vars=['Number of taxa','Lake','Target'], value_vars=['R2_CV'], var_name='Functional group', value_name='R2')
    g = sns.lmplot(x='Number of taxa',y='R2',data=df, hue='Target', col='Lake', fit_reg=False, sharex=False, legend=False, scatter_kws=dict(edgecolor="k", linewidth=1))
    col_order = ['A','B','C']
    for ax, title in zip(g.axes.flat, col_order):
        ax.text(-25.0, 0.98, title, fontsize=18, weight='bold')
    x_coord = [56,13,102]
    y_coord = [0.92,0.52,0.85]
    for ax, x, y in zip(g.axes.flat, x_coord, y_coord):
        ax.plot([x,x], [-0.2,y], linewidth=2, linestyle='-', c='royalblue')
    x_coord = [84,10,102]
    y_coord = [0.87,0.79,0.91]
    for ax, x, y in zip(g.axes.flat, x_coord, y_coord):
        ax.plot([x,x], [-0.2,y], linewidth=2, linestyle='--', c='orange')
    plt.subplots_adjust(top=0.8)
    plt.legend(loc='lower right')
    g.set_titles(size=20)
    g.set_axis_labels('Number of taxa',r'$R^2_{CV}$')
    g.set_xlabels(fontsize=18)
    g.set_ylabels(fontsize=18)
    plt.savefig('Analysis_Figures/R2CV_HNA_LNA_Lasso_RL.png',bbox_inches='tight', dpi=500)
    plt.show()
    return df

mfe_lassoRL_INL_hna = calc_otus_ranking(fs_INL_hna,mfe_lassoRL_INL_hna,'RL ranking','RL ranking')
mfe_lassoRL_MICH_hna = calc_otus_ranking(fs_MICH_hna,mfe_lassoRL_MICH_hna,'RL ranking','RL ranking')
mfe_lassoRL_MUS_hna = calc_otus_ranking(fs_MUS_hna,mfe_lassoRL_MUS_hna,'RL ranking','RL ranking')
mfe_lassoRL_INL_lna = calc_otus_ranking(fs_INL_lna,mfe_lassoRL_INL_lna,'RL ranking','RL ranking')
mfe_lassoRL_MICH_lna = calc_otus_ranking(fs_MICH_lna,mfe_lassoRL_MICH_lna,'RL ranking','RL ranking')
mfe_lassoRL_MUS_lna = calc_otus_ranking(fs_MUS_lna,mfe_lassoRL_MUS_lna,'RL ranking','RL ranking LNA')

mfe_hna = pd.concat([mfe_lassoRL_INL_hna,mfe_lassoRL_MICH_hna,mfe_lassoRL_MUS_hna])
mfe_hna = mfe_hna.rename(columns={'RL HNA': 'R2_CV'})
mfe_lna = pd.concat([mfe_lassoRL_INL_lna,mfe_lassoRL_MICH_lna,mfe_lassoRL_MUS_lna])
mfe_lna = mfe_lna.rename(columns={'RL LNA': 'R2_CV'})
mfe = pd.concat([mfe_hna,mfe_lna],axis=0)
df = plot_R2CV_Lake(mfe)