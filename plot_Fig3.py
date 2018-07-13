#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:47:32 2018

@author: prubbens
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_color_codes()
tips = sns.load_dataset("tips")

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
sns.set_style("ticks")

fs_phylum_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_phylum.csv', index_col=0, header=0)
fs_class_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_class.csv', index_col=0, header=0)
fs_order_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_order.csv', index_col=0, header=0)
fs_family_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_family.csv', index_col=0, header=0)
fs_genus_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_genus.csv', index_col=0, header=0)
fs_species_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_species.csv', index_col=0, header=0)
fs_OTU_hna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10.csv', index_col=0, header=0)

mfe_lassoRL_phylum_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking_phylum.csv', index_col=0, header=0)
mfe_lassoRL_phylum_hna['Taxonomic rank'] = 'Phylum'
mfe_lassoRL_class_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking_class.csv', index_col=0, header=0)
mfe_lassoRL_class_hna['Taxonomic rank'] = 'Class'
mfe_lassoRL_order_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking_order.csv', index_col=0, header=0)
mfe_lassoRL_order_hna['Taxonomic rank'] = 'Order'
mfe_lassoRL_family_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking_family.csv', index_col=0, header=0)
mfe_lassoRL_family_hna['Taxonomic rank'] = 'Family'
mfe_lassoRL_genus_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking_genus.csv', index_col=0, header=0)
mfe_lassoRL_genus_hna['Taxonomic rank'] = 'Genus'
mfe_lassoRL_species_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking_species.csv', index_col=0, header=0)
mfe_lassoRL_species_hna['Taxonomic rank'] = 'Species'
mfe_lassoRL_OTU_hna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_HNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_OTU_hna['Taxonomic rank'] = 'OTU'

fs_phylum_lna = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_phylum.csv', index_col=0, header=0)
fs_class_lna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_class_LNA.csv', index_col=0, header=0)
fs_order_lna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_order_LNA.csv', index_col=0, header=0)
fs_family_lna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_family_LNA.csv', index_col=0, header=0)
fs_genus_lna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_genus_LNA.csv', index_col=0, header=0)
fs_species_lna = pd.read_csv('FS_final/Muskegon_fs_scores_HNA_5seq10_species_LNA.csv', index_col=0, header=0)
fs_OTU_lna = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10.csv', index_col=0, header=0)

mfe_lassoRL_phylum_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking_phylum.csv', index_col=0, header=0)
mfe_lassoRL_phylum_lna['Taxonomic rank'] = 'Phylum'
mfe_lassoRL_class_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking_class.csv', index_col=0, header=0)
mfe_lassoRL_class_lna['Taxonomic rank'] = 'Class'
mfe_lassoRL_order_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking_order.csv', index_col=0, header=0)
mfe_lassoRL_order_lna['Taxonomic rank'] = 'Order'
mfe_lassoRL_family_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking_family.csv', index_col=0, header=0)
mfe_lassoRL_family_lna['Taxonomic rank'] = 'Family'
mfe_lassoRL_genus_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking_genus.csv', index_col=0, header=0)
mfe_lassoRL_genus_lna['Taxonomic rank'] = 'Genus'
mfe_lassoRL_species_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking_species.csv', index_col=0, header=0)
mfe_lassoRL_species_lna['Taxonomic rank'] = 'Species'
mfe_lassoRL_OTU_lna = pd.read_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking.csv', index_col=0, header=0)
mfe_lassoRL_OTU_lna['Taxonomic rank'] = 'OTU'

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
            
def plot_R2CV_HNA_RL_TAX(df): 
    df.loc[:,'Number of taxa'] = df.loc[:,'Number of taxa'].astype(int)
    df = pd.melt(df,id_vars=['Number of taxa','Taxonomic rank','Target'], value_vars=['R2_CV'], var_name='Functional group', value_name='R2')
    g = sns.lmplot(x='Number of taxa',y='R2',data=df, hue='Target', col='Taxonomic rank', col_wrap=4, fit_reg=False, sharey=True, sharex=False, legend=False)
    plt.subplots_adjust(top=0.9)
    g.set_xlabels(fontsize=20)
    g.set_ylabels(fontsize=20)
    col_order = ['A','B','C','D','E','F','G']
    for ax, title in zip(g.axes.flat, col_order):
        ax.text(-0.0, 1.0, title, fontsize=16, weight='bold')
    x_coord = [11,31,21,39,34,33,102] #HNAcc
    y_coord = [0.540,0.666,0.6,0.746,0.655,0.732,0.853] #HNAcc
    for ax, x, y in zip(g.axes.flat, x_coord, y_coord):
        ax.plot([x,x], [-0.2,y], linewidth=2, linestyle='-', c='royalblue')
    x_coord = [8,25,10,39,21,35,102] #LNAcc
    y_coord = [0.648,0.716,0.766,0.903,0.878,0.898,0.911] #LNAcc
    for ax, x, y in zip(g.axes.flat, x_coord, y_coord):
        ax.plot([x,x], [-0.2,y], linewidth=2, linestyle='--', c='orange')
    g.set_titles(size=19)
    for ax in g.axes.flat: 
        ax.set_xticklabels(labels=ax.get_xticklabels(), size=14)
        #ax.set_yticklabels(labels=np.arange(-0.2,1,0.2), size=14)
    g.set_axis_labels('Number of taxa',r'$R^2_{CV}$')
    plt.legend(loc='lower right')
    plt.savefig('Analysis_Figures/R2CV_MUS_Lasso_RL_TAX_5seq10.png',bbox_inches='tight', dpi=300)
    plt.show()
    return df

mfe_lassoRL_phylum_hna = calc_otus_ranking(fs_phylum_hna,mfe_lassoRL_phylum_hna,'RL ranking','RL ranking')
mfe_lassoRL_phylum_hna['Target'] = 'HNAcc'
mfe_lassoRL_class_hna = calc_otus_ranking(fs_class_hna,mfe_lassoRL_class_hna,'RL ranking','RL ranking')
mfe_lassoRL_class_hna['Target'] = 'HNAcc'
mfe_lassoRL_order_hna = calc_otus_ranking(fs_order_hna,mfe_lassoRL_order_hna,'RL ranking','RL ranking')
mfe_lassoRL_order_hna['Target'] = 'HNAcc'
mfe_lassoRL_family_hna = calc_otus_ranking(fs_family_hna,mfe_lassoRL_family_hna,'RL ranking','RL ranking')
mfe_lassoRL_family_hna['Target'] = 'HNAcc'
mfe_lassoRL_genus_hna = calc_otus_ranking(fs_genus_hna,mfe_lassoRL_genus_hna,'RL ranking','RL ranking')
mfe_lassoRL_genus_hna['Target'] = 'HNAcc'
mfe_lassoRL_species_hna = calc_otus_ranking(fs_species_hna,mfe_lassoRL_species_hna,'RL ranking','RL ranking')
mfe_lassoRL_species_hna['Target'] = 'HNAcc'
mfe_lassoRL_OTU_hna = calc_otus_ranking(fs_OTU_hna,mfe_lassoRL_OTU_hna,'RL ranking','RL ranking')
mfe_lassoRL_OTU_hna['Target'] = 'HNAcc'
mfe_hna = pd.concat([mfe_lassoRL_phylum_hna,mfe_lassoRL_class_hna,mfe_lassoRL_order_hna,mfe_lassoRL_family_hna,mfe_lassoRL_genus_hna,mfe_lassoRL_species_hna,mfe_lassoRL_OTU_hna])
mfe_hna = mfe_hna.rename(columns={'RL HNA': 'R2_CV'})

mfe_lassoRL_phylum_lna = calc_otus_ranking(fs_phylum_lna,mfe_lassoRL_phylum_lna,'RL ranking','RL ranking LNA')
mfe_lassoRL_phylum_lna['Target'] = 'LNAcc'
mfe_lassoRL_class_lna = calc_otus_ranking(fs_class_lna,mfe_lassoRL_class_lna,'RL ranking LNA','RL ranking LNA')
mfe_lassoRL_class_lna['Target'] = 'LNAcc'
mfe_lassoRL_order_lna = calc_otus_ranking(fs_order_lna,mfe_lassoRL_order_lna,'RL ranking LNA','RL ranking LNA')
mfe_lassoRL_order_lna['Target'] = 'LNAcc'
mfe_lassoRL_family_lna = calc_otus_ranking(fs_family_lna,mfe_lassoRL_family_lna,'RL ranking LNA','RL ranking LNA')
mfe_lassoRL_family_lna['Target'] = 'LNAcc'
mfe_lassoRL_genus_lna = calc_otus_ranking(fs_genus_lna,mfe_lassoRL_genus_lna,'RL ranking LNA','RL ranking LNA')
mfe_lassoRL_genus_lna['Target'] = 'LNAcc'
mfe_lassoRL_species_lna = calc_otus_ranking(fs_species_lna,mfe_lassoRL_species_lna,'RL ranking LNA','RL ranking LNA')
mfe_lassoRL_species_lna['Target'] = 'LNAcc'
mfe_lassoRL_OTU_lna = calc_otus_ranking(fs_OTU_lna,mfe_lassoRL_OTU_lna,'RL ranking','RL ranking LNA')
mfe_lassoRL_OTU_lna['Target'] = 'LNAcc'
mfe_lna = pd.concat([mfe_lassoRL_phylum_lna,mfe_lassoRL_class_lna,mfe_lassoRL_order_lna,mfe_lassoRL_family_lna,mfe_lassoRL_genus_lna,mfe_lassoRL_species_lna,mfe_lassoRL_OTU_lna])
mfe_lna = mfe_lna.rename(columns={'RL LNA': 'R2_CV'})

mfe = pd.concat([mfe_hna,mfe_lna])

df = plot_R2CV_HNA_RL_TAX(mfe)