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

fs_INL_phylum = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10_phylum.csv', index_col=0, header=0)
fs_INL_phylum['Taxonomic level'] = 'Phylum'
fs_INL_phylum['Lake'] = 'Inland'
fs_INL_class = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10_class.csv', index_col=0, header=0)
fs_INL_class['Taxonomic level'] = 'Class'
fs_INL_class['Lake'] = 'Inland'
fs_INL_order = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10_order.csv', index_col=0, header=0)
fs_INL_order['Taxonomic level'] = 'Order'
fs_INL_order['Lake'] = 'Inland'
fs_INL_family = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10_family.csv', index_col=0, header=0)
fs_INL_family['Taxonomic level'] = 'Family'
fs_INL_family['Lake'] = 'Inland'
fs_INL_genus = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10_genus.csv', index_col=0, header=0)
fs_INL_genus['Taxonomic level'] = 'Genus'
fs_INL_genus['Lake'] = 'Inland'
fs_INL_species = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10_species.csv', index_col=0, header=0)
fs_INL_species['Taxonomic level'] = 'Species'
fs_INL_species['Lake'] = 'Inland'
fs_INL_OTU = pd.read_csv('FS_final/Inland_fs_scores_LNA_5seq10.csv', index_col=0, header=0)
fs_INL_OTU['Taxonomic level'] = 'OTU'
fs_INL_OTU['Lake'] = 'Inland'
fs_MICH_phylum = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10_phylum.csv', index_col=0, header=0)
fs_MICH_phylum['Taxonomic level'] = 'Phylum'
fs_MICH_phylum['Lake'] = 'Michigan'
fs_MICH_class = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10_class.csv', index_col=0, header=0)
fs_MICH_class['Taxonomic level'] = 'Class'
fs_MICH_class['Lake'] = 'Michigan'
fs_MICH_order = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10_order.csv', index_col=0, header=0)
fs_MICH_order['Taxonomic level'] = 'Order'
fs_MICH_order['Lake'] = 'Michigan'
fs_MICH_family = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10_family.csv', index_col=0, header=0)
fs_MICH_family['Taxonomic level'] = 'Family'
fs_MICH_family['Lake'] = 'Michigan'
fs_MICH_genus = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10_genus.csv', index_col=0, header=0)
fs_MICH_genus['Taxonomic level'] = 'Genus'
fs_MICH_genus['Lake'] = 'Michigan'
fs_MICH_species = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10_species.csv', index_col=0, header=0)
fs_MICH_species['Taxonomic level'] = 'Species'
fs_MICH_species['Lake'] = 'Michigan'
fs_MICH_OTU = pd.read_csv('FS_final/Michigan_fs_scores_LNA_5seq10.csv', index_col=0, header=0)
fs_MICH_OTU['Taxonomic level'] = 'OTU'
fs_MICH_OTU['Lake'] = 'Michigan'
fs_MUS_phylum = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_phylum.csv', index_col=0, header=0)
fs_MUS_phylum['Taxonomic level'] = 'Phylum'
fs_MUS_phylum['Lake'] = 'Muskegon'
fs_MUS_class = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_class.csv', index_col=0, header=0)
fs_MUS_class['Taxonomic level'] = 'Class'
fs_MUS_class['Lake'] = 'Muskegon'
fs_MUS_order = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_order.csv', index_col=0, header=0)
fs_MUS_order['Taxonomic level'] = 'Order'
fs_MUS_order['Lake'] = 'Muskegon'
fs_MUS_family = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_family.csv', index_col=0, header=0)
fs_MUS_family['Taxonomic level'] = 'Family'
fs_MUS_family['Lake'] = 'Muskegon'
fs_MUS_genus = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_genus.csv', index_col=0, header=0)
fs_MUS_genus['Taxonomic level'] = 'Genus'
fs_MUS_genus['Lake'] = 'Muskegon'
fs_MUS_species = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10_species.csv', index_col=0, header=0)
fs_MUS_species['Taxonomic level'] = 'Species'
fs_MUS_species['Lake'] = 'Muskegon'
fs_MUS_OTU = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10.csv', index_col=0, header=0)
fs_MUS_OTU['Taxonomic level'] = 'OTU'
fs_MUS_OTU['Lake'] = 'Muskegon'

def patch_violinplot():
    """Patch seaborn's violinplot in current axis to workaround matplotlib's bug ##5423."""
    from matplotlib.collections import PolyCollection
    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, PolyCollection):
            art.set_edgecolor((0.1, 0.1, 0.1)) 
            
def plot_R2CV_HNA_RL_TAX(df): 
    g = sns.factorplot(x='Taxonomic level',y='RL score',data=df, col='Lake', kind='box', size=6, aspect=0.8, sharey=True)
    plt.subplots_adjust(top=0.88)
    g.fig.suptitle('Distribution RL score LNAcc ', size=22)
    g.set_titles(size=19)
    g.set_axis_labels('Taxonomic level','Score')
    g.set_xlabels(fontsize=18)
    g.set_ylabels(fontsize=18)
    plt.savefig('Analysis_Figures/Box_LNA_RLscore_TAX_5seq10.png',bbox_inches='tight', dpi=300)
    plt.show()
    return df

fs = pd.concat([fs_INL_phylum,fs_INL_class,fs_INL_order,fs_INL_family,fs_INL_genus,fs_INL_species,fs_INL_OTU,fs_MICH_phylum,fs_MICH_class,fs_MICH_order,fs_MICH_family,fs_MICH_genus,fs_MICH_species,fs_MICH_OTU,fs_MUS_phylum,fs_MUS_class,fs_MUS_order,fs_MUS_family,fs_MUS_genus,fs_MUS_species,fs_MUS_OTU],ignore_index=True,axis=0)
df = plot_R2CV_HNA_RL_TAX(fs)