#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:25:48 2018

@author: prubbens
"""

'''Import packages'''
'''Requires numpy, pandas, scikit-learn, and matplotlib/seaborn'''
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import LabelEncoder

'''Import script which contains functions'''
from analysis_functions import standardize_df, get_r2_scores_Lasso

data_rel = pd.read_csv('data/Chloroplasts_removed/ByLake_Filtering/5in10/muskegon/muskegon_relative_otu_5in10.tsv', sep=' ', index_col=None, header=0, float_precision='high')
target = pd.read_csv('ByLake_Filtering/5in10/muskegon/muskegon_sampledata_5in10.tsv', sep= ' ', index_col=0, header=0)
index = target.Lake[target.Lake == 'Muskegon'].index

data_rel = data_rel.loc[index,:]
target = target.loc[index,:]

#Create target columns of HNA-values: 
hna = target.loc[index,'HNA.cells']
hna_rel = hna/target.loc[index,'Total.cells']
hna = pd.Series(hna, index=hna.index)
hna_rel = pd.Series(hna_rel, index=hna.index)

#Create target columns of LNA-values: 
lna = target.loc[index,'LNA.cells']
lna_rel = lna/target.loc[index,'Total.cells']
lna = pd.Series(lna, index=lna.index)
lna_rel = pd.Series(lna_rel, index=lna.index)

otus = list(data_rel.columns)

data_stand, scaler = standardize_df(data_rel,otus)

fs = pd.read_csv('FS_final/Muskegon_fs_scores_LNA_5seq10.csv', index_col=0, header=0)

target.loc[index,'spatiotemporal'] = target.loc[index,'Year'].astype(str) + target.loc[index,'Site']
le = LabelEncoder()
le_values = le.fit_transform(target.loc[index,'spatiotemporal'].values)
logo = LeaveOneGroupOut().split(data_stand, groups=le_values)
cv=logo

unique_scores_hna_RL_, r2_cv_hna_RL_ = get_r2_scores_Lasso(data_stand.loc[index,otus], hna, otus, fs.loc[otus,'RL ranking'], cv, le_values)
unique_scores_lna_RL_, r2_cv_lna_RL_ = get_r2_scores_Lasso(data_stand.loc[index,otus], lna, otus, fs.loc[otus,'RL ranking'], cv, le_values)

final_RL = pd.DataFrame({'RL ranking': unique_scores_hna_RL_, 'RL HNA': r2_cv_hna_RL_, 'RL LNA': r2_cv_lna_RL_})
#final_RL.to_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_RLRanking.csv')

target.loc[index,'spatiotemporal'] = target.loc[index,'Year'].astype(str) + target.loc[index,'Site']
le = LabelEncoder()
le_values = le.fit_transform(target.loc[index,'spatiotemporal'].values)
logo = LeaveOneGroupOut().split(data_stand, groups=le_values)
cv=logo

unique_scores_hna_, r2_cv_hna_ = get_r2_scores_Lasso(data_stand.loc[index,otus], hna, otus, fs.loc[otus,'Boruta ranking'], cv, le_values)
unique_scores_lna_, r2_cv_lna_ = get_r2_scores_Lasso(data_stand.loc[index,otus], lna, otus, fs.loc[otus,'Boruta ranking'], cv, le_values)

final_Boruta = pd.DataFrame({'Boruta ranking': unique_scores_hna_, 'Boruta HNA': r2_cv_hna_, 'Boruta LNA': r2_cv_lna_})
#final_Boruta.to_csv('RFE/R2CV_RFE_Lasso_MUS_5seq10_LNA_BorutaRanking.csv')