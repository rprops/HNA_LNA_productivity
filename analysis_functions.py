# -*- coding: utf-8 -*-
"""
Spyder Editor
"""

import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.cross_validation import train_test_split 
from sklearn.ensemble import RandomForestRegressor
#from sklearn.feature_selection import f_regression
#from sklearn.feature_selection import mutual_info_regression
#from sklearn import cross_validation
#from sklearn import metrics
#from sklearn.decomposition import PCA
from sklearn.linear_model import Lasso
from sklearn.linear_model import lasso_stability_path
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_predict
#from sklearn.model_selection import cross_val_score
from sklearn.model_selection import validation_curve
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
#from sklearn.preprocessing import LabelEncoder
#from sklearn.preprocessing import OneHotEncoder
#from sklearn.preprocessing import StandardScaler
#import matplotlib.pyplot as plt
#import seaborn as sns

import time
start_time = time.time()

def preprocess_df(df, threshold, ABS): 
    if ABS == True:  
        df_sum = df.sum(axis=1)
        df_relabs = df.divide(df_sum, axis='rows')
        return df.ix[:,(df_relabs > threshold).any(axis=0)]
    else: 
        return df.ix[:,(df > threshold).any(axis=0)]

'''Import data'''
#data_abs = pd.read_csv('data/Chloroplasts_removed/nochloro_absolute_otu.tsv', sep=' ', index_col=None, header=0)
#data_rel = pd.read_csv('data/Chloroplasts_removed/nochloro_relative_otu.tsv', sep=' ', index_col=None, header=0)
#target = pd.read_csv('data/Chloroplasts_removed/nochloro_HNA_LNA.tsv', sep=' ', index_col=0, header=0)
#productivity = pd.read_csv('data/Chloroplasts_removed/productivity_data.tsv', sep=' ', index_col=0, header=0)


''' Preprocessing '''
#data_abs = preprocess_df(data_abs, abun, True)
#data_rel = preprocess_df(data_rel, abun, False)
#features = list(data_rel.columns)
#y = target.loc[:,'HNA.cells']

''' Normalize data '''
#data_abs = pd.DataFrame(scaler.fit_transform(data_abs), columns=features)
#data_rel = pd.DataFrame(scaler.fit_transform(data_rel), columns=features)
#y = pd.Series(scaler.fit_transform(y), index=target.index)

#abs_corr = data_abs.corr('pearson')
#rel_corr = data_rel.corr('pearson')
#mask = np.zeros_like(abs_corr, dtype=np.bool)
#mask[np.triu_indices_from(mask)] = False
     
#f, ax = plt.subplots(figsize=(14, 14))
#cmap = sns.diverging_palette(220, 10, as_cmap=True)
#sns.heatmap(rel_corr, mask=mask, cmap=cmap, vmin=-1, vmax=1, square=True, annot=False, xticklabels=True, yticklabels=True, linewidths=0., cbar_kws={"shrink": .5}, ax=ax,  annot_kws={"size": 16})     
#ax.tick_params(labelsize=0)
#plt.savefig('pearson_rel_all.png', bbox_inches='tight')
#plt.show()

def get_train_test(df, target): 
    x_train, x_test, y_train, y_test = train_test_split(df, target, test_size = 0.3, random_state = 26)
    return x_train, x_test, y_train, y_test

''' 2. lassoCV '''        
def perform_lassoCV(df, target):
    lassoCV = linear_model.LassoCV(eps=0.001, n_alphas=400, max_iter=20000, cv = 5, normalize = False, random_state=6)
    lassoCV.fit(df, target)
    return lassoCV, lassoCV.mse_path_, lassoCV.alpha_
    
def get_lassoCV_scores(df, target, alpha): 
    lasso = Lasso(alpha=alpha, normalize=False, max_iter=20000)
    return cross_val_predict(lasso, df, target, cv=5)
    
def perform_randomizedLasso(df, target, alpha): 
    randomLasso = linear_model.RandomizedLasso(alpha=alpha, sample_fraction=0.75, n_resampling=500, normalize=False, random_state=33)
    randomLasso.fit(df, target)
    return randomLasso.scores_
    
def get_RF(): 
    return RandomForestRegressor(n_estimators=500, criterion='mse')#, max_features=0.5)  
    
def get_lassoCV(): 
    lassoCV =  linear_model.LassoCV(eps=0.0001, n_alphas=400, max_iter=20000, cv=5, normalize=False)
    return lassoCV

def get_lasso10CV(): 
    lassoCV = linear_model.LassoCV(eps=0.0001, n_alphas=200, max_iter=200000, cv=10, normalize=False)
    return lassoCV
    
def get_ridgeCV(): 
    return linear_model.RidgeCV(alphas=(0.0001,1,400), cv=5, normalize=False)  

def get_ridge10CV(): 
    return linear_model.RidgeCV(alphas=(0.0001,1,400), cv=10, normalize=False)   
    
def perform_nested_RF_cv(df, target): 
    rf = get_RF()
    clf = GridSearchCV(rf, param_grid={'max_features': [0.05,0.1,0.3,0.5,0.7,0.9,0.95]}, cv=5)
    nested_pred = cross_val_predict(clf, df, target, cv=4)
    return nested_pred

def perform_nested_RF_loocv(df, target): 
    rf = get_RF()
    clf = GridSearchCV(rf, param_grid={'max_features': [0.05,0.1,0.3,0.5,0.7,0.9,0.95]}, cv=10)
    nested_pred = cross_val_predict(clf, df, target, cv=df.shape[0])
    return nested_pred
      
def perform_nested_lasso_cv(df, target): 
    lassoCV = get_lassoCV()
    outer_cv = KFold(n_splits=4, shuffle=False)
    alphas = np.zeros(4)
    t = 0
    pred = pd.Series(index=df.index)
    for idx_train, idx_test in outer_cv.split(df, target): 
        lassoCV.fit(df.iloc[idx_train,:], target[idx_train])
        #lassoCV.predict(df.iloc[idx_test,:], target[idx_test])
        alphas[t] = lassoCV.alpha_
        pred.iloc[idx_test] = lassoCV.predict(df.iloc[idx_test,:])
        t+=1
    return alphas, pred  

def perform_nested_lasso_loocv(df, target): 
    lassoCV = get_lasso10CV()
    outer_cv = KFold(n_splits=df.shape[0], shuffle=False)
    alphas = np.zeros(df.shape[0])
    t = 0
    pred = pd.Series(index=df.index)
    for idx_train, idx_test in outer_cv.split(df, target): 
        lassoCV.fit(df.iloc[idx_train,:], target[idx_train])
        #lassoCV.predict(df.iloc[idx_test,:], target[idx_test])
        alphas[t] = lassoCV.alpha_
        pred.iloc[idx_test] = lassoCV.predict(df.iloc[idx_test,:])
        t+=1
    return alphas, pred 

def perform_nested_randomized_lasso_cv(df, target, alpha): 
    outer_cv = KFold(n_splits=4, shuffle=False)
    scores = pd.DataFrame(columns=list(df.columns))
    t = 0
    for idx_train, idx_test in outer_cv.split(df, target): 
        scores.loc[t] = perform_randomizedLasso(df.iloc[idx_train,:], target[idx_train], alpha)
        t+=1
    return scores
        
def perform_lasso_stability_path(df, target): 
    return  lasso_stability_path(df, target, scaling=0.5, random_state=2703, n_resampling=1000, n_grid=300, sample_fraction=0.75)
    
def perform_nested_ridge_cv(df, target): 
    ridgeCV = get_ridgeCV()
    outer_cv = KFold(n_splits=4, shuffle=False, random_state=2)
    alphas = np.zeros(4)
    t = 0
    pred = pd.Series(index=df.index)
    for idx_train, idx_test in outer_cv.split(df, target): 
        ridgeCV.fit(df.iloc[idx_train,:], target[idx_train])
        #lassoCV.predict(df.iloc[idx_test,:], target[idx_test])
        alphas[t] = ridgeCV.alpha_
        pred.iloc[idx_test] = ridgeCV.predict(df.iloc[idx_test,:])
        t+=1
    return alphas, pred  

def perform_nested_ridge_loocv(df, target): 
    ridgeCV = get_ridge10CV()
    outer_cv = KFold(n_splits=df.shape[0], shuffle=False, random_state=5)
    alphas = np.zeros(df.shape[0])
    t = 0
    pred = pd.Series(index=df.index)
    for idx_train, idx_test in outer_cv.split(df, target): 
        ridgeCV.fit(df.iloc[idx_train,:], target[idx_train])
        #lassoCV.predict(df.iloc[idx_test,:], target[idx_test])
        alphas[t] = ridgeCV.alpha_
        pred.iloc[idx_test] = ridgeCV.predict(df.iloc[idx_test,:])
        t+=1
    return alphas, pred 
    
def perform_validation_curve(df, target): 
    rf = get_RF()
    param_range = [0.05,0.1,0.2,0.4,0.5,0.6,0.8,0.95]
    #clf = GridSearchCV(rf, param_grid={'max_features': [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1]}, cv=10)
    return validation_curve(rf, df, target, 'max_features', param_range, cv=10)
    
def get_r2(true, pred): 
    return r2_score(true, pred)

def get_r2_adj(r2, n, k): 
    return 1. - ((1.-r2)*(n-1.)/(n-k-1.))
    
def get_mse(true, pred): 
    return mean_squared_error(true, pred)

''' Univariate feature selection methods '''
'''
mi = pd.Series(mutual_info_regression(data_abs, y, discrete_features=False, n_neighbors=3), index=features)
f, p = f_regression(data_abs, y)
f = pd.Series(f, index=features)
p = pd.Series(p, index=features)
idx = list(p[p < 0.001].index.values)
'''

''' Get nested estimate of R-squared '''
'''
alphas, nested_preds = perform_nested_lasso_cv(data_abs[features], y)
r2 = get_r2(y, nested_preds)
print(r2)
plt.figure(figsize=(18,12))
#plt.plot(rf.feature_importances_)
plt.scatter(y, nested_preds)
plt.xlabel('HNA.cells True', size=24)
plt.ylabel('HNA.cells Predicted', size=24)
plt.title('r2: ' + str(r2), size=28)
plt.show()

#plt.savefig('analysis_RF_abs.png')
'''

''' Lasso analysis '''
'''
#preds = perform_nested_cv(data_rel[features], y)
alpha = alphas.mean()
scores = pd.Series(perform_randomizedLasso(data_abs[features], y, alpha), index=features)
scores_nested = pd.DataFrame(perform_nested_randomized_lasso_cv(data_abs[features], y, alpha))
scores_nested_mean = scores_nested.mean()
scores.sort_values(ascending=False,inplace=True)
pd.DataFrame(scores).to_excel('scores_nochloro_OTUs.xlsx')
pd.DataFrame(scores_nested_mean).to_excel('scoresnested_nochloro_OTUs.xlsx')
scores_nested_mean.sort_values(ascending=False,inplace=True)
scores = scores[scores.values > 0.24]
scores_nested_mean = scores_nested_mean[scores_nested_mean > 0.16]
features_new = scores.index
features_mean_new = scores_nested_mean.index
print(len(features_mean_new))
alphas, preds = perform_nested_lasso_cv(data_abs[features_mean_new],y)
r2_new = get_r2(y,preds)

#scores = perform_nested_randomized_lasso_cv(data_abs[features], y, alpha)
#scores_mean = scores.mean()
#scores_std = scores.std()
#scores_relstd = scores_std/scores_mean


i=0
r2_new_scores = np.zeros(51)
for thresh in np.arange(0.,1.02,0.02): 
    if len(scores_nested_mean) > 0: 
        scores = scores[scores > thresh]
        features_new = scores.index
        alphas, preds = perform_nested_ridge_cv(data_abs[features_new],y)
        r2_new_scores[i] = get_r2(y,preds)
        i+=1
    else: 
        break
    
plt.subplots(figsize=(18,12))
plt.axis([0,1,0.,1])
sns.pointplot(np.arange(0,1.02,0.02),r2_new_scores, linestyles='--')
plt.xlabel('Threshold score', size=22)
plt.ylabel(r'$R^2$', size=22)
plt.savefig('thresh_vs_r2_randomizedlasso.png',bbox_inches='tight')
plt.show()
'''

''' PCA analysis '''
'''
pca = PCA()
data_abs_pca = pd.DataFrame(pca.fit_transform(data_abs[features_new]), index=data_abs.index)
var = pca.explained_variance_ratio_
plt.figure(figsize=(18,12))
plt.scatter(data_abs_pca.loc[idx_inland,1].values, data_abs_pca.loc[idx_inland,2].values, c='r')
plt.scatter(data_abs_pca.loc[idx_michigan,1].values, data_abs_pca.loc[idx_michigan,2].values, c='g')
plt.scatter(data_abs_pca.loc[idx_muskegon,1].values, data_abs_pca.loc[idx_muskegon,2].values, c='b')
plt.show()
'''

'''
plt.figure(figsize=(18,12))
#plt.plot(rf.feature_importances_)
plt.scatter(y, preds)
plt.xlabel('HNA.cells True', size=24)
plt.ylabel('HNA.cells Predicted', size=24)
plt.title('r2: ' + str(r2_new), size=28)
#plt.show()
plt.savefig('Lasso_selectedOTUsrandomLasso.png')
'''

''' Randomized Lasso analysis '''
'''
alpha_grid, scores_path = perform_lasso_stability_path(data_abs[features].values, y)
plt.figure()
plt.xlabel(r'${\alpha/\alpha_{max}^0.333$')
plt.ylabel('Stability')
plt.plot((alpha_grid[1:])**.333, scores_path.T[1:], 'k')
plt.show()
scores_path = pd.DataFrame(scores_path, index=features)
'''

''' Split up according to lake '''
'''
idx_inland = target.loc[target.Lake == 'Inland'].index
idx_michigan = target[target.Lake == 'Michigan'].index
idx_muskegon = target[target.Lake == 'Muskegon'].index
data_abs_inland = data_abs.loc[idx_inland,:]
data_abs_michigan = data_abs.loc[idx_michigan,:]
data_abs_muskegon = data_abs.loc[idx_muskegon,:]
y_inland = target.loc[idx_inland,'HNA.cells']#/target['Total.cells']
y_michigan = target.loc[idx_michigan,'HNA.cells']#/target['Total.cells']
y_muskegon = target.loc[idx_muskegon,'HNA.cells']#/target['Total.cells']

alphas1, preds1 = perform_nested_lasso_cv(data_abs_inland[features],y_inland)
alphas2, preds2 = perform_nested_lasso_cv(data_abs_michigan[features],y_michigan)
alphas3, preds3 = perform_nested_lasso_cv(data_abs_muskegon[features],y_muskegon)
alpha1 = alphas1.mean()
alpha2 = alphas2.mean()
alpha3 = alphas3.mean()

y = pd.concat([y_inland,y_michigan,y_muskegon])
preds = pd.concat([preds1,preds2,preds3])

r2 = get_r2(y, preds)

scores1 = perform_randomizedLasso(data_abs_inland[features],y_inland,alpha1)
scores1 = pd.Series(perform_randomizedLasso(data_abs_inland[features],y_inland,alpha1),index=features)
scores2 = perform_randomizedLasso(data_abs_michigan[features],y_michigan,alpha2)
scores2 = pd.Series(perform_randomizedLasso(data_abs_michigan[features],y_michigan,alpha2),index=features)
scores3 = perform_randomizedLasso(data_abs_muskegon[features],y_muskegon,alpha3)
scores3 = pd.Series(perform_randomizedLasso(data_abs_muskegon[features],y_muskegon,alpha3),index=features)

'''

print("--- %s seconds ---" % (time.time() - start_time)) 