# -*- coding: utf-8 -*-
"""
Spyder Editor
"""

import numpy as np
import pandas as pd

from sklearn import linear_model
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso
from sklearn.linear_model import lasso_stability_path
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_predict, RandomizedSearchCV
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler

from boruta_py.boruta.boruta_py import BorutaPy

def preprocess_df(df, threshold, ABS): 
    if ABS == True:  
        df_sum = df.sum(axis=1)
        df_relabs = df.divide(df_sum, axis='rows')
        return df.ix[:,(df_relabs > threshold).any(axis=0)]
    else: 
        return df.ix[:,(df > threshold).any(axis=0)]
    
def preprocess_df_meanabun(df, threshold, ABS): 
    if ABS == True:  
        df_sum = df.sum(axis=1)
        df_relabs = df.divide(df_sum, axis='rows')
        return df.ix[:,(df_relabs.mean() > threshold)]
    else: 
        return df.ix[:,(df.mean() > threshold)]

def get_train_test(df, target): 
    x_train, x_test, y_train, y_test = train_test_split(df, target, test_size = 0.3, random_state = 26)
    return x_train, x_test, y_train, y_test

''' 2. lassoCV '''        
def perform_lassoCV(df, target, cv):
    lassoCV = linear_model.LassoCV(eps=0.0001, n_alphas=200, max_iter=10000, cv = cv, normalize = False, random_state=6)
    lassoCV.fit(df, target)
    return lassoCV, lassoCV.mse_path_, lassoCV.alpha_
    
def get_lassoCV_scores(df, target, alpha): 
    lasso = Lasso(alpha=alpha, normalize=False, max_iter=20000)
    return cross_val_predict(lasso, df, target, cv=5)
    
def perform_randomizedLasso(df, target): 
    randomLasso = linear_model.RandomizedLasso(alpha=np.logspace(-3,3,100), sample_fraction=0.5, n_resampling=500, normalize=False, random_state=36, scaling=0.5)
    randomLasso.fit(df, target)
    return randomLasso.scores_#, randomLasso.all_scores_
      
    
def get_lassoCV(cv): 
    lassoCV = linear_model.LassoCV(eps=0.0001, n_alphas=400, max_iter=200000, cv=cv, normalize=False, random_state=9)
    return lassoCV

        
def perform_lasso_stability_path(df, target): 
    return  lasso_stability_path(df, target, scaling=0.5, random_state=2703, n_resampling=1000, n_grid=300, sample_fraction=0.75)
    
def get_r2(true, pred): 
    return r2_score(true, pred)

def get_r2_adj(r2, n, k): 
    return 1. - ((1.-r2)*(n-1.)/(n-k-1.))
    
def get_mse(true, pred): 
    return mean_squared_error(true, pred)

def standardize_df(df,features): 
    scaler = StandardScaler()
    scaler.fit(df.loc[:,features])
    df_stand = pd.DataFrame(scaler.transform(df.loc[:,features]),index=df.index,columns=features)    
    return df_stand, scaler    

def get_lassoCV_alpha(df,target,features,cv): 
    lassoCV = get_lassoCV(cv)
    lassoCV.fit(df.loc[:,features],target)
    #mse = np.sum(lassoCV.mse_path_, axis=1)
    return lassoCV.alpha_#, lassoCV.alphas_.max()

def get_lassoCV_alpha_max(df,target,features,cv): 
    lassoCV = get_lassoCV(cv)
    lassoCV.fit(df.loc[:,features],target)#, groups)
    #mse = np.sum(lassoCV.mse_path_, axis=1)
    return lassoCV.alpha_, lassoCV.alphas_.max()

def get_r2_scores_Lasso(df,target,features,scores,cv,groups):
    r2 = []
    features = scores.index
    unique_scores = np.unique(scores.values)
    for score in unique_scores: 
        cv = LeaveOneGroupOut().split(df, groups=groups)
        alpha = get_lassoCV_alpha(df,target,features,cv)
        lasso = Lasso(alpha,max_iter=20000,normalize=False)
        pred = cross_val_predict(lasso, df.loc[:, features], target, cv=LeaveOneGroupOut(), groups=groups)
        r2.append(get_r2(target,pred))
        scores = scores[scores.values > score]
        features = scores.index
    return unique_scores, np.array(r2)

def get_r2_scores_RFR(df,target,features,scores,groups):
    thr_score = 0.
    thr = []
    r2 = []
    features = scores.index
    while (thr_score < 1) and (len(features) > 1):
        cv = LeaveOneGroupOut().split(df, groups=groups)
        if len(features)*5 < 200: 
            rfr = RandomizedSearchCV(RandomForestRegressor(n_estimators=200, criterion='mse'), scoring='r2', param_distributions={'max_features': np.arange(1,len(features)), 'min_samples_leaf': np.arange(1,6)}, cv=cv, n_iter=len(features))
        else: 
            rfr = RandomizedSearchCV(RandomForestRegressor(n_estimators=200, criterion='mse'), scoring='r2', param_distributions={'max_features': np.arange(1,len(features)), 'min_samples_leaf': np.arange(1,6)}, cv=cv, n_iter=100)
        rfr.fit(df.loc[:, features], target, groups=groups)
        max_features = rfr.best_params_['max_features']
        min_samples_leaf = rfr.best_params_['min_samples_leaf']
        rfr_final = RandomForestRegressor(n_estimators=200, max_features=max_features, min_samples_leaf=min_samples_leaf)
        pred = cross_val_predict(rfr_final, df.loc[:, features], target, cv=LeaveOneGroupOut(), groups=groups)
        r2.append(get_r2(target,pred))
        #thr_score = scores.values[scores.shape[0]-1]
        thr_score += 0.01
        thr.append(thr_score)
        scores = scores[scores.values > thr_score]
        features = scores.index
    return np.array(thr), np.array(r2)

def perform_Boruta(n_trees, p, n_leaf, df, target, features): 
    rf = RandomForestRegressor(n_estimators=n_trees, max_features=p, min_samples_leaf=n_leaf)
    feat_selector = BorutaPy(rf, n_estimators=n_trees, verbose=0, random_state=1, alpha=0.05, two_step=False, max_iter=300)
    feat_selector.fit(df.loc[:,features].values, target.values)
    result = pd.DataFrame(feat_selector.ranking_, index=features, columns = ['Boruta ranking'])
    result.loc[features,'Boruta score'] = pd.DataFrame(feat_selector.imp_history_).mean().values
    return result