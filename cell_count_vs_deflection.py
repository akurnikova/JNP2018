#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 14:48:15 2018

@author: asya
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import os


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
import import_data_for_plots

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

mouse_list = {'GluReaCh23',
              'GluReaCh24',
              'GluReaCh25',
              'GluReaCh27',
              'GluReaCh28',
              'GluReaCh29',
              'GluReaCh30',
              'GluReaCh31',
              'GluReaCh32',
              'GluReaCh33',
              'GluReaCh34',
              'GluReaCh35',
              'GluReaCh36',
              'GluReaCh37',
              'GluReaCh38',
              'GluReaCh39',
              'GluReaCh41',
              'GluReaCh42',
                       }



L_deflection = dict()
V_deflection = dict()
for ms_name in mouse_list:
    stack = 'GR' + ms_name[8:10] + '_L'
        
    dfX = get_10ms_dataframes(mouse_name = ms_name)

    if len(dfX) == 0: continue
    dfX = dfX.reset_index(drop=True)
    dfXfilt = dfX.loc[dfX['varLpre']<0.1]
    dfXfilt = dfXfilt.reset_index(drop=True)
    dfXfilt = dfXfilt.dropna()
    i_const =np.abs(dfXfilt['br_rate_at_stim']-dfXfilt['br_rate_at_100ms'])<=0.5
    dfXfilt = dfXfilt[i_const]


    L_deflection[stack] = np.mean(dfXfilt['max_deflection_L'])
    V_deflection[stack] = np.mean(dfXfilt['max_deflection_V'])
    
#%%
    
def fit_line_with_ci(ax,axn,x,y):
    x = sm.add_constant(x) # constant intercept term
    model = sm.OLS(y, x)
    fitted = model.fit()
    x_pred = np.linspace(x.min(), x.max(), 50)
    x_pred2 = sm.add_constant(x_pred)
    y_pred = fitted.predict(x_pred2)
    ax[axn].plot(x_pred, y_pred, '-', color='darkorchid', linewidth=2)
    y_hat = fitted.predict(x) # x is an array from line 12 above
    y_err = y - y_hat
    mean_x = x.T[1].mean()
    n = len(x)
    dof = n - fitted.df_model - 1

    t = stats.t.ppf(1-0.025, df=dof)
    
    s_err = np.sum(np.power(y_err, 2))
    conf = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((x_pred-mean_x),2) / 
     ((np.sum(np.power(x_pred,2))) - n*(np.power(mean_x,2))))))
    upper = y_pred + abs(conf)
    lower = y_pred - abs(conf)
    ax[axn].fill_between(x_pred, lower, upper, color='darkorchid', alpha=0.4)
    
def fit_line_with_ci_2factor(x,y):
    fig, ax = plt.subplots(2,1, sharex = True,sharey = True)
    
    x = sm.add_constant(x) # constant intercept term
    model = sm.OLS(y, x)
    fitted = model.fit()
    x_pred1 = np.linspace(x.min(0)[1], x.max(0)[1], 50)
    x_pred2 = np.linspace(x.min(0)[2], x.max(0)[2], 50)
    x_pred = sm.add_constant(zip(x_pred1,x_pred2))
    y_pred = fitted.predict(x_pred)
    
    ax[0].plot(x[:,1],y,'o')
    ax[0].plot(x_pred1, y_pred, '-', color='darkorchid', linewidth=2)
    ax[2].plot(x[:,2],y,'o')
    ax[1].plot(x_pred2, y_pred, '-', color='darkorchid', linewidth=2)
    
    #ax[axn].plot(x_pred1, y_pred, '-', color='darkorchid', linewidth=2)
    #ax[axn].plot(x_pred2, y_pred, '-', color='darkorchid', linewidth=2)
    
    
    y_hat = fitted.predict(x) # x is an array from line 12 above
    y_err = y - y_hat
    #mean_x = x.T[1].mean()
    #n = len(x)
    #dof = n - fitted.df_model - 1

    #t = stats.t.ppf(1-0.025, df=dof)
    
    #s_err = np.sum(np.power(y_err, 2))
    #conf = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((x_pred-mean_x),2) / 
    # ((np.sum(np.power(x_pred,2))) - n*(np.power(mean_x,2))))))
    #upper = y_pred + abs(conf)
    #lower = y_pred - abs(conf)
    #ax[axn].fill_between(x_pred, lower, upper, color='darkorchid', alpha=0.4)

fig, ax = plt.subplots(3,1, sharex = True,sharey = True)
 
pd1 = pd.DataFrame.from_dict(L_deflection,orient='index')
pd1.columns = ['L_deflection']
pd1['V_deflection'] = pd.DataFrame.from_dict(V_deflection,orient='index')
pd1['point_count0'] = pd.DataFrame.from_dict(point_count0,orient='index')
pd1['point_count1'] = pd.DataFrame.from_dict(point_count1,orient='index')
pd1['pt_count_sum'] = pd1['point_count0']+pd1['point_count1']
pd1 = pd1.fillna(0)
ax[0].plot(pd1['pt_count_sum'],pd1['L_deflection'],'bo')

x = [pd1.point_count0,pd1.point_count1]
y = pd1.L_deflection
p, V = np.polyfit(x, y, 1, cov=True)

x = pd1[['point_count0','point_count1']].values #pd1.pt_count_sum.values

fit_line_with_ci(ax,0,x,y)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       


#ax[0].plot(x,p[1]+p[0]*x,'b')
#ax[0].plot(x,p[1]+(p[0]+np.sqrt(V[0][0]))*x,'b--')
#ax[0].plot(x,p[1]+(p[0]-np.sqrt(V[0][0]))*x,'b--')

ax[1].plot(pd1['point_count0'],pd1['L_deflection'],'ro')
x = pd1.point_count0.values
y = pd1.L_deflection
fit_line_with_ci(ax,1,x,y)

ax[2].plot(pd1['point_count1'],pd1['L_deflection'],'go')
x = pd1.point_count1.values
y = pd1.L_deflection

fit_line_with_ci(ax,2,x,y)
    