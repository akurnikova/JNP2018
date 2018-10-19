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

pd1 = pd.DataFrame.from_dict(L_deflection,orient='index')
pd1.columns = ['L_deflection']
pd1['V_deflection'] = pd.DataFrame.from_dict(V_deflection,orient='index')
pd1['point_count0'] = pd.DataFrame.from_dict(point_count0,orient='index')
pd1['point_count1'] = pd.DataFrame.from_dict(point_count1,orient='index')
pd1 = pd1.fillna(0)

Vol0 = (4./3.)*np.pi*np.prod(2*s0)
Vol1 = (4./3.)*np.pi*np.prod(3*s1/2)
pd1['point_D0'] = 3*pd1.point_count0/Vol0
pd1['point_D1'] = 3*pd1.point_count1/Vol1

y = pd1.L_deflection
x = 3*pd1[['point_count0','point_count1']].values #pd1.pt_count_sum.values



fig, ax = plt.subplots(2,1, sharex = True,sharey = True)
x = sm.add_constant(x) # constant intercept term
model = sm.OLS(y, x)
fitted = model.fit()
x_pred1 = np.linspace(x.min(0)[1], x.max(0)[1], 50)
x_pred2 = np.linspace(x.min(0)[2], x.max(0)[2], 50)
x_pred = sm.add_constant(zip(x_pred1,x_pred2))
y_pred = fitted.predict(x_pred)
    
print fitted.summary()


ax[0].plot(x[:,1],y,'o')
ax[0].plot(x_pred1, y_pred, '-', color='darkorchid', linewidth=2)
ax[1].plot(x[:,2],y,'o')
ax[1].plot(x_pred2, y_pred, '-', color='darkorchid', linewidth=2)

#ax[axn].plot(x_pred1, y_pred, '-', color='darkorchid', linewidth=2)
#ax[axn].plot(x_pred2, y_pred, '-', color='darkorchid', linewidth=2)


y_hat = fitted.predict(x) # x is an array from line 12 above
y_err = y - y_hat
mean_x = x.T[1].mean()
n = len(x)
dof = n - fitted.df_model - 1
t = stats.t.ppf(1-0.025, df=dof)
s_err = np.sum(np.power(y_err, 2))

conf1 = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((x_pred1-mean_x),2) / 
 ((np.sum(np.power(x_pred1,2))) - n*(np.power(mean_x,2))))))
upper = y_pred + abs(conf1)
lower = y_pred - abs(conf1)
ax[0].fill_between(x_pred1, lower, upper, color='darkorchid', alpha=0.4)

conf2 = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((x_pred2-mean_x),2) / 
 ((np.sum(np.power(x_pred2,2))) - n*(np.power(mean_x,2))))))
upper = y_pred + abs(conf2)
lower = y_pred - abs(conf2)
ax[1].fill_between(x_pred2, lower, upper, color='darkorchid', alpha=0.4)

## F-test results:
x = pd1[['point_count0','point_count1']].values #pd1.pt_count_sum.values
x = sm.add_constant(x) # constant intercept term
model = sm.OLS(y, x)
results = model.fit()
A = np.identity(len(results.params))
A = A[1:,:]
results.f_test(A)


## F-test results:
x = pd1[['point_count0']].values #pd1.pt_count_sum.values
x = sm.add_constant(x) # constant intercept term
model = sm.OLS(y, x)
results = model.fit()
A = np.identity(len(results.params))
A = A[1:,:]
results.f_test(A)

## F-test results:
x = pd1[['point_count1']].values #pd1.pt_count_sum.values
x = sm.add_constant(x) # constant intercept term
model = sm.OLS(y, x)
results = model.fit()
A = np.identity(len(results.params))
A = A[1:,:]
results.f_test(A)

#%%
import pandas as pd
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
model = ols('L_deflection ~ point_count0+point_count1', data=pd1).fit()
anovaResults = anova_lm(model)
print(anovaResults)
