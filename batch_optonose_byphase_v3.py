#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 09:46:44 2018

@author: asya
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import os
import pandas as pd
import seaborn as sns

import import_data_for_plots ## Sets the file to work on and gets the structures

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

#%% Fit cosine to stim phase
def func_cos(X,a,b,c):
    return a*np.cos(X+b)+c
def func_cos_nomean(X,a,b,c):
    return a*np.cos(X+b)+c

mouse_list = {'GluReaCh24', #RF
              'GluReaCh25', 
              'GluReaCh25b', 
              'GluReaCh27',
              'GluReaCh27b',##RF/rostral IRT
              'GluReaCh27c',
              'GluReaCh29',
              'GluReaCh29b',
              'GluReaCh31',
              'GluReaCh31b',
              'GluReaCh35', 
              'GluReaCh35b',
              'GluReaCh38',
              'GluReaCh42',
                       }


list_coords = pd.read_pickle('/home/asya/Documents/data/mouse_nose/list_coords')
list_coords2 = list_coords.loc[['GR24_L','GR25_L','GR27_L','GR29_L','GR31_L','GR35_L','GR38_L','GR42_L'],:]
df_sort = list_coords2.sort_values(by = 'cx', axis=0, ascending=True)
pal = sns.diverging_palette(10, 220, s=85,l=60, n=8)
#pal_rev = np.vstack((pal[0:3][::-1],pal[9:18][::-1]))
df_sort['color'] = pal#pal_rev.tolist()

#%%
fig = plt.figure()
ax0 = fig.add_subplot(221, projection='polar')
ax1 = fig.add_subplot(222, projection='polar')
ax2 = fig.add_subplot(223, projection='polar')


for ms_name in mouse_list:
    for i in [5,10,20]:
        RC = list_coords.loc['GR'+ms_name[8:10]+'_L']['cx']
        
        if i == 5:
            dfX = get_5ms_dataframes(mouse_name = ms_name)
        if i == 10:
            dfX = get_10ms_dataframes(mouse_name = ms_name)
        if i == 20:
            dfX = get_20ms_dataframes(mouse_name = ms_name)
            
        if len(dfX) == 0: continue
        dfXfilt = dfX.loc[dfX['varLpre']<0.001]
        dfXfilt = dfXfilt.reset_index(drop=True)
        dfX = dfXfilt
        dfX.replace([np.inf, -np.inf], np.nan)
        dfX = dfX.dropna()
        if len(dfX)<20: continue
    
        ## Sort trials where movement was too small
        dfX = dfX[dfX['max_deflection_V']>0.05]
        
        #Sort out trials_with a breath rate change
        i_const =np.abs(dfX['br_rate_at_stim']-dfX['br_rate_at_100ms'])<=0.1
        dfBrConst = dfX[i_const]
        
        popt_T, pcov_T = curve_fit(func_cos, dfBrConst['projected_phase_at_stim'], dfBrConst['max_deflection_time'])  
        if abs(popt_T[0]) - np.sqrt(pcov_T[0,0]) > 0:
            DT = -popt_T[1]+np.pi*(popt_T[0]>0)
            errT = np.sqrt(pcov_T[1,1])
            ax0.plot([DT,DT],[0,abs(popt_T[0])], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'])
            ax0.fill_between([DT+errT,DT-errT],[abs(popt_T[0]),abs(popt_T[0])],[0,0], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'],alpha = 0.3)
        
        
        popt_V, pcov_V = curve_fit(func_cos, dfBrConst['projected_phase_at_stim'], dfBrConst['max_deflection_V'])    
        if abs(popt_V[0]) - np.sqrt(pcov_V[0,0]) > 0:
            DV = -popt_V[1]+np.pi*(popt_V[0]<0)
            errV = np.sqrt(pcov_V[1,1])
            ax1.plot([DV,DV],[0,abs(popt_V[0])], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'])
            ax1.fill_between([DV+errV,DV-errV],[abs(popt_V[0]),abs(popt_V[0])],[0,0], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'],alpha = 0.3)


        popt_L, pcov_L = curve_fit(func_cos, dfBrConst['projected_phase_at_stim'], dfBrConst['max_deflection_L'])    
        if abs(popt_L[0]) - np.sqrt(pcov_L[0,0]) > 0:
            DL = -popt_L[1]+np.pi*(popt_L[0]<0)
            errL = np.sqrt(pcov_L[1,1])
            ax2.plot([DL,DL],[0,abs(popt_L[0])], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'])
            ax2.fill_between([DL+errL,DL-errL],[abs(popt_L[0]),abs(popt_L[0])],[0,0], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'],alpha = 0.3)

            
ax0.set_rmin(0)
ax1.set_rmin(0)
ax2.set_rmin(0)

fig_title = '/home/asya/Documents/data/mouse_nose/phases_all_stim.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
#%%
fig = plt.figure()
ax0 = fig.add_subplot(111, projection='polar')
#ax0 = fig.add_subplot(212, projection='polar')
for ms_name in mouse_list:
    for i in [10]:
        RC = list_coords.loc['GR'+ms_name[8:10]+'_L']['cx']
        
        if i == 5:
            dfX = get_5ms_dataframes(mouse_name = ms_name)
        if i == 10:
            dfX = get_10ms_dataframes(mouse_name = ms_name)
        if i == 20:
            dfX = get_20ms_dataframes(mouse_name = ms_name)
            
        if len(dfX) == 0: continue
        dfXfilt = dfX.loc[dfX['varLpre']<0.001]
        dfXfilt = dfXfilt.reset_index(drop=True)
        dfX = dfXfilt
        dfX.replace([np.inf, -np.inf], np.nan)
        dfX = dfX.dropna()
        if len(dfX)<20: continue
    
        #Sort out trials_with a breath rate change
        i_const =np.abs(dfX['br_rate_at_stim']-dfX['br_rate_at_100ms'])<=0.5
        dfBrConst = dfX[i_const]
        
        h = np.histogram(dfBrConst['br_phase_at_peak'],linspace(0,2*np.pi,10))
        ax0.plot(h[1][:-1] + np.diff(h[1]),h[0].astype(float)/np.sum(h[0]), color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'])
        
      #  h = np.histogram(dfBrConst['projected_phase_at_stim'],linspace(0,2*np.pi,16))
      #  ax1.plot(h[1][:-1] + np.diff(h[1]),h[0], color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'])