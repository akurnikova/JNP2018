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

#%% Fit cosine to stim phase
def func_cos(X,a,b,c):
    return a*np.cos(X+b)+c
def func_cos_nomean(X,a,b,c):
    return a*np.cos(X+b)

mouse_list = {'GluReaCh23', 
              'GluReaCh24',
              'GluReaCh25', #RF
           #  'GluReaCh26',  #no
              'GluReaCh27',##RF/rostral IRT
              'GluReaCh29', #RF?
              'GluReaCh30', ##IRT
              'GluReaCh31', ##IRT
              'GluReaCh33', ##IRT?? sparse
              'GluReaCh35', #RF
              'GluReaCh37', ##IRT? very caudal
              'GluReaCh38', ##IRT
              #'GluReaCh41', ##IRT? very caudal
              'GluReaCh42', ##IRT
                       }

mouse_ctrl_list = {'GluReaCh28'}#, 'GluReaCh34', 'GluReaCh36','GluReaCh39'} #control
df_sort = list_coords.sort_values(by = 'cx', axis=0, ascending=True)
pal = sns.diverging_palette(10, 220, s=95,l=80, n=16, center="dark")
pal_rev = np.vstack((pal[0:7][::-1],pal[9:18][::-1]))
df_sort['color'] = pal_rev.tolist()
df_sort = df_sort.set_index('ms')

fig = plt.figure()
ax1 = fig.add_subplot(111)

for ms_name in mouse_ctrl_list.union({'GluReaCh25'}):
    dfX = get_10ms_dataframes(mouse_name = ms_name)
    dfXfilt = dfX.loc[dfX['varLpre']<0.02]
    dfXfilt = dfXfilt.reset_index(drop=True)
    dfX = dfXfilt
    dfX = dfX.dropna()

    dfX_fitT = dfX.loc[np.isfinite(dfX['max_deflection_time'])]
    
    if ms_name in mouse_ctrl_list: col = (0,0,0)
    else: col = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'][0:3]
    
    phase_plotphase_plot = dfX_fitT['time_from_last_br']
    Delta_phase = dfX_fitT['br_phase_at_50ms']-dfX_fitT['br_phase_at_stim']
    Delta_phase[Delta_phase<0] = Delta_phase+2*np.pi
    
    phase_plot = dfX_fitT[['br_phase_at_stim','br_phase_at_50ms']]
    phase_plot = phase_plot.sort_values(by = 'br_phase_at_stim')
    phase_plot['br_phase_at_50_unwrap'] = np.unwrap(phase_plot['br_phase_at_50ms'])#,discont = 2*np.pi)
    
    ax1.scatter(dfX_fitT['time_from_last_br'],dfX_fitT['time_to_next_br'],color = col,s = 3)#,color = col)

#    x_phase = dfX_fitT[dfX_fitT['br_phase_at_stim']<2.5]['br_phase_at_stim']
#    y_phase = phase_plot[dfX_fitT['br_phase_at_stim']<2.5]
#    f = np.polyfit(x_phase, y_phase, 1)
#    fitC = np.poly1d(f)
#    ax1.plot(np.linspace(0,2.5,30),fitC(np.linspace(0,2.5,30)),color = col)
    
#    x_phase = dfX_fitT[dfX_fitT['br_phase_at_stim']>3.5]['br_phase_at_stim']
#    y_phase = phase_plot[dfX_fitT['br_phase_at_stim']>3.5]
#    f = np.polyfit(x_phase, y_phase, 1)
#    fitC = np.poly1d(f)
#    ax1.plot(np.linspace(3.5,6,30),fitC(np.linspace(3.5,6,30)),color = col)
   

ax1.axis('equal')
