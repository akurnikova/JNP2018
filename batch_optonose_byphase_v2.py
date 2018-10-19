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

mouse_list = {#'GluReaCh23', #RF
              #'GluReaCh23b', #RF
              #'GluReaCh24', #RF
              'GluReaCh33', 
              #'GluReaCh37',
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

#%%
#plt.figure()

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
        dfXfilt = dfX.loc[dfX['varLpre']<0.02]
        dfXfilt = dfXfilt.reset_index(drop=True)
        dfX = dfXfilt
        dfX.replace([np.inf, -np.inf], np.nan)
        dfX = dfX.dropna()
        if len(dfX)<20: continue
        # Fit curve
        dfX_fitT = dfX.loc[np.isfinite(dfX['max_deflection_time'])]
        popt_deflection, pcov_deflection = curve_fit(func_cos, dfX['projected_phase_at_stim'], dfX['max_deflection_V'])
        popt_T, pcov_T = curve_fit(func_cos, dfX['projected_phase_at_stim'], dfX['max_deflection_time'])
        
        D = -popt_deflection[1]+np.pi*(popt_deflection[0]>0)
        DT = -popt_T[1]+np.pi*(popt_T[0]>0)
        
        
        diff = D-DT
        if diff < -np.pi: diff= diff+2*np.pi
#        D = D%(2*pi)
        
        
        errT = np.sqrt(pcov_T[1,1])
        
        plt.plot(RC,np.mean(dfX['max_deflection_V']/dfX['max_deflection_L']),'go')
        '''
        D = D%(2*pi)
#        if D < np.pi: D = D+ 2*np.pi
        if abs(popt_T[0])<0.5: continue
        if DT < 2.5: DT = DT+ 2*np.pi
        DT = DT%(2*pi)
        subplot(211)
        plt.plot(RC,abs(popt_T[0]), 'ko')
        subplot(212)
        plt.plot(RC,DT, 'ko')
        plt.plot([RC,RC],[DT-errT,DT+errT], 'k-')
        
        
        D = D%(2*pi)
        if D < np.pi: D = D+ 2*np.pi
        subplot(211)
        plt.plot(RC,abs(popt_deflection[0]), 'ko')
        subplot(212)
        plt.plot(RC,D, 'ko')
        '''
#numpy.polynomial.polynomial.polyfit(x, y, deg, rcond=None, full=False, w=weights)
