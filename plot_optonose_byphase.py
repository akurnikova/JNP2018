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

dfX = get_10ms_dataframes(mouse_name = 'GluReaCh38')
dfXfilt = dfX.loc[dfX['varLpre']<0.02]
dfXfilt = dfXfilt.reset_index(drop=True)
dfX = dfXfilt
dfX = dfX.dropna()

#Sort out trials_with a breath rate change
i_const =np.abs(dfX['br_rate_at_stim']-dfX['br_rate_at_100ms'])<=0.5
dfX = dfX[i_const]

#%% Fit cosine to stim phase
def func_cos(X,a,b,c):
    return a*np.cos(X+b)+c

##########################
#%% Make plots by breathing phase at stimulus
##########

# Fit curve
dfX_fitT = dfX.loc[np.isfinite(dfX['max_deflection_time'])]
popt_time, pcov_time = curve_fit(func_cos, dfX_fitT['projected_phase_at_stim'], dfX_fitT['max_deflection_time'])
popt_deflection, pcov_deflection = curve_fit(func_cos, dfX['projected_phase_at_stim'], dfX['max_deflection_L'])

# Display plot
fig, ax_new = plt.subplots(2,1, sharey=False)
dfX.plot(x = 'projected_phase_at_stim',y = 'max_deflection_L',style = 'o',xlim = [0,6.28], ax = ax_new[0])
dfX.plot(x = 'projected_phase_at_stim',y = 'max_deflection_time',style = 'o',xlim = [0,6.28], ax = ax_new[1])

ax_new[0].plot(np.arange(0,2*np.pi,np.pi/15), func_cos(np.arange(0,2*np.pi,np.pi/15), *popt_deflection), 'g--')
ax_new[1].plot(np.arange(0,2*np.pi,np.pi/15), func_cos(np.arange(0,2*np.pi,np.pi/15), *popt_time), 'g--')

plot_phase_boxplots = 0
if plot_phase_boxplots:
    rounding_N = 3
    dfX['phase_pi'] = np.round((dfX['projected_phase_at_stim']/np.pi-1)*rounding_N,0)/rounding_N
    sns.boxplot(x='phase_pi', y='max_deflection_L', data=dfX, ax = ax_new[0], showfliers=False,color = 'w')
    sns.boxplot(x='phase_pi', y='max_deflection_time', data=dfX, ax = ax_new[1], showfliers=False,color = 'w')
ax_new[0].set_ylim([-0.5,0.5])
ax_new[1].set_ylim([25,75])


if plot_vertical:
    fig, ax_new = plt.subplots(2,1, sharey=False)
    dfX.plot(x = 'projected_phase_at_stim',y = 'max_deflection_V',style = 'o',xlim = [0,6.28], ax = ax_new[0])
    dfX.plot(x = 'projected_phase_at_stim',y = 'max_deflection_timeV',style = 'o',xlim = [0,6.28], ax = ax_new[1])
    sns.boxplot(x='phase_pi', y='max_deflection_V', data=dfX, ax = ax_new[0], showfliers=False,color = 'w')
    sns.boxplot(x='phase_pi', y='max_deflection_timeV', data=dfX, ax = ax_new[1], showfliers=False,color = 'w')
    ax_new[0].set_ylim([0,0.5])
    ax_new[1].set_ylim([25,75])
    

##########################
#%% Make plots by yes or no deflection
##########################

fig, ax_new = plt.subplots(2,1, sharey=False)
dfX['is_deflection'] = np.isfinite(dfX['max_deflection_time'])
dfX.hist('projected_phase_at_stim', by = 'is_deflection', ax = ax_new[0])



##########################
#%% Make plots by breathing phase at peak
##########################
  
## Peak by phase analysis

# Fit curve
dfX_fitT = dfX.loc[np.isfinite(dfX['max_deflection_time'])]
popt_time, pcov_time = curve_fit(func_cos, dfX_fitT['br_phase_at_peak'], dfX_fitT['max_deflection_time'])
popt_deflection, pcov_deflection = curve_fit(func_cos, dfX['br_phase_at_peak'], dfX['max_deflection_L'])

# Display plot
fig, ax_new = plt.subplots(2,1, sharey = False)
dfX_fitT.plot(x = 'br_phase_at_peak',y = 'max_deflection_L',style = 'o',xlim = [0,6.28], ax = ax_new[0])
dfX_fitT.plot(x = 'br_phase_at_peak',y = 'max_deflection_time',style = 'o',xlim = [0,6.28], ax = ax_new[1])

ax_new[0].plot(np.arange(0,2*np.pi,np.pi/15), func_cos(np.arange(0,2*np.pi,np.pi/15), *popt_deflection), 'g--')
ax_new[1].plot(np.arange(0,2*np.pi,np.pi/15), func_cos(np.arange(0,2*np.pi,np.pi/15), *popt_time), 'g--')

plot_phase_boxplots = 1
if plot_phase_boxplots:
    rounding_N = 3
    dfX['phase_pi'] = np.round((dfX['br_phase_at_peak']/np.pi-1)*rounding_N,0)/rounding_N
    sns.boxplot(x='phase_pi', y='max_deflection_L', data=dfX, ax = ax_new[0], showfliers=False,color = 'w')
    sns.boxplot(x='phase_pi', y='max_deflection_time', data=dfX, ax = ax_new[1], showfliers=False,color = 'w')
#ax_new[0].set_ylim([-0.5,0.5])
ax_new[1].set_ylim([25,75])

if 0: ## plot vertical
    fig, ax_new = plt.subplots(2,1, sharey = False)
    dfX.plot(x = 'br_phase_at_stim',y = 'max_deflection_V',style = 'o',xlim = [0,6.28], ax = ax_new[0])
    dfX.plot(x = 'br_phase_at_stim',y = 'max_deflection_timeV',style = 'o',xlim = [0,6.28], ax = ax_new[1])
    sns.boxplot(x='phase_pi', y='max_deflection_V', data=dfX, ax = ax_new[0], showfliers=False,color = 'w')
    sns.boxplot(x='phase_pi', y='max_deflection_timeV', data=dfX, ax = ax_new[1], showfliers=False,color = 'w')
    ax_new[0].set_ylim([0,0.5])
    ax_new[1].set_ylim([25,75])
    
##########################
#%% Make hist breathing phase
##########################

fig, ax_new = plt.subplots(2,1, sharex=True)
#dfX['is_deflection'] = np.isfinite(dfX['max_deflection_time'])
dfX_fitT.hist(['br_phase_at_peak','br_phase_at_stim'], ax = ax_new)

