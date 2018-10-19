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


dfStimCurrent, dfStimPower = get_stim_current_dataframes(mouse_name = 'GluReaCh35')

#dfXfilt = dfX.loc[dfX['varLpre']<0.02]
#dfXfilt = dfXfilt.reset_index(drop=True)
#dfX = dfXfilt

plot_deflection_time = 0
#%%##  Box plots by stim current at constant duration (10 ms)
fig, ax_new = plt.subplots(1,2, sharey=True)
dfStimCurrent = dfStimCurrent.reset_index()
sns.boxplot(x='stim_current', y='max_deflection_L', data=dfStimCurrent,color = 'b', showfliers=False,ax = ax_new[0])
sns.stripplot(x='stim_current', y='max_deflection_L', data=dfStimCurrent,ax = ax_new[0])

sns.boxplot(x='stim_current', y='max_deflection_V', data=dfStimCurrent,color = 'w', showfliers=False,ax = ax_new[1])
sns.stripplot(x='stim_current', y='max_deflection_V', data=dfStimCurrent,ax = ax_new[1])

fig.suptitle('Max nose deflection by stim current (10 ms)')
ax_new[0].set_ylabel('Max deflection (mm)')
ax_new[1].set_ylabel('Max deflection (mm)')
ax_new[0].set_title('Lateral')
ax_new[1].set_title('Vertical')

fig_title = '/home/asya/Documents/data/mouse_nose/GR35_Stim_current_nose.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

if plot_deflection_time:
    fig, ax_new = plt.subplots(1,2, sharey=True)
    sns.boxplot(x='stim_current', y='max_deflection_time', data=dfStimCurrent,color = 'b', showfliers=False,ax = ax_new[0])
    sns.stripplot(x='stim_current', y='max_deflection_time', data=dfStimCurrent,ax = ax_new[0])

    sns.boxplot(x='stim_current', y='max_deflection_timeV', data=dfStimCurrent,color = 'w', showfliers=False,ax = ax_new[1])
    sns.stripplot(x='stim_current', y='max_deflection_timeV', data=dfStimCurrent,ax = ax_new[1])

    fig.suptitle('Time of max nose deflection by stim current (10 ms)')
    ax_new[0].set_ylabel('Deflection time (ms)')
    ax_new[1].set_ylabel('Deflection time (ms)')
    ax_new[0].set_title('Lateral')
    ax_new[1].set_title('Vertical')
    ax_new[0].set_ylim([30,80])
    fig_title = '/home/asya/Documents/data/mouse_nose/GR35_Stim_current_time.pdf'
    plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


#%%##  Box plots by stim duration at constant stim power
dfStimPower = dfStimPower.reset_index()
fig, ax_new = plt.subplots(1,2, sharey=False)
sns.boxplot(x='stim_duration', y='max_deflection_L', data=dfStimPower,color = 'b', showfliers=False,ax = ax_new[0])
sns.stripplot(x='stim_duration', y='max_deflection_L', data=dfStimPower,ax = ax_new[0])
sns.boxplot(x='stim_duration', y='max_deflection_V', data=dfStimPower,color = 'w', showfliers=False,ax = ax_new[1])
sns.stripplot(x='stim_duration', y='max_deflection_L', data=dfStimPower,ax = ax_new[1])
fig.suptitle('Max nose deflection by stim duration (constant power)')
ax_new[0].set_ylabel('Max deflection (mm)')
ax_new[1].set_ylabel('Max deflection (mm)')
ax_new[0].set_title('Lateral')
ax_new[1].set_title('Vertical')

fig_title = '/home/asya/Documents/data/mouse_nose/GR35_Stim_power_nose.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


if plot_deflection_time:
    fig, ax_new = plt.subplots(1,2, sharey=False)
    sns.boxplot(x='stim_duration', y='max_deflection_time', data=dfStimPower,color = 'b', showfliers=False,ax = ax_new[0])
    sns.stripplot(x='stim_duration', y='max_deflection_time', data=dfStimPower,ax = ax_new[0])

    sns.boxplot(x='stim_duration', y='max_deflection_timeV', data=dfStimPower,color = 'w', showfliers=False,ax = ax_new[1])
    sns.stripplot(x='stim_duration', y='max_deflection_timeV', data=dfStimPower,ax = ax_new[1])
    fig.suptitle('Time of max nose deflection by stim duration (constant power)')
    ax_new[0].set_ylabel('Deflection time (ms)')
    ax_new[1].set_ylabel('Deflection time (ms)')
    ax_new[0].set_title('Lateral')
    ax_new[1].set_title('Vertical')
    
    fig_title = '/home/asya/Documents/data/mouse_nose/GR35_Stim_power_time.pdf'
    plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

#%%##  plot trace at constant stim power
#%%##  Box plots by stim duration at constant stim power

dfX = dfStimPower

dfX = dfX.reset_index(drop=False)
dfXfilt = dfX.loc[dfX['varLpre']<0.001]
dfXfilt = dfXfilt.reset_index(drop=True)

df5 = dfXfilt.loc[dfXfilt['stim_duration']==5]
A5 = np.vstack(df5['clipX'].as_matrix())
df10 = dfXfilt.loc[dfXfilt['stim_duration']==10]
A10 = np.vstack(df10['clipX'].as_matrix())
df20 = dfXfilt.loc[dfXfilt['stim_duration']==20]
A20 = np.vstack(df20['clipX'].as_matrix())

fig, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(np.arange(-200,350),-np.mean(A5[:,300:850],0), color = sns.color_palette("Reds")[5])
ax[0].plot(np.arange(-200,350),-np.mean(A10[:,300:850],0),color = sns.color_palette("Reds")[3])
ax[0].plot(np.arange(-200,350),-np.mean(A20[:,300:850],0),color = sns.color_palette("Reds")[1])
legend([5,10,20])

dfX = dfStimCurrent

dfX = dfX.reset_index(drop=False)
dfXfilt = dfX.loc[dfX['varLpre']<0.001]
dfXfilt = dfXfilt.reset_index(drop=True)

df1 = dfXfilt.loc[dfXfilt['stim_current']==225]
A1 = np.vstack(df1['clipX'].as_matrix())

df2 = dfXfilt.loc[dfXfilt['stim_current']==300]
A2 = np.vstack(df2['clipX'].as_matrix())

df3 = dfXfilt.loc[dfXfilt['stim_current']==345]
A3 = np.vstack(df3['clipX'].as_matrix())

df4 = dfXfilt.loc[dfXfilt['stim_current']==390]
A4 = np.vstack(df4['clipX'].as_matrix())


df5 = dfXfilt.loc[dfXfilt['stim_current']==445]
A5 = np.vstack(df4['clipX'].as_matrix())


ax[1].plot(np.arange(-200,350),-np.mean(A1[:,300:850],0),color = sns.color_palette("Blues")[5])
ax[1].plot(np.arange(-200,350),-np.mean(A2[:,300:850],0),color = sns.color_palette("Blues")[4])
ax[1].plot(np.arange(-200,350),-np.mean(A3[:,300:850],0),color = sns.color_palette("Blues")[3])
ax[1].plot(np.arange(-200,350),-np.mean(A4[:,300:850],0),color = sns.color_palette("Blues")[2])
legend([225,300,345,390])
#ax[1].plot(np.arange(-200,350),np.mean(A5[:,300:850],0),color = sns.color_palette("Blues")[1])
fig_title = '/home/asya/Documents/data/mouse_nose/GR35_Traces_by_const_power_and_current.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

