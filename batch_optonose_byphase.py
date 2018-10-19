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

mouse_list = {'GluReaCh23', #RF
              'GluReaCh24', #RF
              'GluReaCh33', #RF
              'GluReaCh37', #RF
              'GluReaCh25', #RF
              'GluReaCh27',##RF/rostral IRT
              'GluReaCh29', #RF?
              'GluReaCh31', ##IRT
              'GluReaCh35', #RF
              'GluReaCh38', ##IRT
              'GluReaCh42', ##IRT
                       }
cRF = (1,0,1)
cNIRT = (0,1,0)
stack_to_color = {'GluReaCh23':cRF, ###+-+-
               'GluReaCh24':cRF, ###+-+-
               'GluReaCh25':cRF, ###++++
               'GluReaCh27':cRF, ###+-+-
               'GluReaCh29':cRF,
               'GluReaCh30':cNIRT, 
               'GluReaCh31':cNIRT, ###++++
               'GluReaCh33':cNIRT, 
               'GluReaCh35':cRF, ###++++
               'GluReaCh37':cNIRT, 
               'GluReaCh38':cNIRT, ###++++
               'GluReaCh41':cNIRT, 
               'GluReaCh42':cNIRT, ###+-+-
                       }


cLesion = (1.,0.8,0.)
cLesionIRT = (0.95,  0.35,  0.5)
#cLess = (0.6,0.3,0.17)#(0.56,0.3,0.05)
cLeast = (0.1,0.1,0.1)
cControl = 	(0.2, 0.4, 0.8)

stack_to_color = {'GluReaCh23':cLeast, ###+-+-
               'GluReaCh24':cLeast, ###+-+-
               'GluReaCh25':cLesion, ###++++
               'GluReaCh26':cLeast, ###????
               'GluReaCh27':cLesion, ###+-+-
               'GluReaCh28':cControl, ###NO LESION
               'GluReaCh29':cLess,
               'GluReaCh30':cLeast, 
               'GluReaCh31':cLesionIRT, ###++++
               'GluReaCh32':cControl, ###NO LESION
               'GluReaCh33':cLeast, 
               'GluReaCh34':cControl, ###NO LESION
               'GluReaCh35':cLesion, ###++++
               'GluReaCh36':cControl, ###NO LESION
               'GluReaCh37':cLeast, 
               'GluReaCh38':cLesionIRT, ###++++
               'GluReaCh39':cControl, ###NO LESION
               'GluReaCh41':cLeast, 
               'GluReaCh42':cLess, ###+-+-
                       }


'''
#list_coords = list_coords.set_index('ms')
list_coords2 = list_coords.loc[['GR25_L','GR27_L','GR29_L','GR31_L','GR35_L','GR38_L','GR42_L'],:]
df_sort = list_coords2.sort_values(by = 'cx', axis=0, ascending=True)
pal = sns.diverging_palette(10, 220, s=85,l=60, n=7)
#pal_rev = np.vstack((pal[0:3][::-1],pal[9:18][::-1]))
df_sort['color'] = pal#pal_rev.tolist()
'''

#fig, ax_new = plt.subplots(2,1, sharey=False)
fig = plt.figure()
ax0 = fig.add_subplot(121, projection='polar')
ax1 = fig.add_subplot(122, projection='polar')

def vector_sum(phase,r):
    x = r*np.cos(phase)
    y = r*np.sin(phase)
    cx = np.mean(x)
    cy = np.mean(y)
    r_mean = (cx**2+cy**2)**(0.5)
    p_mean = arctan(cx/cy)+np.pi
    return r_mean, p_mean
    
    
for ms_name in mouse_list:
    dfX = get_10ms_dataframes(mouse_name = ms_name)
    dfXfilt = dfX.loc[dfX['varLpre']<0.02]
    dfXfilt = dfXfilt.reset_index(drop=True)
    dfX = dfXfilt
    dfX = dfX.dropna()

    ##########################
    #%% Make plots by breathing phase at stimulus
    ##########
    
    # Fit curve
    dfX_fitT = dfX.loc[np.isfinite(dfX['max_deflection_time'])]
    popt_time, pcov_time = curve_fit(func_cos, dfX_fitT['projected_phase_at_stim'], dfX_fitT['max_deflection_time'])
    popt_deflection, pcov_deflection = curve_fit(func_cos, dfX['projected_phase_at_stim'], dfX['max_deflection_L'])
    
 
    phases = np.linspace(0,2*np.pi,16)
    r = func_cos(phases, *popt_deflection)
    cosine_fit_time =  func_cos_nomean(phases, *popt_time)
    
    #ax0.scatter([-popt_deflection],max(cosine_fit_deflection), color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'][0:3],s=500*popt_deflection[2])
    #ax1.scatter(phases[cosine_fit_time==max(cosine_fit_time)][0],max(cosine_fit_time),  color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'][0:3],s=500*popt_deflection[2])
    
    #ax0.plot(phases,func_cos(phases, *popt_deflection), color = stack_to_color[ms_name])
    #ax1.plot(phases,func_cos(phases, *popt_time), color = stack_to_color[ms_name])
    #color = df_sort.loc['GR'+ms_name[8:10]+'_L']['color'][0:3]
    
    col = stack_to_color[ms_name]
    
    M_deflection = (func_cos(phases, *popt_deflection) == np.max(func_cos(phases, *popt_deflection)))
    M_deflection = np.where(M_deflection)[0][-1]
    
    #ax0.plot([0,phases[M_deflection]],[0,func_cos(phases, *popt_deflection)[M_deflection]], color = stack_to_color[ms_name])
    
    ax0.plot([0,-popt_deflection[1]+np.pi*(popt_deflection[0]<0)],[0,abs(popt_deflection[0])], color = col)

    #ax0.scatter(dfX['projected_phase_at_stim'],dfX['max_deflection_L'], color = stack_to_color[ms_name])
    #ax1.scatter(dfX['projected_phase_at_stim'],dfX['max_deflection_time'], color = stack_to_color[ms_name])
    
    #r,p = vector_sum(dfX['projected_phase_at_stim'],dfX['max_deflection_L'])
    #ax0.plot([0,p],[0,r], color = stack_to_color[ms_name])
ax0.set_rmin(0)

fig_title = '/home/asya/Documents/data/mouse_nose/by_phase_modulation_depth.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

#%%
plt.figure()
#list_coords = list_coords.set_index('ms')
for ms_name in mouse_list:
    if ms_name == 'GluReaCh24': continue
    dfX = get_20ms_dataframes(mouse_name = ms_name)
    dfXfilt = dfX.loc[dfX['varLpre']<0.02]
    dfXfilt = dfXfilt.reset_index(drop=True)
    dfX = dfXfilt
    dfX = dfX.dropna()
    # Fit curve
    dfX_fitT = dfX.loc[np.isfinite(dfX['max_deflection_time'])]
    popt_time, pcov_time = curve_fit(func_cos, dfX_fitT['projected_phase_at_stim'], dfX_fitT['max_deflection_time'])
    popt_deflection, pcov_deflection = curve_fit(func_cos, dfX['projected_phase_at_stim'], dfX['max_deflection_L'])
    
 
    RC = list_coords.loc['GR'+ms_name[8:10]+'_L']['cx']
    
    
#    M_deflection = (func_cos(phases, *popt_deflection) == np.max(func_cos(phases, *popt_deflection)))
#    M_deflection = np.where(M_deflection)[0][-1]
    D = -popt_deflection[1]+np.pi*(popt_deflection[0]<0)
    #if D < 2: D = D+2*np.pi
    #if D < 2: D = D+2*np.pi
    #D = D-2*pi
    err = sqrt(pcov_deflection[1,1])

    plt.plot(RC,D, 'ko')
    plt.plot([RC,RC],[D-err,D+err], 'k-')

#numpy.polynomial.polynomial.polyfit(x, y, deg, rcond=None, full=False, w=weights)
