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
import pandas as pd
import seaborn as sns

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

cLesion = (1.,0.8,0.)
cLesionNIRT = (0.2,0.8,0.8)
cLess = (0.6,0.3,0.17)
cLeast = (0.1,0.1,0.1)
cControl = 	(0.2, 0.4, 0.8)

stack_to_color = {'GluReaCh23':cLeast, ###+-+-
                  'GluReaCh24':cLess, ###+-+-
                  'GluReaCh25':cLesion, ###++++
                  'GluReaCh26':cLeast, ###????
                  'GluReaCh27':cLesion, ###+-+-
                  'GluReaCh28':cControl, ###NO LESION
                  'GluReaCh29':cLess,
                  'GluReaCh30':cLeast, 
                  'GluReaCh31':cLesion, ###++++
                  'GluReaCh32':cControl, ###NO LESION
                  'GluReaCh33':cLeast, 
                  'GluReaCh34':cControl, ###NO LESION
                  'GluReaCh35':cLesion, ###++++
                  'GluReaCh36':cControl, ###NO LESION
                  'GluReaCh37':cLeast, 
                  'GluReaCh38':cLesion, ###++++
                  'GluReaCh39':cControl, ###NO LESION
                  'GluReaCh41':cLeast, 
                  'GluReaCh42':cLesion, ###+-+-
                       }

list_coords = pd.read_pickle('/home/asya/Documents/data/mouse_nose/list_coords')
list_coords2 = list_coords.loc[['GR24_L','GR25_L','GR27_L','GR29_L','GR31_L','GR35_L','GR38_L','GR42_L'],:]
df_sort = list_coords2.sort_values(by = 'cx', axis=0, ascending=True)
pal = sns.diverging_palette(10, 220, s=85,l=60, n=8)
#pal_rev = np.vstack((pal[0:3][::-1],pal[9:18][::-1]))
df_sort['color'] = pal#pal_rev.tolist()

#fig, ax_new = plt.subplots(3,1, sharex=True)
#fig, ax_new = plt.subplots(2,1, sharex=False)


for ms_name in mouse_list:
    if ms_name in ['GluReaCh41','GluReaCh32','GluReaCh34']: continue ## br signal bad

    dfX = get_10ms_dataframes(mouse_name = ms_name)

    if len(dfX) == 0: continue
    dfX = dfX.reset_index(drop=True)
    dfXfilt = dfX.loc[dfX['varLpre']<0.1]
    dfXfilt = dfXfilt.reset_index(drop=True)
    dfXfilt = dfXfilt.dropna()
    stack = 'GR' + ms_name[8:10] + '_L'
    if stack in df_sort.index:
        col = df_sort.loc[stack]['color'][0:3]
    else: 
        col = stack_to_color[ms_name]

    plt.scatter(dfXfilt['br_rate_at_stim'],dfXfilt['br_rate_at_100ms'],color = col,s = 1)
    MN1 = np.mean(dfXfilt['br_rate_at_stim'])
    MN2 = np.mean(dfXfilt['br_rate_at_100ms'])
    ERR1 = np.std(dfXfilt['br_rate_at_stim'])/np.sqrt(dfXfilt['br_rate_at_stim'].shape[0])
    ERR2 = np.mean(dfXfilt['br_rate_at_100ms'])/np.sqrt(dfXfilt['br_rate_at_100ms'].shape[0])
    plt.scatter(MN1,MN2,color = col,s = 30, marker = '>')
    plt.plot([MN1-ERR1,MN1+ERR1],[MN2,MN2],color = col)
    plt.plot([MN1,MN1],[MN2-ERR2,MN2+ERR2],color = col)
    
plt.axis('equal')
    
fig_title = '/home/asya/Documents/data/mouse_nose/breathing_rate_change_100ms.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
