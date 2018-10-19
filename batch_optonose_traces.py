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
              'GluReaCh26',
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
cControl = 	(0.7, 0.7, 0.7)

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
'''
mouse_list = {'GluReaCh24', 
              'GluReaCh25', #RF
              'GluReaCh27', #RF/Rostral IRT
              'GluReaCh29', #RF
              'GluReaCh31', #IRT
              'GluReaCh35', #RF
              'GluReaCh38', #IRT
              'GluReaCh42', #IRT
                       }

stack_to_color = {'GluReaCh23':(0.5,0.,1.), ###+-+-
               'GluReaCh24':(0.,0.5,1.), ###+-+-
               'GluReaCh25':'r',#(1.,0.8,0.), ###+++
               'GluReaCh27':'k',#(1.,0.8,0.), ###+-+-
               'GluReaCh29':'y',#(1.,0.8,0.),
               'GluReaCh31':'b',#(1.,0.3,0.), ###++++
               'GluReaCh35':'m',#(1.,0.8,0.), ###++++
               'GluReaCh38':'c',#(1.,0.3,0.), ###++++
               'GluReaCh42':'g',#(1.,0.3,0.), ###+-+-
                       }
'''

fig, ax_new = plt.subplots(3,1, sharex=True)
#fig, ax_new = plt.subplots(2,1, sharex=False)


for ms_name in mouse_list:
    dfX = get_10ms_dataframes(mouse_name = ms_name)
    
    if len(dfX) == 0: continue
    dfX = dfX.reset_index(drop=True)
    dfXfilt = dfX.loc[dfX['varLpre']<0.1]
    dfXfilt = dfXfilt.reset_index(drop=True)
    dfXfilt = dfXfilt.dropna()
    
    i_const =np.abs(dfXfilt['br_rate_at_stim']-dfXfilt['br_rate_at_100ms'])<=0.5
    dfXfilt = dfXfilt[i_const]
        
    A = np.vstack(dfXfilt['clipX'].as_matrix())
    V = np.vstack(dfXfilt['clipV'].as_matrix())
    Br = np.vstack(dfXfilt['clip_br_rate'].as_matrix())
    
    if ms_name == 'GluReaCh26': continue

    MN = np.mean(A[:,250:850],0)-np.mean(A[:,450:500])
    STDERR = (np.std(A[:,250:850],0))/np.sqrt(A.shape[0])
    ax_new[0].plot(np.arange(-250,350),MN, color =  stack_to_color[ms_name])
   # ax_new[0].fill_between(np.arange(-250,350),MN+STDERR,MN-STDERR, color =  stack_to_color[ms_name],alpha = 0.3)
    
    MN_V = np.mean(V[:,250:850],0)-np.mean(V[:,450:500])
    STDERR_V = (np.std(V[:,250:850],0))/np.sqrt(V.shape[0])
    ax_new[1].plot(np.arange(-250,350),MN_V, color =  stack_to_color[ms_name])
   # ax_new[1].fill_between(np.arange(-250,350),MN_V+STDERR_V,MN_V-STDERR_V, color =  stack_to_color[ms_name],alpha = 0.3)
   
    if ms_name in ['GluReaCh41','GluReaCh32','GluReaCh34']: continue
    MN_Br = np.mean(Br[:,250:850],0)
    STDERR_Br = (np.std(Br[:,250:850],0))/np.sqrt(Br.shape[0])
    ax_new[2].plot(np.arange(-250,350),MN_Br, color =  stack_to_color[ms_name])
    ax_new[2].fill_between(np.arange(-250,350),MN_Br+STDERR_Br,MN_Br-STDERR_Br, color =  stack_to_color[ms_name],alpha = 0.3)
    
   # ax_new[2].scatter(np.mean(Br[:,450:490],0)-np.mean(Br[:,530:570],0),np.mean(A[:,530:570],0)-np.mean(A[:,450:500]),color =  stack_to_color[ms_name])

fig_title = '/home/asya/Documents/data/mouse_nose/traces_10ms_stim.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
