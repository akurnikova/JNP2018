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

mouse_list = {#'GluReaCh23', #RF
              'GluReaCh24', #RF
              'GluReaCh25', #RF
           #  'GluReaCh26',  #no
              'GluReaCh27',##RF/rostral IRT
              'GluReaCh29', #RF?
           #   'GluReaCh30', ##IRT
              'GluReaCh31', ##IRT
           #   'GluReaCh33', ##IRT?? sparse
              'GluReaCh35', #RF
           #   'GluReaCh37', ##IRT? very caudal
              'GluReaCh38', ##IRT
           #   'GluReaCh41', ##IRT? very caudal
            #  'GluReaCh42', ##IRT
                       }

dfX = get_10ms_dataframes(mouse_name = 'GluReaCh31')

dfX = dfX.reset_index(drop=True)
dfXfilt = dfX.loc[dfX['varLpre']<0.001]
#ind_br_filt = np.ones((len(dfX),))
#for i in range(0,len(ind_br_filt)):
#    ind_br_filt[i] = dfX['clip_br_rate'][i][400] < 3.5
#dfXfilt = dfX#.loc[ind_br_filt]
dfXfilt = dfXfilt.reset_index(drop=True)

A = np.vstack(dfXfilt['clipX'].as_matrix())
Br = np.vstack(dfXfilt['clip_br_rate'].as_matrix())

#fig, ax_new = plt.subplots(2,1, sharex=True)
#ax_new[0].imshow(np.vstack(A[:,450:650]),vmin=-0.5 , vmax=0.2,cmap = 'bone', interpolation='none', extent=[0,dfXfilt['clipX'].count(),-50,150])
#ax_new[0].plot(dfXfilt['max_deflection_time'],dfXfilt.index,'r.', markersize = 2)
#ax_new[1].plot(np.arange(-50,150),np.mean(A[:,450:650],0))

fig, ax = plt.subplots(2,1, sharex=True)
ax[0].plot(np.arange(-200,350),np.mean(A[:,300:850],0))
ax[1].plot(np.arange(-200,350),np.mean(Br[:,300:850],0))

