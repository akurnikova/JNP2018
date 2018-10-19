# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 21:38:19 2016

@author: asya
"""

from __future__ import division
import pickle
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import seaborn as sns

'''
'''

filelist = {'/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO01_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO02_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO03_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO04_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO05_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO06_odor_wk5.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO07_odor_wk5.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO09_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO10_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DIO/dat_preproc/DIO11_odor_wk6.pckl'
            }


type_by_file = {'DIO01':1,'DIO02':2,'DIO03':3,'DIO04':4,'DIO05':5,'DIO06':6,'DIO07':7,'DIO09':9,'DIO10':10,'DIO11':11}
type_by_file = {'DIO01':1,'DIO02':2,'DIO03':-1,'DIO04':-2,'DIO05':3,'DIO06':4,'DIO07':5,'DIO09':-3,'DIO10':6,'DIO11':-4}

plot_individual = 0
sort_breathe = 0
BR_rate_cutoff = 40
var_cutoff = 0.1
trials_max = 25

maxima_all = []
trialtype_all = []
sidetype_all = []

for fle in filelist:
    f = open(fle)
    dat = pickle.load(f)
    f.close()
    
    #%% Make triggers
    
    TS = np.arange(0, len(dat['stim1']))/dat['fpsStim']
    
    stim_trig1 = np.where(np.diff(dat['stim1'])>3)
    stim_trig_T1 = TS[stim_trig1]
    
    stim_trig2 = np.where(np.diff(dat['stim2'])>3)
    stim_trig_T2 = TS[stim_trig2]
    
    TX = np.arange(0, len(dat['L_filt']))/dat['fpsX']
    
    TX_s1 = np.empty(np.shape(stim_trig_T1)) 
    for i,s in enumerate(stim_trig_T1):
        TX_s1[i] = np.where((TX >= s))[0][0]
    
    TX_s2 = np.empty(np.shape(stim_trig_T2)) 
    for i,s in enumerate(stim_trig_T2):
        TX_s2[i] = np.where((TX >= s))[0][0]
        
    #%% Make averaging structures
    
    t_pre = 0.5
    t_post = 1.5


    nX = int(np.floor(dat['fpsX']*(t_pre+t_post)))
    L_av1 = np.zeros((nX,len(TX_s1)))
    L_av2 = np.zeros((nX,len(TX_s2)))
    V_av1 = np.zeros((nX,len(TX_s1)))
    V_av2 = np.zeros((nX,len(TX_s2)))
    
    for i, t in enumerate(TX_s1):
        i_pre = round(t_pre*dat['fpsX'])
        i_post = round(t_post*dat['fpsX'])
        if t+i_post > len(dat['L_filt']): continue
        if t-i_pre < 0: continue
        
        pre_clip = dat['L_filt'][int(t-dat['fpsX']):int(t)]
        pre_clip_V = dat['V_filt'][int(t-dat['fpsX']):int(t)]
        
        if sort_breathe:
            if dat['breathe']['instrate'][int(t)] > BR_rate_cutoff: continue
        else:
            if np.var(pre_clip)>var_cutoff: continue
    
            
        L_av1[:,i] = dat['L_filt'][int(t-i_pre):int(t+i_post)]
        V_av1[:,i] = dat['V_filt'][int(t-i_pre):int(t+i_post)]
        
    for i, t in enumerate(TX_s2):
        i_pre = round(t_pre*dat['fpsX'])
        i_post = round(t_post*dat['fpsX'])
        if t+i_post > len(dat['L_filt']): continue
        if t-i_pre < 0: continue
            
        pre_clip = dat['L_filt'][int(t-dat['fpsX']):int(t)]
        pre_clip_V = dat['V_filt'][int(t-dat['fpsX']):int(t)]
        
        if sort_breathe:
            if dat['breathe']['instrate'][int(t)] > BR_rate_cutoff: continue
        else:
            if np.var(pre_clip)>var_cutoff: continue
    
        
        L_av2[:,i] = dat['L_filt'][int(t-i_pre):int(t+i_post)]
        V_av2[:,i] = dat['V_filt'][int(t-i_pre):int(t+i_post)]


    #%% Create array of peaks for boxplot
    L1_pre = L_av1[0:int(dat['fpsX']*t_pre),:]
    L2_pre = L_av2[0:int(dat['fpsX']*t_pre),:]
    M1 = np.mean(L_av1[750:800,:],axis = 0)-np.mean(L1_pre,axis = 0)
    M2 = np.mean(-L_av2[750:800,:],axis = 0)-np.mean(L2_pre,axis = 0)
    maxima=np.hstack((M1,M2))
    
    ## generate labels
    ch1 = np.chararray(np.shape(M1), itemsize=10)
    ch1[:] = 'control'
    ch2 = np.chararray(np.shape(M2), itemsize=10)
    ch2[:] = 'lesion'
    sidetype = np.hstack((ch1,ch2))
    
    trialtype=np.hstack((M1*0+type_by_file[fle[54:59]],M2*0+type_by_file[fle[54:59]]))
    
    maxima_all = np.hstack((maxima_all,maxima))
    trialtype_all = np.hstack((trialtype_all,trialtype))
    sidetype_all = np.hstack((sidetype_all,sidetype))

    
    #%%
    range_av = range(0,min(trials_max,np.shape(L_av1)[1]-1))
    print range_av
    
    L_av1 = L_av1[:,range_av]
    V_av1 = V_av1[:,range_av]
    L_av2 = L_av2[:,range_av]
    V_av2 = V_av2[:,range_av]
    
    
    TX_clip = range(-int(dat['fpsX']*t_pre),int(dat['fpsX']*t_post))/dat['fpsX']
    TX_clip = TX_clip-0.0097
    
    #%%
    if plot_individual == 1:
        plt.figure()
        ax1 = plt.subplot(211)
        plt.plot(TX_clip,np.mean(L_av1,axis = 1),label='Lat Rstim',color = 'b')
        plt.fill_between(TX_clip,np.mean(L_av1,axis = 1)+st.sem(L_av1,axis = 1),np.mean(L_av1,axis = 1)-st.sem(L_av1,axis = 1),color = 'b', alpha = 0.3)
        plt.plot(TX_clip,-np.mean(L_av2,axis = 1),label='-Lat Lstim',color = 'r')
        plt.fill_between(TX_clip,-np.mean(L_av2,axis = 1)+st.sem(L_av2,axis = 1),-np.mean(L_av2,axis = 1)-st.sem(L_av2,axis = 1),color = 'r', alpha = 0.3)
        plt.xlim(-0.5,1.5)
        plt.ylim(-0.5,1.5)
        plt.legend()
        
        ax1 = plt.subplot(212)
        plt.plot(TX_clip,np.mean(V_av1,axis = 1),label='Vert Rstim',color = 'b')
        plt.fill_between(TX_clip,np.mean(V_av1,axis = 1)+st.sem(V_av1,axis = 1),np.mean(V_av1,axis = 1)-st.sem(V_av1,axis = 1),color = 'b', alpha = 0.3)
        plt.plot(TX_clip,np.mean(V_av2,axis = 1),label='Vert Lstim',color = 'r')
        plt.fill_between(TX_clip,np.mean(V_av2,axis = 1)+st.sem(V_av2,axis = 1),np.mean(V_av2,axis = 1)-st.sem(V_av2,axis = 1),color = 'r', alpha = 0.3)
        plt.legend()
        plt.xlim(-0.5,1.5)
        plt.ylim(-0.5,1.5)
        plt.show()    
    
#%%
plt.figure()
sns.set_style("whitegrid")
sns.boxplot(x = maxima_all,y = trialtype_all,hue = sidetype_all,orient='h',showmeans=True)
plt.figure()
sns.violinplot(x = sidetype_all,y = maxima_all,hue = trialtype_all,showmeans=True)#,split=False)#,bw=.5)
#ax = sns.distplot(maxima_all[trialtype_all == 0], rug=True, hist=False,color = 'b')

#%%

all_means = []
all_labels = []
plt.figure()
for i in type_by_file:
    print type_by_file[i]
    maxima = maxima_all[trialtype_all==type_by_file[i]]
    sides = sidetype_all[trialtype_all==type_by_file[i]]
    mean_lesion = np.mean(maxima[sides=='lesion'])
    mean_control = np.mean(maxima[sides=='control'])
    all_means.append(mean_control)
    all_labels.append('control')
    if type_by_file[i]>0:
        all_means.append(mean_lesion)
        all_labels.append('lesion')
    else:
        all_means.append(mean_lesion)
        all_labels.append('nt_lesions')
plt.figure()
sns.boxplot(x = all_means,y = all_labels,orient='h')
sns.stripplot(x = all_means,y = all_labels,orient='h',size=10, marker="o", edgecolor="gray", alpha=.4)
