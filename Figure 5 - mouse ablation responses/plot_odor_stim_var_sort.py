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

    
'''
'''

fle = '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA06_odor_wk5.pckl'
plotsave = 0
sort_breathe = 0
BR_rate_cutoff = 40
var_cutoff = 0.1
trials_max = 15

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
    
    pre_clip = dat['L_filt'][int(t-0.5*dat['fpsX']):int(t)]
    pre_clip_V = dat['V_filt'][int(t-0.5*dat['fpsX']):int(t)]
    
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
        
    pre_clip = dat['L_filt'][int(t-0.5*dat['fpsX']):int(t)]
    pre_clip_V = dat['V_filt'][int(t-0.5*dat['fpsX']):int(t)]
    
    if sort_breathe:
        if dat['breathe']['instrate'][int(t)] > BR_rate_cutoff: continue
    else:
        if np.var(pre_clip)>var_cutoff: continue

    
    L_av2[:,i] = dat['L_filt'][int(t-i_pre):int(t+i_post)]
    V_av2[:,i] = dat['V_filt'][int(t-i_pre):int(t+i_post)]


#%%
if plotsave == 0:
    plt.figure()
    plt.subplot(121)
    M1 = np.amax(L_av1,axis = 0)
    M2 = np.amax(-L_av2,axis = 0)
    
    plt.plot(M1,'bo-')
    plt.plot(M2,'ro-')
    
    plt.ylim([0,4])
    
    plt.subplot(122)
    
    histM1 = np.histogram(M1,bins = np.arange(0,4,0.4))
    histM2 = np.histogram(M2,bins = np.arange(0,4,0.4))
    
    plt.barh(histM1[1][1:],histM1[0],0.1,color='b')
    plt.barh(histM1[1][1:]-0.07,histM2[0],0.1,color='r')
    plt.ylim([0,4])

#%%
range_av = range(0,min(trials_max,np.shape(L_av1)[1]))

L_av1 = L_av1[:,range_av]
V_av1 = V_av1[:,range_av]
L_av2 = L_av2[:,range_av]
V_av2 = V_av2[:,range_av]


TX_clip = range(-int(dat['fpsX']*t_pre),int(dat['fpsX']*t_post))/dat['fpsX']
TX_clip = TX_clip-0.0097

#%%
if plotsave == 0:
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
if plotsave == 1:
    fig = plt.figure()
    ax1 = plt.subplot(211)
    plt.plot(TX_clip,np.mean(L_av1,axis = 1),label='Lat Rstim',color = 'b')
    plt.plot(TX_clip,np.mean(L_av1,axis = 1)+st.sem(L_av1,axis = 1),color = 'c')
    plt.plot(TX_clip,np.mean(L_av1,axis = 1)-st.sem(L_av1,axis = 1),color = 'c')
    
    plt.plot(TX_clip,-np.mean(L_av2,axis = 1),label='Lat -Lstim',color = 'r')
    plt.plot(TX_clip,-np.mean(L_av2,axis = 1)+st.sem(L_av2,axis = 1),color = 'm')
    plt.plot(TX_clip,-np.mean(L_av2,axis = 1)-st.sem(L_av2,axis = 1),color = 'm')
    plt.xlim(-0.5,1.5)
    plt.ylim(-0.5,1.5)
    
    plt.legend()
    
    plt.subplot(212, sharex=ax1)
    plt.plot(TX_clip,np.mean(V_av1,axis = 1),label='Vert Rstim',color = 'b')
    plt.plot(TX_clip,np.mean(V_av1,axis = 1)+st.sem(V_av1,axis = 1),color = 'c')
    plt.plot(TX_clip,np.mean(V_av1,axis = 1)-st.sem(V_av1,axis = 1),color = 'c')
    
    plt.plot(TX_clip,np.mean(V_av2,axis = 1),label='Vert Lstim',color = 'r')
    plt.plot(TX_clip,np.mean(V_av2,axis = 1)+st.sem(V_av2,axis = 1),color = 'm')
    plt.plot(TX_clip,np.mean(V_av2,axis = 1)-st.sem(V_av2,axis = 1),color = 'm')
    plt.xlim(-0.5,1.5)
    plt.ylim(-0.5,1.5)
    fig_title = '/home/asya/Documents/data/odor_lesion_DTA/traces_example.pdf'
    plt.savefig(fig_title, format='pdf',dpi=fig.dpi)



plt.legend()
plt.show()


fig = plt.figure()
ax1 = fig.add_subplot(211)
plt.plot(TS,dat['L_filt'])
ax2 = fig.add_subplot(212,sharex = ax1)
plt.plot(TS,dat['stim1'])
plt.plot(TS,dat['stim2'])

fig_title = '/home/asya/Documents/data/odor_lesion_DTA/trace_example1.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)