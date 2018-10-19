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
import pandas as pd

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

filelist = {'/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA01_odor_wk5.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA02_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA03_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA04_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA05_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA06_odor_wk5.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA07_odor_wk5.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA09_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA10_odor_wk6.pckl',
            '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA11_odor_wk6.pckl'
            }


type_by_file = {'DTA01':1,'DTA02':2,'DTA03':3,'DTA04':4,'DTA05':5,'DTA06':6,'DTA07':7,'DTA09':9,'DTA10':10,'DTA11':11}
type_by_file = {'DTA01':1,'DTA02':2,'DTA03':-3,'DTA04':-4,'DTA05':5,'DTA06':6,'DTA07':7,'DTA09':-9,'DTA10':10,'DTA11':-11}


plot_individual = 0
sort_breathe = 0
var_cutoff = 0.1
trials_max = 15

maxima_all = []
maxima_all_V = []
trialtype_all = []
sidetype_all = []

L_all = np.empty((2000,0), float)
V_all = np.empty((2000,0), float)

L_C = np.empty((2000,0), float)
L_L = np.empty((2000,0), float)

ALL_L = dict()
for fle in filelist:
    f = open(fle)
    dat = pickle.load(f)
    f.close()
    
    dat['L_filt'] = 1.1*dat['L_filt']/1.3
    
    ALL_L[fle[54:59]] = dat['L_filt']
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
        
        if np.var(pre_clip)>var_cutoff: continue

        L_av2[:,i] = dat['L_filt'][int(t-i_pre):int(t+i_post)]
        V_av2[:,i] = dat['V_filt'][int(t-i_pre):int(t+i_post)]

    
    #%%
    range_av = range(0,min(trials_max,np.shape(L_av1)[1]-1))
    print range_av
    
    L_av1 = L_av1[:,range_av]
    V_av1 = V_av1[:,range_av]
    L_av2 = L_av2[:,range_av]
    V_av2 = V_av2[:,range_av]
    
    
    #%% Create array of peaks for boxplot
    L1_pre = L_av1[0:int(dat['fpsX']*t_pre),:]
    L2_pre = L_av2[0:int(dat['fpsX']*t_pre),:]
    
    M1 = np.mean(L_av1[850:950,:],axis = 0)#-np.mean(L1_pre,axis = 0)
    M2 = np.mean(-L_av2[850:950,:],axis = 0)#-np.mean(L2_pre,axis = 0)
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


    L_all = np.hstack((L_all,L_av1))
    L_all = np.hstack((L_all,L_av2))
    V_all = np.hstack((V_all,V_av1))
    V_all = np.hstack((V_all,V_av2))
    
    L_C = np.hstack((L_all,L_av1))
    L_L = np.hstack((L_all,L_av2))
    
    
    TX_clip = range(-int(dat['fpsX']*t_pre),int(dat['fpsX']*t_post))/dat['fpsX']
    TX_clip = TX_clip-0.0097
   

#%%
'''
plt.figure()
ax1 = plt.subplot(211)
MN1 = np.mean(L_all[sides=='control'],axis = 1)
plt.plot(TX_clip,label='Lat Rstim',color = 'b')
plt.fill_between(TX_clip,np.mean(L_av1,axis = 1)+st.sem(L_av1,axis = 1),np.mean(L_av1,axis = 1)-st.sem(L_av1,axis = 1),color = 'b', alpha = 0.3)
plt.xlim(-0.5,1.5)
plt.ylim(-0.5,1.5)
plt.legend()
   '''     
#%%
plt.figure()
sns.set_style("whitegrid")
sns.boxplot(x = maxima_all,y = trialtype_all,hue = sidetype_all,orient='h',showmeans=True)


#%%#

#%%
dat = pd.DataFrame()
dat_pairs = pd.DataFrame()
dat_means = pd.DataFrame()
all_means = []
all_stderr = []
all_labels = []

control_side_trace = []
true_control_trace = []
lesion_trace = []
nt_lesion_trace = []

fig = plt.figure()
ax = fig.add_subplot(111)

for i in type_by_file:
    print type_by_file[i]
    maxima = maxima_all[trialtype_all==type_by_file[i]]
    sides = sidetype_all[trialtype_all==type_by_file[i]]
    mean_lesion = np.mean(maxima[sides=='lesion'])
    mean_control = np.mean(maxima[sides=='control'])
    stderr_lesion = st.sem(maxima[sides=='lesion'])
    stderr_control = st.sem(maxima[sides=='control'])
    
    
    if type_by_file[i]>100: ## control
        plt.plot(0,mean_control,'bo')
        plt.plot([0,0],[mean_control+stderr_control,
                         mean_control-stderr_control],'b-')
    
        plt.plot(0,mean_control,'bo')
        plt.plot([0,0],[mean_lesion+stderr_lesion,
                         mean_lesion-stderr_lesion],'b-')
        
        
    elif np.abs(type_by_file[i])>20 and type_by_file[i]<100:
        if type_by_file[i]>0: ## el_lesion
            continue
        else: ## el_nt_lesion
            plt.plot(2,mean_lesion,'yo')
            plt.plot([2,2],[mean_lesion+stderr_lesion,
                         mean_lesion-stderr_lesion],'y-')    
            
            plt.plot(1,mean_control,'go')
            plt.plot([1,1],[mean_control+stderr_control,
                         mean_control-stderr_control],'g-')
            
            
    elif abs(type_by_file[i])<20:
        all_means.append(mean_control)
        all_stderr.append(stderr_control)
        all_labels.append('control')
        if type_by_file[i]>0:
            plt.plot(3,mean_lesion,'mo')
            plt.plot([3,3],[mean_lesion+stderr_lesion,
                         mean_lesion-stderr_lesion],'m-')    
            
            plt.plot(1,mean_control,'go')
            plt.plot([1,1],[mean_control+stderr_control,
                         mean_control-stderr_control],'g-')
            
        else:
            
            dat_means = dat_means.append(pd.DataFrame([['lesion',mean_control,'C'],
                                           ['lesion',mean_lesion,'NT']],
                            columns = ['Type','mean','side']))
            
            plt.plot(2,mean_lesion,'co')
            plt.plot([2,2],[mean_lesion+stderr_lesion,
                         mean_lesion-stderr_lesion],'c-')    
            
            plt.plot(1,mean_control,'go')
            plt.plot([1,1],[mean_control+stderr_control,
                         mean_control-stderr_control],'g-')
            

fig_title = '/home/asya/Documents/data/odor_lesion_DTA/averages_all_dots.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

#%%
dat = pd.DataFrame()
dat_pairs = pd.DataFrame()
dat_means = pd.DataFrame()
all_means = []
all_stderr = []
all_labels = []

control_side_trace = []
true_control_trace = []
lesion_trace = []
nt_lesion_trace = []

for i in type_by_file:
    print type_by_file[i]
    maxima = maxima_all[trialtype_all==type_by_file[i]]
    sides = sidetype_all[trialtype_all==type_by_file[i]]
    mean_lesion = np.mean(maxima[sides=='lesion'])
    mean_control = np.mean(maxima[sides=='control'])
    stderr_lesion = st.sem(maxima[sides=='lesion'])
    stderr_control = st.sem(maxima[sides=='control'])
    
    
    if type_by_file[i]>100:
        all_means.append(mean_control)
        all_means.append(mean_lesion)
        all_stderr.append(stderr_control)
        all_stderr.append(stderr_lesion)
        all_labels.append('true_control')
        all_labels.append('true_control')
        
        dat_means =dat_means.append(pd.DataFrame([['true_control',mean_control,'C'],
                                  ['true_control',mean_lesion,'C']],
                        columns = ['Type','mean','side']))
        df2 = pd.DataFrame()
        df2['maxima'] = maxima[sides=='control']
        df2['Type'] = 'true_control'
        df2['side'] = 'C'
        dat = dat.append(df2)
        
        df2['maxima'] = maxima[sides=='lesion']
        df2['Type'] = 'true_control'
        df2['side'] = 'C'
        dat = dat.append(df2)
        
    elif np.abs(type_by_file[i])>20 and type_by_file[i]<100:
        all_means.append(mean_control)
        all_stderr.append(stderr_control)
        all_labels.append('control')
        if type_by_file[i]>0:
            all_means.append(mean_lesion)
            all_stderr.append(stderr_lesion)
            all_labels.append('el_lesion')            
        else:
            all_means.append(mean_lesion)
            all_stderr.append(stderr_lesion)
            all_labels.append('el_nt_lesions')
            dat_means = dat_means.append(pd.DataFrame([['lesion',mean_control,'C'],
                                           ['lesion',mean_lesion,'EL']],
                            columns = ['Type','mean','side']))
            
            df2 = pd.DataFrame()
            df2['maxima'] = maxima[sides=='lesion']
            df2['Type'] = 'lesion'
            df2['side'] = 'EL'
            dat = dat.append(df2)
            
            df2['maxima'] = maxima[sides=='control']
            df2['Type'] = 'lesion'
            df2['side'] = 'C'
            dat = dat.append(df2)
        
            
            dat_pairs = dat_pairs.append(pd.DataFrame([['el_lesion',mean_control,mean_lesion]],
                            columns = ['Type','mean_control','mean_lesion']),ignore_index = True)
    elif abs(type_by_file[i])<10:
        all_means.append(mean_control)
        all_stderr.append(stderr_control)
        all_labels.append('control')
        if type_by_file[i]>0:
            all_means.append(mean_lesion)
            all_stderr.append(stderr_lesion)
            all_labels.append('lesion')
            dat_means = dat_means.append(pd.DataFrame([['lesion',mean_control,'C'],
                                           ['lesion',mean_lesion,'L']],
                            columns = ['Type','mean','side']))
            
            dat_pairs = dat_pairs.append(pd.DataFrame([['lesion',mean_control,mean_lesion]],
                            columns = ['Type','mean_control','mean_lesion']),ignore_index = True)
            
            df2 = pd.DataFrame()
            df2['maxima'] = maxima[sides=='lesion']
            df2['Type'] = 'lesion'
            df2['side'] = 'L'
            dat = dat.append(df2)
            
            df2['maxima'] = maxima[sides=='control']
            df2['Type'] = 'lesion'
            df2['side'] = 'C'
            dat = dat.append(df2)
            
        else:
            all_means.append(mean_lesion)
            all_stderr.append(stderr_lesion)
            all_labels.append('nt_lesions')
            dat_means = dat_means.append(pd.DataFrame([['lesion',mean_control,'C'],
                                           ['lesion',mean_lesion,'NT']],
                            columns = ['Type','mean','side']))
            
            dat_pairs = dat_pairs.append(pd.DataFrame([['nt_lesions',mean_control,mean_lesion]],
                            columns = ['Type','mean_control','mean_lesion']),ignore_index = True)
            
            df2 = pd.DataFrame()
            df2['maxima'] = maxima[sides=='lesion']
            df2['Type'] = 'lesion'
            df2['side'] = 'NT'
            dat = dat.append(df2)
            
            df2['maxima'] = maxima[sides=='control']
            df2['Type'] = 'lesion'
            df2['side'] = 'C'
            dat = dat.append(df2)


fig = plt.figure()
ax = fig.add_subplot(111)

sns.boxplot(x = all_means,y = all_labels,orient='h')
sns.stripplot(x = all_means,y = all_labels,orient='h',size=10, marker="o", edgecolor="gray", alpha=.4)

#fig_title = '/home/asya/Documents/data/odor_lesion_DTA/boxplots_v2.pdf'
#plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
#%%
fig = plt.figure()
ax = fig.add_subplot(111)

#sns.boxplot(x = 'mean',y = 'Type',hue='side',data = dat_means,orient='h')
sns.boxplot(x = 'maxima',y = 'Type',hue='side',data = dat,orient='h',dodge=True)

#sns.stripplot(x = all_means,y = all_labels,orient='h',size=10, marker="o", edgecolor="gray", alpha=8)


x1 = [dat_pairs.loc[0]['mean_control'],dat_pairs.loc[0]['mean_lesion']]

for i in dat_pairs.index:
    if dat_pairs.loc[i].Type == 'nt_lesions':
        plt.plot([dat_pairs.loc[i]['mean_control'],dat_pairs.loc[i]['mean_lesion']],[-0.25,0.25],'b')        
        
    if dat_pairs.loc[i].Type == 'lesion':
        plt.plot([dat_pairs.loc[i]['mean_control'],dat_pairs.loc[i]['mean_lesion']],[-0.25,0.25],'r')
        
    if dat_pairs.loc[i].Type == 'el_lesion':
        plt.plot([dat_pairs.loc[i]['mean_control'],dat_pairs.loc[i]['mean_lesion']],[-0.25,0.25],'g')

#fig_title = '/home/asya/Documents/data/odor_lesion_DTA/boxplots_v1.pdf'
#plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


#%%

fig = plt.figure()
ax = fig.add_subplot(111)
## True Control
T1 = -L_all[:,(trialtype_all>100)*(sidetype_all == 'lesion')]
T2 = L_all[:,(trialtype_all>100)*(sidetype_all == 'control')]
trace_true_control = np.hstack((T1,T2))

T1 = -L_all[:,(trialtype_all>100)*(sidetype_all == 'lesion')]
T2 = L_all[:,(trialtype_all>100)*(sidetype_all == 'control')]
trace_true_control = np.hstack((T1,T2))

MN = np.mean(trace_true_control,axis = 1)
STERR = np.std(trace_true_control,axis = 1)/(trace_true_control.shape[1])**0.5
plt.plot(MN,'b')
plt.plot(MN+STERR,'b--')
plt.plot(MN-STERR,'b--')

## EL lesion
trace_el_lesion = -L_all[:,(np.abs(trialtype_all)>20)*(trialtype_all<100)*(sidetype_all == 'lesion')]
MN = np.mean(trace_el_lesion,axis = 1)
STERR = np.std(trace_el_lesion,axis = 1)/(trace_el_lesion.shape[1])**0.5
plt.plot(MN,'y')
plt.plot(MN+STERR,'y--')
plt.plot(MN-STERR,'y--')

## lesion
trace_lesion = -L_all[:,(np.abs(trialtype_all)<20)*(trialtype_all>0)*(sidetype_all == 'lesion')]
MN = np.mean(trace_lesion,axis = 1)
STERR = np.std(trace_lesion,axis = 1)/(trace_lesion.shape[1])**0.5
plt.plot(MN,'m')
plt.plot(MN+STERR,'m--')
plt.plot(MN-STERR,'m--')

##nt lesion
trace_nt_lesion = -L_all[:,(np.abs(trialtype_all)<20)*(trialtype_all<0)*(sidetype_all == 'lesion')]
MN = np.mean(trace_nt_lesion,axis = 1)
STERR = np.std(trace_nt_lesion,axis = 1)/(trace_nt_lesion.shape[1])**0.5
plt.plot(MN,'c')
plt.plot(MN+STERR,'c--')
plt.plot(MN-STERR,'c--')

## All control
TC1 = L_all[:,(np.abs(trialtype_all)>20)*(trialtype_all<100)*(sidetype_all == 'control')]
TC2 = L_all[:,(np.abs(trialtype_all)<20)*(sidetype_all == 'control')]
#TC3 = L_all[:,(np.abs(trialtype_all)<10)*(trialtype_all>0)*(sidetype_all == 'control')]
trace_control = np.hstack((TC1,TC2))
MN = np.mean(trace_control,axis = 1)
STERR = np.std(trace_control,axis = 1)/(trace_control.shape[1])**0.5
plt.plot(MN,'g')
plt.plot(MN+STERR,'g--')
plt.plot(MN-STERR,'g--')
## Yellow = Control
## c = nt.lesion
## g = el lesion
## m = lesion
## b = true control

fig_title = '/home/asya/Documents/data/odor_lesion_DTA/traces_all_average.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


#%%
fig = plt.figure()
#plt.plot(np.arange(0,0.35,0.01),np.arange(0,0.35,0.01),color = (.1,.1,.1))
for stack in ALL_L:
    a = histogram(ALL_L[stack],bins = np.arange(-4,4,0.25))
    a_tot = np.float(sum(a[0]))
    if type_by_file[stack] <0:
        color = [1,0,0]
    if type_by_file[stack] >0:
        color = [0,1,1]
    #plt.plot(sum(a[0][-5:])/a_tot,sum(a[0][0:6])/a_tot,'o',color = color)
    plt.plot(a[1][1:],a[0]/a_tot,color = color)

fig_title = '/home/asya/Documents/data/odor_lesion_DTA/dist_pos_DTA.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


#%%
fig = plt.figure()
#plt.plot(np.arange(0,0.35,0.01),np.arange(0,0.35,0.01),color = (.1,.1,.1))
n = 31
L = np.zeros((n,))
NT = np.zeros((n,))
TimeL = 0 
TimeNT = 0 
for stack in ALL_L.keys():
    a = histogram(ALL_L[stack],bins = np.arange(-4,4,0.25))
    a_tot = np.float(sum(a[0]))
    if type_by_file[stack] > 0:
        color = 'm'
        L = L+a[0]/a_tot
        TimeL = TimeL + len(ALL_L[stack])
    if type_by_file[stack] <0:
        color = 'c'
        NT = NT + a[0]/a_tot
        TimeNT = TimeNT + len(ALL_L[stack])
    #plt.plot(sum(a[0][-5:])/a_tot,sum(a[0][0:6])/a_tot,'o',color = color)
    #plt.plot(a[1][1:],(a[0]/a_tot),color = color)

plt.plot(a[1][1:],log10(L/float(6)),'m')
plt.plot(a[1][1:],log10(NT/float(4)),'c')
fig_title = '/home/asya/Documents/data/odor_lesion/dist_pos_DTA.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
    