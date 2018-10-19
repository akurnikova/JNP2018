#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 14:41:17 2018

@author: asya
"""

import pandas as pd


## Lines to determine optimal for a given mouse
#%%
#df = pd.read_pickle('/home/asya/Documents/data/mouse_nose/pandas_data/'+'GluReaCh35')
#df.index.get_level_values(0).unique()
#%%
#print df.loc['dy1',10].index.unique()
#%%
#print df.loc['dy4', 10,355].count()
#%%
plot_deflection_time  = 0
plot_vertical = 0
plot_phase_boxplots  = 0
plot_vertical = 0

def get_5ms_dataframes(mouse_name = 'GluReaCh35'):  
    df = pd.read_pickle('/home/asya/Documents/data/mouse_nose/pandas_data/'+mouse_name[:10])
    dfX = []
    
    if mouse_name == 'GluReaCh27':
        dfX = df.loc['dy2', 5,675] ##
    if mouse_name == 'GluReaCh29':
        dfX = df.loc['dy2', 5,675] ##
    if mouse_name == 'GluReaCh29b':
        dfX = df.loc['dy3', 5,675] ##
    if mouse_name == 'GluReaCh30':
        dfX = df.loc['dy2', 5,675] ##
    if mouse_name == 'GluReaCh31':
        dfX = df.loc['dy1', 5,675] ##
    if mouse_name == 'GluReaCh31b':
        dfX = df.loc['dy2', 5,675] ##
    if mouse_name == 'GluReaCh33':
        dfX = df.loc['dy1', 5,690] ##
    if mouse_name == 'GluReaCh35':
        dfX = df.loc['dy1', 5,690] ##
    if mouse_name == 'GluReaCh37':
        dfX = df.loc['dy1', 5,690] ##
    if mouse_name == 'GluReaCh38':
        dfX = df.loc['dy1', 5,685] ##
    if mouse_name == 'GluReaCh42':
        dfX = df.loc['dy1', 5,675] ##
    return dfX
    
def get_20ms_dataframes(mouse_name = 'GluReaCh25'):  
    df = pd.read_pickle('/home/asya/Documents/data/mouse_nose/pandas_data/'+mouse_name[:10])
    dfX = []
  #  if mouse_name == 'GluReaCh23':
  #     dfX = df.loc['dy3', 20,355] ## 
      #  dfX = df.loc['dy4', 20,350] ##
    
    if mouse_name == 'GluReaCh24':
        dfX = df.loc['dy1', 20,355] ##
    #    dfX = df.loc['dy2', 20,355] ## 
    if mouse_name == 'GluReaCh25':
        dfX = df.loc['day', 20,350] ##
    if mouse_name == 'GluReaCh25b':
        dfX = df.loc['dy1', 20,355] ##
    if mouse_name == 'GluReaCh26':
        dfX = df.loc['dy2', 20,350] ##
    if mouse_name == 'GluReaCh27':
        dfX = df.loc['dy2', 20,280] ##
    if mouse_name == 'GluReaCh28':
        dfX = df.loc['dy3', 20,285] ##
    if mouse_name == 'GluReaCh29':
       # dfX = df.loc['dy3', 20,280] ## 
        dfX = df.loc['dy2', 20,280] ## 
    if mouse_name == 'GluReaCh30':
        dfX = df.loc['dy2', 20,280] ## 
    if mouse_name == 'GluReaCh31':
        dfX = df.loc['dy1', 20,285] ##
    if mouse_name == 'GluReaCh32':
        dfX = df.loc['dy1', 20,280]
    if mouse_name == 'GluReaCh33':
        dfX =  df.loc['dy1', 20,280]
    if mouse_name == 'GluReaCh34':
        dfX = df.loc['pt4', 20,175]
    if mouse_name == 'GluReaCh35':
        dfX = df.loc['dy1', 20,175] ## 253
    if mouse_name == 'GluReaCh36':
        dfX = df.loc['dy1', 20,175] ## 192
    if mouse_name == 'GluReaCh37':
        dfX = df.loc['dy1', 20,175]
  #  if mouse_name == 'GluReaCh38':
  #      dfX = df.loc['dy1', 20,175]
    if mouse_name == 'GluReaCh39':
        dfX = df.loc['dy1', 20,175]
    if mouse_name == 'GluReaCh42':
        dfX = df.loc['dy1', 20,165] ## 434
        
    return dfX

def get_10ms_dataframes(mouse_name = 'GluReaCh25'):  
    df = pd.read_pickle('/home/asya/Documents/data/mouse_nose/pandas_data/'+mouse_name[:10])
    dfX = []
    if mouse_name == 'GluReaCh23':
        #dfX = df.loc['dy3', 10,355] ## 165 stimuli
        dfX = df.loc['dy4', 10,350] ## 331 stimuli
    if mouse_name == 'GluReaCh24':
        dfX = df.loc['dy1', 10,355] ## 100 stimuli    
    if mouse_name == 'GluReaCh25':
        dfX = df.loc['day', 10,355] ## 201 stimuli
    if mouse_name == 'GluReaCh25b':
        dfX = df.loc['dy4', 10,355] ## 201 stimuli
    if mouse_name == 'GluReaCh23b':
        dfX = df.loc['dy3', 10,355] ## 165 stimuli        
    if mouse_name == 'GluReaCh27b':
        dfX = df.loc['dy3', 10,360] ## 218 stimuli    
    if mouse_name == 'GluReaCh27c':
        dfX = df.loc['dy3', 10,445] ## 164 stimuli    
    if mouse_name == 'GluReaCh27':
        dfX = df.loc['dy2', 10,390] ## 182 stimuli    
    if mouse_name == 'GluReaCh28':
        dfX = df.loc['dy3', 10,390] ## 188 stimuli    
    if mouse_name == 'GluReaCh29':
        dfX = df.loc['dy2', 10,390] ## 214 stimuli
    if mouse_name == 'GluReaCh29b':
        dfX = df.loc['dy3', 10,390] ## 161 stimuli   
    if mouse_name == 'GluReaCh30':
        dfX = df.loc['dy2', 10,390.0] ## 286 stimuli
        #dfX = dfX[0:200]    
    if mouse_name == 'GluReaCh31':
        dfX = df.loc['dy1', 10,390] ## 320 stimuli
    if mouse_name == 'GluReaCh32':
        dfX = df.loc['dy1', 10,390] ## 228 stimuli ????
    if mouse_name == 'GluReaCh33':
        dfX = df.loc['dy1', 10,390] ## 295 stimuli    
    if mouse_name == 'GluReaCh34':
        dfX = df.loc['pt1', 10,390] ## 160
    if mouse_name == 'GluReaCh35':
        dfX = df.loc['dy1', 10,390] ## 253  
    if mouse_name == 'GluReaCh35b':
        dfX = df.loc['dy1', 10,345] ## 253  
    if mouse_name == 'GluReaCh36':
        dfX = df.loc['dy1', 10,345] ## 192    
    if mouse_name == 'GluReaCh37':
        dfX = df.loc['dy1', 10,345] ## 377    
    if mouse_name == 'GluReaCh38':
        dfX = df.loc['dy2', 10,345] ## 276    
    if mouse_name == 'GluReaCh39':
        dfX = df.loc['dy1', 10,345] ## 320
#    if mouse_name == 'GluReaCh40':
#        dfX = df.loc['dy1', 10,390] ## 
    if mouse_name == 'GluReaCh41':
        dfX = df.loc['dy1', 10,345] ## 299
    if mouse_name == 'GluReaCh42':
        dfX = df.loc['dy1', 10,390] ## 434    
    
    return dfX
    
def get_stim_current_dataframes(mouse_name = 'GluReaCh35'):  
    df = pd.read_pickle('/home/asya/Documents/data/mouse_nose/pandas_data/'+mouse_name)
    df2 = df.reorder_levels(['file_day','stim_power','stim_duration','stim_current'], axis=0)
    
    if mouse_name == 'GluReaCh27':
        dfStimCurrent = df.loc['dy3', 10]
        dfStimPower = df2.loc['dy2',7.1]
        
    if mouse_name == 'GluReaCh29':
        dfStimCurrent = df.loc['dy2', 10]
        dfStimPower = df2.loc['dy2',7.1]

    
    if mouse_name == 'GluReaCh30':
        dfStimCurrent = df.loc['dy2', 10]
        dfStimPower = df2.loc['dy2',7.1]

        
    if mouse_name == 'GluReaCh31':
        dfStimCurrent = df.loc['dy1', 10]
        dfStimPower = df2.loc['dy2']
    
    if mouse_name == 'GluReaCh35':
        dfStimCurrent = df.loc['dy1', 10]
        dfStimPower = df2.loc['dy1',7.2].append(df2.loc['dy1',7.4])
    
    return dfStimCurrent, dfStimPower


## Select by breathing rate    
#dfStimCurrent = dfStimCurrent.loc[dfStimCurrent['br_rate_at_stim']<3.5]
#dfStimPower = dfStimPower.loc[dfStimPower['br_rate_at_stim']<3.5]

#dfStimCurrent.reset_index(inplace=True)
#dfStimPower.reset_index(inplace=True)