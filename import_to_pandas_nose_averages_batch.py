# -*- coding: utf-8 -*-
"""
Created on Thu Mar 2 13:19:09 2017

@author: asya
"""

from __future__ import division
import numpy as np
import pickle

import os
import pandas as pd
import matplotlib.pyplot as plt
import peakutils
import sys
sys.path.append('/home/asya/Documents/data/utilities/')
from utils_filtering import util_lowpass

def round_to_base(x, base=5):
    return int(base * round(float(x)/base))


mouse_list = {'GluReaCh23','GluReaCh24','GluReaCh25','GluReaCh26',
        'GluReaCh27','GluReaCh28','GluReaCh29','GluReaCh30',
        'GluReaCh32','GluReaCh31','GluReaCh33','GluReaCh34',
                   'GluReaCh35','GluReaCh36','GluReaCh37','GluReaCh38',
                   'GluReaCh39','GluReaCh40','GluReaCh41','GluReaCh42'}
#mouse_name = 'GluReaCh33'
for mouse_name in mouse_list:#
     #{'GluReaCh25'}:
    input_dir = '/home/asya/Documents/data/mouse_nose/preproc_data/'+mouse_name+'/'
    input_dir_list = os.listdir(input_dir)
    
    file_list = []
    for names in input_dir_list:
        if names.startswith("GluR"):
            file_list.append(names)
    
    current_setting = {'GluReaCh23':350,'GluReaCh24':350,'GluReaCh25':350,'GluReaCh26':350,
                       'GluReaCh27':750,'GluReaCh28':750,'GluReaCh29':750,'GluReaCh30':750,
                       'GluReaCh31':750,'GluReaCh32':750,'GluReaCh33':750,'GluReaCh34':750,
                       'GluReaCh35':750,'GluReaCh36':750,'GluReaCh37':750,'GluReaCh38':750,
                       'GluReaCh39':750,'GluReaCh40':750,'GluReaCh41':750,'GluReaCh42':750,
                       }
    fps = 1000.
    ib = int(0.5*fps) ## pre stim timeduration
    ie = int(2*fps)  ## Post stim time
    norm_time = 100 # ms
    Tclip = np.arange(-ib,ie)/fps
        
    
    df = pd.DataFrame()
    data_loaded = {'max_deflection_L':[],'max_deflection_time':[],'varLpre':[],
                   'max_deflection_V':[],'max_deflection_timeV':[],
                   'stim_duration':[], 'stim_current':[],
                   'br_rate_at_stim':[],'br_phase_at_stim':[],'br_phase_at_peak':[],
                   'br_phase_at_50ms':[],'br_phase_at_100ms':[],
                   'br_rate_at_peak':[], 'br_rate_at_50ms':[], 'br_rate_at_100ms':[],
                   'br_rate_at_N50ms':[], 'br_rate_at_N100ms':[],
                   'time_from_last_br':[],'time_to_next_br':[],
                   'projected_phase_at_stim':[],
                   'file_day':[],'clipX':[],'clipV':[],'clip_br_rate':[]}
    
    

    
    for fle in file_list:
        f = open(input_dir+fle)
        dat = pickle.load(f)
        f.close()
        
        if dat.has_key('breathe'):
            adjusted_phases = dat['breathe']['phases']
            BR_Rate = util_lowpass(dat['breathe']['instrate'],1000,10,1)
        T_all = np.arange(0,dat['noseX'].shape[0])/fps
        ## start the alignment
        for i,s in enumerate(dat['stim_trig_idx']):
            
            df2 = dict()
            
            clipX = dat['noseX'][(s-ib):(s+ie)]#-np.mean(dat['noseX'][(s-norm_time):s])
            clipV = dat['noseV'][(s-ib):(s+ie)]#-np.mean(dat['noseV'][(s-norm_time):s])
            clipX = 0.4*clipX # to mm - assume maximum camera zoom for all mice
            clipV = 0.4*clipV # to mm
            
            if len(clipX)<ie+ib: continue
            
            if dat.has_key('breathe'):
                clip_br_rate = BR_Rate[(s-ib):(s+ie)]
                phase_stim = dat['breathe']['phases'][s+1]
                data_loaded['br_rate_at_stim'].append(BR_Rate[s])
                data_loaded['br_phase_at_stim'].append(phase_stim)
                data_loaded['clip_br_rate'].append(clip_br_rate)
                            
                data_loaded['br_rate_at_100ms'].append(BR_Rate[s+100])
                data_loaded['br_rate_at_50ms'].append(BR_Rate[s+50])
                data_loaded['br_rate_at_N100ms'].append(BR_Rate[s-100])
                data_loaded['br_rate_at_N50ms'].append(BR_Rate[s-50])
                data_loaded['br_phase_at_100ms'].append(dat['breathe']['phases'][s+100])
                data_loaded['br_phase_at_50ms'].append(dat['breathe']['phases'][s+50])
                
                ##Find previous breath:
                time_at_stim = T_all[s]
                rise = dat['breathe']['rise']/fps
                prev_rise = np.max(rise[rise < time_at_stim])
                next_rise = np.min(rise[rise > time_at_stim])
                data_loaded['time_from_last_br'].append(time_at_stim-prev_rise)
                data_loaded['time_to_next_br'].append(next_rise-time_at_stim)
                
                ## Projected br phase
                prev_rises = (fps*(rise[rise < time_at_stim][-2:]))
                index_to_prev = s-prev_rises[1]
                projected_phase_at_stim = dat['breathe']['phases'][int(prev_rises[0]+index_to_prev)]
                data_loaded['projected_phase_at_stim'].append(projected_phase_at_stim)
#                adjusted_phases[index_to_prev:s] = dat['breathe']['phases'][int(prev_rises[0]):int(prev_rises[0]+index_to_prev)]
                
                #next_rise = int(fps*np.min(rise[rise > time_at_stim]))
                #adjusted_phases[s:next_rise] = dat['breathe']['phases'][prev_rises[0]:(prev_rises[0]+index_to_prev)]
                
            else:
                
                data_loaded['br_rate_at_stim'].append(np.nan)
                data_loaded['br_phase_at_stim'].append(np.nan)
                data_loaded['clip_br_rate'].append(np.nan)
                data_loaded['br_rate_at_100ms'].append(np.nan)
                data_loaded['br_rate_at_50ms'].append(np.nan)
                data_loaded['br_rate_at_N100ms'].append(np.nan)
                data_loaded['br_rate_at_N50ms'].append(np.nan)
                data_loaded['br_phase_at_100ms'].append(np.nan)
                data_loaded['br_phase_at_50ms'].append(np.nan)
                data_loaded['time_from_last_br'].append(np.nan)
                data_loaded['time_to_next_br'].append(np.nan)
                data_loaded['projected_phase_at_stim'].append(np.nan)
    
            data_loaded['stim_duration'].append(dat['stim_time'][i])
            data_loaded['stim_current'].append(round_to_base(current_setting[mouse_name]*dat['stim_Volt'][i]/5.))
            if mouse_name == 'GluReaCh27' or mouse_name == 'GluReaCh28' or \
                mouse_name == 'GluReaCh29' or mouse_name == 'GluReaCh24' or\
                mouse_name == 'GluReaCh25' or mouse_name == 'GluReaCh26':
                 data_loaded['file_day'].append(fle[13:16])
            else:
                 data_loaded['file_day'].append(fle[11:14])
            
            data_loaded['clipX'].append(clipX)
            data_loaded['clipV'].append(clipV)
            
            data_loaded['varLpre'].append(np.var(dat['noseX'][(s-norm_time):s]))
    
    
            ## Perform the peak_detection
    #        threshold = 3*np.std(clipX[(ib-100):ib]) #Set the threshold as standard deviation movement 100ms before stim
            prior_mean = np.mean(clipX[(ib-100):ib]) #Mean position 100ms before stim
    #        
    #        index_left = peakutils.indexes(-clipX, thres=threshold, min_dist = 10)
    #        index_right = peakutils.indexes(clipX, thres=threshold, min_dist = 10)
    #        indices =  np.union1d(index_left,index_right)
    #        indices = indices[indices<(ib+70)]
            
            ## Find the local minimum or maximum within the preset range of times. If none, set time as nan
            clip_trace = clipX[(ib+30):(ib+70)]
            max_deviation = np.where(clip_trace == np.min(clip_trace))[0][0]+ib+30
            max_deflection = clipX[max_deviation]-prior_mean
            max_deflection_time = max_deviation - ib
            
            if max_deviation== (ib+30) or max_deviation == (ib+69):
                max_deviation = np.where(clip_trace == np.max(clip_trace))[0][0]+ib+30
                if max_deviation== (ib+30) or max_deviation == (ib+69):
                    max_deflection = np.mean(clipX[(ib+30):(ib+70)])- prior_mean
                    max_deflection_time = np.nan
                else:
                    max_deflection = clipX[max_deviation]-prior_mean
                    max_deflection_time = max_deviation - ib
            
            data_loaded['max_deflection_L'].append(-max_deflection)
            data_loaded['max_deflection_time'].append(max_deflection_time)
            
            if dat.has_key('breathe') and not np.isnan(max_deflection_time):
                br_phase_at_peak = dat['breathe']['phases'][s+max_deflection_time]
                br_rate_at_peak = BR_Rate[s+max_deflection_time]

            else:
                br_phase_at_peak = np.nan
                br_rate_at_peak = np.nan
                
    
            data_loaded['br_phase_at_peak'].append(br_phase_at_peak)
            data_loaded['br_rate_at_peak'].append(br_rate_at_peak)

            ## Perform the peak_detection (Vertical)
            thresholdV = 3*np.std(clipV[(ib-100):ib]) #Set the threshold as standard deviation movement 100ms before stim
            prior_meanV = np.mean(clipV[(ib-100):ib]) #Mean position 100ms before stim
            
            indices = peakutils.indexes(clipV, thres=thresholdV, min_dist = 10)
            if len(indices) == 0 or max(indices)<ib:
                max_deflectionV = np.nan
                max_deflection_timeV = np.nan
            else:
                first_peakV = indices[np.where(indices>ib)[0][0]]
                max_deflectionV = clipV[first_peakV] - prior_meanV
                max_deflection_timeV = first_peakV-ib
            
            data_loaded['max_deflection_V'].append(max_deflectionV)
            data_loaded['max_deflection_timeV'].append(max_deflection_timeV)
    
            
    #%%
    df = pd.DataFrame(data=data_loaded)
    df['stim_power'] = df['stim_duration']*df['stim_current']*2.1/1000
    df = df.round({'stim_power': 1,'stim_current': 0})
    df = df.set_index(['file_day','stim_duration','stim_current','stim_power'])
    
    df.to_pickle('/home/asya/Documents/data/mouse_nose/pandas_data/'+mouse_name)
    
    #df.loc['dy2'].plot(x='max_deflection_L', y='max_deflection_time', style='o')
