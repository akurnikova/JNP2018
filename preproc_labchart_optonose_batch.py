# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 09:27:23 2016

@author: asya
"""
from __future__ import division
import scipy.io as si
import numpy as np
import sys
sys.path.append('/home/asya/Documents/data/utilities')
from utils_filtering import util_hipass,util_lowpass
from breathe_proc_mouse import util_proc_breathe
import scipy.interpolate as scint
import os

import pickle

mouse_name = 'GluReaCh26'

input_dir = '/home/asya/Documents/data/mouse_nose/raw_data/'+mouse_name+'/'
output_dir = '/home/asya/Documents/data/mouse_nose/preproc_data/'+mouse_name+'/'


if not os.path.exists(input_dir):
    print "INPUT DIRECTORY NOT FOUND"
    quit()
    
file_list = os.listdir(input_dir)   
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def parse_by_channels(dat):
    blockI = 0
    CH = dict()
    for i,st in enumerate(dat["datastart"]):
        endidx = np.ndarray.astype(dat["dataend"][i],'int')[blockI]-1
        stidx = np.ndarray.astype(st,'int')[blockI]-1
        CH[i] = dat["data"][0,stidx:endidx]
    return CH

def round_to_base(x, base=5):
    return int(base * round(float(x)/base))

for fle in file_list:
    datraw = si.loadmat(input_dir+fle)
    save_filename = output_dir+fle
    save_filename = save_filename[:-4]
    
    CH = parse_by_channels(datraw)
    fps = datraw["samplerate"]
    
    if len(CH) == 4:
        chX = 0 # Lateral nose
        chV = 1 # Vertical nose
        chBr= 2 # Breathing
        chS = 3 # LED stim
        doBr = 1 # Flag to process breathing signal on
    else:
        chX = 0 # Lateral nose
        chV = 1 # Vertical nose
        chS = 2 # LED stim
        doBr = 0 # Flag to process breathing signal off
    ## cut 1st second
    #CH = cropsection(CH,10000)
    
    dat = dict()
    # Preprocess and downsample the nose position
    if fps[0] == 10000:
        nose_rawX = -(CH[0][::10]-np.mean(CH[0]))
        dat['noseX'] = util_lowpass(nose_rawX,1000,30,1)
        nose_rawV = -(CH[1][::10]-np.mean(CH[1]))
        dat['noseV'] = util_lowpass(nose_rawV,1000,30,1)
        dat['stim_trace'] = CH[chS][::10]
    else:
        nose_rawX = -(CH[0]-np.mean(CH[0]))
        dat['noseX'] = util_lowpass(nose_rawX,1000,30,1)
        nose_rawV = -(CH[1]-np.mean(CH[1]))
        dat['noseV'] = util_lowpass(nose_rawV,1000,30,1)
        dat['stim_trace'] = CH[chS]#[::10]
    #%%
    # Preprocess and downsample the stim triggers
    dat['stim_trig'] = np.diff(dat['stim_trace'])>0.9
    dat['stim_trig_idx'] = np.where(dat['stim_trig'])[0]
    
    dat['stim_trig_off'] = np.diff(-dat['stim_trace'])>0.9
    dat['stim_trig_idx_off'] = np.where(dat['stim_trig_off'])[0]
    
    
    stim_Volt = []
    stim_time = []
    for i,on in enumerate(dat['stim_trig_idx']):
        t_on = np.round(dat['stim_trig_idx_off'][i]-on,0)
        Volt_clip = dat['stim_trace'][(on+2):(on+t_on-2)]
        stim_time.append(t_on)
        stim_Volt.append(0.01*round_to_base(100*np.mean(Volt_clip)))
    dat['stim_time'] = stim_time
    dat['stim_Volt'] = stim_Volt
    #dat['stim_trig_idx'] = np.round(dat['stim_trig_idx']/10)
    
    #%% Breathing proc
    # Preprocess and downsample the breathing
    if doBr:
        br_raw = util_lowpass(CH[chBr],fps[chBr],40,1)
        if fps[chBr] == 10000:
            br_raw_dwn = br_raw[::10]
        else:
            br_raw_dwn = br_raw
        
        if np.mean(br_raw_dwn)<1: #Inverting amplifier not used
            br_raw_dwn = br_raw_dwn-np.mean(br_raw_dwn)
        else:
            br_raw_dwn = -(br_raw_dwn-np.mean(br_raw_dwn))
        dat_temp = dict()
        dat_temp['breathe_raw'] = util_hipass(br_raw_dwn,1000,5,1)
        dat_temp['breathe_raw'] = util_lowpass(dat_temp['breathe_raw'],1000,30,1)
        
        #%% Breathing proc
        # continue breathing preproc - hilbert
        dat['breathe'] = util_proc_breathe(dat_temp,1000)
        
        ## Interpolate breathing phase
        IndX = np.arange(0, np.size(dat['breathe']['breathe']))
        pk_s = np.tile(np.pi,np.size(dat['breathe']['peaks']))
        val_s = np.tile(0,np.size(dat['breathe']['valleys']))
        val_s1 = np.tile(2*np.pi,np.size(dat['breathe']['valleys']))
        
        P = np.vstack((dat['breathe']['peaks'],pk_s))
        V = np.vstack((dat['breathe']['valleys'],val_s))
        V1 = np.vstack((dat['breathe']['valleys']-1,val_s1))
        
        phases_raw = np.hstack((P,V,V1))
        
        phases_raw_interp = scint.interp1d(phases_raw[0],phases_raw[1] ,fill_value = np.nan,bounds_error = False)
        dat['breathe']['phases'] = phases_raw_interp(IndX)
    
    
    ## save the file
    print 'start save'
    f = open(save_filename, 'w')
    pickle.dump(dat,f)
    f.close()