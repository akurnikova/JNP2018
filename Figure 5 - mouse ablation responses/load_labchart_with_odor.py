# -*- coding: utf-8 -*-
"""
Function to load odor stim data from Labchart (converted to matlab)
@author: asya
"""

import scipy.io as si
import numpy as np
sys.path.append('/home/asya/Documents/data/utilities')

from breathe_proc import util_proc_breathe
from utils_filtering import util_lowpass, util_hipass
import pickle


datraw = si.loadmat("/home/asya/Documents/data/odor_lesion_DTA/raw_dat/DIO1_X_odorday2_great.mat"); br_sign = -1; 
save_filename = '/home/asya/Documents/data/odor_lesion_DTA/dat_preproc/DTA01_odor_wk6.pckl'

doBr = 0

def parse_by_channels(dat):
    CH = dict()
    for i,st in enumerate(dat["datastart"]):
        endidx = np.ndarray.astype(dat["dataend"][i],'int')[0]-1
        stidx = np.ndarray.astype(st,'int')[0]-1
        CH[i] = dat["data"][0,stidx:endidx]
    return CH
    
CH = parse_by_channels(datraw)
fps = datraw["samplerate"]

print datraw["titles"]
print fps

#%%
if 0:
    dat = dict()
    chL = 0
    chV = 1
    chBr = 2
    chS1 = 4
    chS2 = 3
    
if 1:
    dat = dict()
    chL = 0
    chV = 1
    chS1 = 3
    chS2 = 2

if 0:
    dat = dict()
    chL = 0
    chV = 1
    chS1 = 2
    chS2 = 3
    
dat['mm_per_volt'] = 1.3

if doBr:
    dat['fpsBR'] = fps[chBr][0]
dat['fpsX'] = fps[chL][0]
dat['fpsStim'] = fps[chS1][0]


#%% now the nose motion
def get_nose_all_mm(CH,chL,chV,dat):
# Returns the lateral movement of the nose
# Centered and converted to mm
    #if dat['fpsX'] == 10000:
    #    nose_L_in_volt_raw = -CH[chL][::10] # extract line
    #    nose_V_in_volt_raw = -CH[chV][::10] # extract line
    #else:
    nose_L_in_volt_raw = -CH[chL]
    nose_V_in_volt_raw = -CH[chV]
    nose_L_in_volt = nose_L_in_volt_raw-np.median(nose_L_in_volt_raw) # center
    nose_L_in_mm = nose_L_in_volt*dat['mm_per_volt'] # conversion factor
    nose_V_in_volt = nose_V_in_volt_raw-np.median(nose_V_in_volt_raw) # center
    nose_V_in_mm = nose_V_in_volt*dat['mm_per_volt'] # conversion factor
    return nose_L_in_mm, nose_V_in_mm
        
dat["L_filt"], dat["V_filt"] = get_nose_all_mm(CH,chL,chV,dat)

#%% Process the stim and stim times

dat["stim1"] = CH[chS1]
dat["stim2"] = CH[chS2]

#%% Get the breathing stuff
if doBr:
    input_dat = dict({"breathe_ALL":br_sign*CH[chBr]})
    dat["breathe"] = util_proc_breathe(input_dat,fps[chBr][0])

#%%
print 'start save'
f = open(save_filename, 'w')
#f = open('NW32_nose_block.pckl', 'w')
#f = open('NW32_PT3A.pckl', 'w')
pickle.dump(dat,f)
f.close()
print 'done save'


#%%
# List of stuff in the labchart data structure
"""
['comtext',
 'unittextmap',
 'com',
 'firstsampleoffset',
 'blocktimes',
 'titles',
 'datastart',
 'unittext',
 'dataend',
 'samplerate',
 'rangemax',
 'data',
 'tickrate',
 'rangemin']
"""
