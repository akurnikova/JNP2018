#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 11:36:58 2017

@author: yuncong
rearranged functions by stacy
"""

import numpy as np
import sys
import os
import bloscpack as bp
sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from metadata import *
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
from vis3d_utilities import *

from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *



DTA_R_stacks = {'DTA03_R',  'DTA06_R',  'DTA05_R', 'DTA10_R','DTA04_R'}
DTA_L_stacks = {'DTA01_L', 'DTA05_L', 'DTA10_L', 'DTA02_L', 'DTA04_L', 'DTA07_L',  'DTA03_L',  'DTA06_L',  'DTA09_L',  'DTA11_L'}

GR_stacks = {'GR23_L', 'GR26_L', 'GR30_L', 'GR35_L', 'GR40_L',
'GR24_L', 'GR27_L', 'GR31_L', 'GR37_L', 'GR41_L',
'GR25_L', 'GR29_L', 'GR33_L', 'GR38_L', 'GR42_L'}

RFles_stacks = {'RFles02', 'RFles05', 'RFles08',  'RFles03', 'RFles06','RFles04', 'RFles07'}

BIRT_les_stacks = {'BIRTles04', 'BIRTles09', 'BIRTles07', 'BIRTles08'}

for stack in DTA_L_stacks:
    cut_orientation = 's'
    
    ## 1) Load in contours
    contours, markers = get_stacy_contours(stack, downscale = 15, cut_orientation = cut_orientation)
    vol_bbox_dict = contours_to_volume(contours,stack,interpolation_direction='z')
    
    ## 2) Contours to volumes
    #if contours.has_key('Brainstem'):
    #    brainstem_contour_to_volume(contours,stack)
    
    ## 3) Save volumes
    save_vol_bboxes(vol_bbox_dict,stack,downscale = 15)
    save_markers(markers,stack)


#%%# Make polydata actors
ox,oy,oz = (0,0,0)
polydata_actor_list = []
for name_s in contours.keys():#['LRT_R','7n_R','7N_R','IO_R']:
    if name_s == 'Brainstem':
        continue
    polydata = volume_to_polydata(vol_bbox_dict[name_s][0], num_simplify_iter=3, smooth=True,)

    xmin, _, ymin, _, zmin, _ = vol_bbox_dict[name_s][1]
    polydata_actor = actor_mesh(polydata, color=(0.,0.,0.),origin=(xmin-ox,ymin-oy,zmin-oz),opacity = 0.4)
    polydata_actor_list.append(polydata_actor)        

launch_vtk(polydata_actor_list)