#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 15:06:19 2018

@author: asya
"""
import sys
import os
#import time

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from metadata import *
# from data_manager import *
from annotation_utilities import *
# from registration_utilities import *
# from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from byhand_alignment import *
from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *

from matplotlib import pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
import seaborn as sns
import matplotlib.patches as mpatches

    
fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/FN_1_S.pdf'

include_stacks = {'RV15','RV16'}

#include_stacks = {'RV7','RV8','RV12'}
include_stacks = {'RV4','RV14','RV13','RV19','RV9','RV10'}

#include_stacks = {'RV12'}
dim1 = 0
dim2 = 1

#dim1 = 0, dim2 = 1

stack_to_color = dict()
for i,stack in enumerate(include_stacks):
    stack_to_color[stack] = sns.color_palette()[i]

fig = figure()
ax = fig.add_subplot(111)

for stack in include_stacks:
    C,_,D = get_sided_contours(stack)
    
    tf_parameter_dict = load_alignment_parameters_v2(stack_f='Rat_brainstem_atlas', stack_m=stack, warp_setting=24, 
                                                     vol_type_f='annotationAsScore', vol_type_m='annotationAsScore',
                                                     downscale=15)
    cf = tf_parameter_dict['centroid_f']
    cm = tf_parameter_dict['centroid_m']
    of = tf_parameter_dict['crop_origin_f']
    om = tf_parameter_dict['crop_origin_m']
    params = tf_parameter_dict['params']
    Rt = np.reshape(params, (3,4))
    R = Rt[:3,:3]
    t = Rt[:3,3]
            
    
    for d in range(0,len(D)):
        if len(D[d]) == 0: continue
        dend = D[d][0]
        dend_aligned2fixed = np.dot(R, (dend - om - cm).T).T + t + of + cf
        plt.plot(dend_aligned2fixed[:,dim1],-dend_aligned2fixed[:,dim2],color = stack_to_color[stack],linewidth = 0.5)
        
    for c in range(0,len(C['CellBody_L'])):
        cell = C['CellBody_L'][c]
        cell_closed = vstack((cell,cell[0]))
        cell_aligned2fixed = np.dot(R, (cell_closed - om - cm).T).T + t + of + cf
        
        ring_mixed = Polygon(zip(cell_aligned2fixed[:,dim1],-cell_aligned2fixed[:,dim2]))
        ring_patch = PolygonPatch(ring_mixed)
        ring_patch.set_color(stack_to_color[stack])
        ax.add_patch(ring_patch)
   
'''     
    for c in range(0,len(C['7N_L'])):
        cell = C['7N_L'][c]
        cell_closed = vstack((cell,cell[0]))
        cell_aligned2fixed = np.dot(R, (cell_closed - om - cm).T).T + t + of + cf
        
        ring_mixed = Polygon(zip(cell_aligned2fixed[:,dim1],cell_aligned2fixed[:,dim2]))
        ring_patch = PolygonPatch(ring_mixed)
        ax.add_patch(ring_patch)
'''

ax.set_aspect('equal')


fp = DataManager.get_mesh_filepath(stack_m='Rat_brainstem_atlas', structure='7N_L')
mesh = load_mesh_stl(fp,return_polydata_only=True)

if dim1 == 0 and dim2 == 1: #Sagittal
    for orig in np.arange(80,200,10):
        contour_points,order, _ = polydata_cut(mesh, origin=(0,0,0),
                                           cut_plane_origin = (0,0,orig), 
                                           cut_plane_normal = (0,0,1))
        if len(order) ==0: continue
        print orig
        FN = contour_points[order]
        FN = vstack((FN,FN[0]))
        plt.plot(FN[:,dim1],-FN[:,dim2],'k')
        
        
if dim1 == 2 and dim2 == 1: #Coronal
    for orig in np.arange(-100,100,10):
        contour_points,order, _ = polydata_cut(mesh, origin=(0,0,0),
                                           cut_plane_origin = (orig,0,0), 
                                           cut_plane_normal = (1,0,0))
        if len(order) ==0: continue
        print orig
        FN = contour_points[order]
        FN = vstack((FN,FN[0]))
        plt.plot(FN[:,dim1],-FN[:,dim2],'k')
        
## create legend for the densities
patch_handles = []
for stack in include_stacks:
     patch_handles.append(mpatches.Patch(color=stack_to_color[stack], label='Density %s'%stack))
plt.legend(handles=patch_handles)

plt.savefig(fig_title, format='pdf',dpi=fig.dpi)