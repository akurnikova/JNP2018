#%load_ext autoreload
#%autoreload 2

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
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
#from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from vis3d_utilities_stacy import *
import vtk

import xml.etree.ElementTree as ET
import pandas as pd
generate_meshes = 1
meshes_from_file = 1
generate_actors = 1
do_slice = 0

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

#Saggital through lesions
#cut_plane_origin = (0,0,-47) ## 1.4 mm lat -30 to -40 = Amb
#cut_plane_normal = (0,0,1)

#Coronall through lesions
#cut_plane_origin = (-1145,0,0)
#cut_plane_normal = (1,0,0)


#Horizontal through lesions
cut_plane_origin = (0,140,0)
cut_plane_normal = (0,1,0)

slice_thickness = 6.66
thickness  = 6.66

include_stacks = {
              'GR25_L',
              'GR27_L',
              'GR29_L',
              'GR31_L',
              'GR35_L',
              'GR38_L',
              'GR42_L',
                       }

include_stacks = {'GR23_L',
              'GR24_L',
              'GR25_L',
           #   'GR26_L',
              'GR27_L',
              'GR29_L',
              'GR30_L',
              'GR31_L',
              'GR33_L',
              'GR35_L',
              'GR37_L',
              'GR38_L',
              'GR41_L',
              'GR42_L',
                       }
stack_fixed = 'GR35_L'

'''
#list_coords = list_coords.set_index('ms')
list_coords2 = list_coords.loc[['GR25_L','GR27_L','GR29_L','GR31_L','GR35_L','GR38_L','GR42_L'],:]
df_sort = list_coords2.sort_values(by = 'cx', axis=0, ascending=True)
pal = sns.diverging_palette(10, 220, s=85,l=60, n=9)
pal_rev = np.vstack((pal[0:4],pal[9:18][::-1]))
df_sort['color'] = pal#pal_rev.tolist()
'''


cLesion = (1.,0.8,0.)
cLesionIRT = (0.95,  0.35,  0.5)
cLess = (0.6,0.3,0.17)#(0.56,0.3,0.05)
cLeast = (0.1,0.1,0.1)
cControl = 	(0.2, 0.4, 0.8)

stack_to_color = {'GluReaCh23':cLeast, ###+-+-
               'GluReaCh24':cLeast, ###+-+-
               'GluReaCh25':cLesion, ###++++
               'GluReaCh26':cLeast, ###????
               'GluReaCh27':cLesion, ###+-+-
               'GluReaCh28':cControl, ###NO LESION
               'GluReaCh29':cLess,
               'GluReaCh30':cLeast, 
               'GluReaCh31':cLesionIRT, ###++++
               'GluReaCh32':cControl, ###NO LESION
               'GluReaCh33':cLeast, 
               'GluReaCh34':cControl, ###NO LESION
               'GluReaCh35':cLesion, ###++++
               'GluReaCh36':cControl, ###NO LESION
               'GluReaCh37':cLeast, 
               'GluReaCh38':cLesionIRT, ###++++
               'GluReaCh39':cControl, ###NO LESION
               'GluReaCh41':cLeast, 
               'GluReaCh42':cLess, ###+-+-
                       }

scaling_factor = 15./1000

volume_actor_list = []
#%%
atlas_slice_dict = dict()
if 1:
    atlas_meshes = generate_ms_meshes_0origin(stack_fixed)
for name_s in atlas_meshes.keys():
        contour_points,order, P = polydata_cut(atlas_meshes[name_s], origin=(0,0,0),
                                   cut_plane_origin = cut_plane_origin, 
                                   cut_plane_normal = cut_plane_normal)
        atlas_slice_dict[name_s] = contour_points[order]*scaling_factor
        A = actor_mesh(atlas_meshes[name_s],color = (0.,0.,0,),opacity = 0.1)
        volume_actor_list.append(A)

#%%

#%%##
marker_volume_actors_all_stacks = {}
all_brain_marker_actors_slice = list()
marker_actors_all_stacks = []
for stack in include_stacks:
        plane = vtk.vtkPlane()
        tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack, warp_setting=24, 
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
    
        moving_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack, structure='All'))
        moving_brain_markers_aligned2fixed = np.dot(R, (moving_brain_markers - om - cm).T).T + t + of + cf
        
        for x,y,z in moving_brain_markers_aligned2fixed:
            dist = np.abs(plane.DistanceToPlane((x,y,z),cut_plane_normal,cut_plane_origin))
            if dist < slice_thickness:
                all_brain_marker_actors_slice.append(np.asarray(x,y,z)*scaling_factor)
    
        
        marker_volume_actors_all_stacks[stack] = moving_brain_markers_aligned2fixed*scaling_factor
        
        #col = df_sort.loc[stack]['color'][0:3]
        col = stack_to_color['GluReaCh'+stack[2:4]]
        marker_actors_all_stacks.extend([actor_sphere((x,y,z),radius = 0.5,color = col) for x,y,z in moving_brain_markers_aligned2fixed])
       
#V = [a for a in marker_volume_actors_all_stacks['GR26'].values()]
#vtk_vol_input_list.extend(V)
#launch_vtk(vtk_vol_input_list)


#%%
new_volume_list = volume_actor_list
from mpl_toolkits.mplot3d import Axes3D
fig3d = plt.figure()
ax = fig3d.add_subplot(111, projection='3d')
for stack in marker_volume_actors_all_stacks.keys():

    CX,CY,CZ  = np.median(marker_volume_actors_all_stacks[stack],0)
    #col = df_sort.loc[stack]['color'][0:3]
    col = stack_to_color['GluReaCh'+stack[2:4]]
    ax.plot([CX], [-CY],[CZ], color=col, marker = 'o', markersize = 5)
    new_volume_list.append(actor_sphere((CX/scaling_factor,CY/scaling_factor,CZ/scaling_factor),color = col,radius = 3))
ax.axis('equal')

new_volume_list.extend([a for a in marker_actors_all_stacks])
launch_vtk(new_volume_list)

#%%
