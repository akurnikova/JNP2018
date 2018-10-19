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

generate_meshes = 1
meshes_from_file = 1
generate_actors = 1
do_slice = 0


#Saggital through lesions
cut_plane_origin = (0,0,97) ## 93 1.4 mm lat
cut_plane_normal = (0,0,1)

#Coronall through lesions
cut_plane_origin = (-45,0,0)
cut_plane_normal = (1,0,0)


#Horizontal through lesions
#cut_plane_origin = (0,5,0)
#cut_plane_normal = (0,1,0)

thickness = 13.3


include_stacks = {'RFles02',
                  'RFles03',
                  'RFles04',
                  'RFles05',
                  'RFles06',
                  'RFles07',
                  'RFles08',
                  'BIRTles07',
                  'BIRTles08',
                  'BIRTles09'
                  }
stack_fixed = 'Rat_brainstem_atlas'

stacks_to_colors = {'RFles02': (0.3422, 0.48, 0.6), #----# blue
                    'RFles03': (0.9, 0.6, 0.1),# +++++# orange
                    'RFles04': (0.3422, 0.48, 0.6), #----# blue
                    'RFles05': (0.9, 0.6, 0.1),# +++++# orange
                    'RFles06':  (0.9, 0.6, 0.1),# +++++# orange
                    'RFles07': (0.3422, 0.48, 0.6), #----# blue
                    'RFles08': (0.3422, 0.48, 0.6), #----# blue
                    'BIRTles07': (0.395, 0.587, 0.389), # green
                    'BIRTles08': (0.395, 0.587, 0.389), # green
                    'BIRTles09': (0.395, 0.587, 0.389), # green
                    'BIRTles04': (0.395, 0.587, 0.389), # green
                    }


#'RFles02.':-2,'RFles03.':3,'RFles04_':-4,'RFles05_':5,'RFles06_':-6,'RFles07_':7,'RFles08_':-8


all_meshes, all_origin, ccentroids, volumes  = generate_lesion_meshes(include_stacks,stack_fixed = stack_fixed,#'RV16',
                                warp_setting = 24)
if generate_actors:
     whole_lesions_actors = generate_whole_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = include_stacks, wireframe=False,opacity= 0.1)

     slice_lesions_actors = generate_slice_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = include_stacks, wireframe=False,\
                                                                 opacity= 0.5, cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness)

## plot centroids # 

atlas_actors = list()

if stack_fixed == 'Rat_brainstem_atlas':
    atlas_name = 'Rat_brainstem_atlas'

    structs_to_plot = ['5N','7N', '7n', 'Amb', 'LRT', 'IO','Brainstem','Pr5','SpVO','SpVI','SpVC']#, 'RN', 'RtTg', 'fr']

    structs_to_plot_sided = list()
    for n in structs_to_plot:
        if n in singular_structures:
            structs_to_plot_sided.append(n)
        else:
            structs_to_plot_sided.append(convert_to_left_name(n))
            structs_to_plot_sided.append(convert_to_right_name(n))
        
    for name_s in structs_to_plot_sided:
        fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
        mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
        color = np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        color = (0.,0.,0.)
        A = actor_mesh(mesh_rel2canon,color=color , wireframe=False,opacity = 0.05) 
        atlas_actors.append(A)
        
    slice_actor_dict_atlas = generate_slice_actors_atlas_brainstem(cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness,\
                                opacity= 0.1, constcolor = (0.5,0.5,0.5))

#%%
if 1:
    GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
    A0,_ = actor_ellipse_vtk(position = (GMM['mu'][0][0],GMM['mu'][0][1],GMM['mu'][0][2]),radius_mat = (3*GMM['sigma0'][0],3*GMM['sigma0'][1],3*GMM['sigma0'][2]),color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)
    A1,_ = actor_ellipse_vtk(position = (GMM['mu'][1][0],GMM['mu'][1][1],GMM['mu'][1][2]),radius_mat = (3*GMM['sigma1'][0],3*GMM['sigma1'][1],3*GMM['sigma1'][2]),color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)


    Origin = (-205, -205, -205)# (200, 78, 200)
    xi = np.arange(-205,205,1)
    yi = np.arange(-205,205,1)
    zi = np.arange(-205,205,1)
    
    total_vol_P = np.zeros((len(xi),len(yi),len(zi)))
    total_vol_N = np.zeros((len(xi),len(yi),len(zi)))
    
    for stack in ['RFles05','RFles03','RFles06']:
        V = volumes[stack]
        x_offset,y_offset,z_offset = np.round((all_origin[stack][[1,0,2]]-Origin).astype(int))
        Vn = np.zeros(total_vol_P.shape)
        Vn[x_offset:V.shape[0]+x_offset,y_offset:V.shape[1]+y_offset,z_offset:V.shape[2]+z_offset] = V
        
        #if stack in ['RFles06']:#'RFles03','RFles05',
        total_vol_P = total_vol_P+Vn
        #else:
        #    total_vol_N = total_vol_N+Vn

    vtk_3d_input_list = []
    vtk_3d_input_list = [A0,A1]
    
    #vtk_3d_input_list.extend([a['Lesion'] for a in whole_lesions_actors.values()])
    vtk_3d_input_list.extend(atlas_actors)
    vtk_3d_input_list = filter(None, vtk_3d_input_list)
            
    TV = (total_vol_P>=3).astype(int)
    total_vol = np.swapaxes(TV,1,0)
  #  total_vol2 = total_vol_P>1

    
    Origin_render = (-205,-205,-205)
    mesh = volume_to_polydata(TV,origin = Origin_render)
    A = actor_mesh(mesh, color = (0.9, 0.6, 0.1), wireframe=False,opacity = 0.5)
    vtk_3d_input_list.append(A)
        
        
    vtk_3d_input_list.append(A)
    vtk_3d_input_list.extend([A0,A1])
    launch_vtk(vtk_3d_input_list)
    launch_vtk(vtk_3d_input_list,interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion/3D_overlap')
    




#%%
'''
GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
A0,_ = actor_ellipse_vtk(position = (GMM['mu'][0][0],GMM['mu'][0][1],GMM['mu'][0][2]),radius_mat = (3*GMM['sigma0'][0],3*GMM['sigma0'][1],3*GMM['sigma0'][2]),color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)
A1,_ = actor_ellipse_vtk(position = (GMM['mu'][1][0],GMM['mu'][1][1],GMM['mu'][1][2]),radius_mat = (3*GMM['sigma1'][0],3*GMM['sigma1'][1],3*GMM['sigma1'][2]),color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)


vtk_slice_input_list = []

vtk_slice_input_list.extend([A0,A1])
vtk_slice_input_list.extend([a[1] for a in slice_lesions_actors.values()])
vtk_slice_input_list.extend([a[2] for a in slice_lesions_actors.values()])
vtk_slice_input_list.extend([a[3] for a in slice_lesions_actors.values()])

for key,a in slice_actor_dict_atlas.iteritems():
    vtk_slice_input_list.extend([a[1]])
    if not key == 'Brainstem':
        vtk_slice_input_list.extend([a[2],a[3]])


if cut_plane_normal == (1,0,0):
#    for stack,c in ccentroids.iteritems():
#        vtk_slice_input_list.extend([actor_sphere(c+[150,0,0],color = stacks_to_colors[stack],radius = 4)])
    launch_vtk(vtk_slice_input_list,init_angle = 'coronal',interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion/C.png',snapshot_magnification=1)
    launch_vtk(vtk_slice_input_list,init_angle = 'coronal')

if cut_plane_normal == (0,1,0):
#    for stack,c in ccentroids.iteritems():
#        vtk_slice_input_list.extend([actor_sphere(c+[0,-50,0],color = stacks_to_colors[stack],radius = 4)]) 
    launch_vtk(vtk_slice_input_list,init_angle = 'horizontal_bottomUp',interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion/H.png',snapshot_magnification=1)
    launch_vtk(vtk_slice_input_list,init_angle = 'horizontal_bottomUp')
if cut_plane_normal == (0,0,1):
#    for stack,c in ccentroids.iteritems():
#        vtk_slice_input_list.extend([actor_sphere(c+[0,0,+50],color = stacks_to_colors[stack],radius = 4)]) #Sagittal
    launch_vtk(vtk_slice_input_list,init_angle = 'sagittal',interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion/S.png',snapshot_magnification=1)
    launch_vtk(vtk_slice_input_list,init_angle = 'sagittal')
'''
