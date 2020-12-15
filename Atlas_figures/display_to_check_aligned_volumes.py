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


include_stacks = {'Rat_brainstem_atlas',
                  'DTA04_L',
                 #'GR35_L',
            }
stack_fixed = 'DTA04_L' #'DTA04_L'



if generate_meshes:
    generate_and_save_RV_meshes(include_stacks,stack_fixed = stack_fixed,
                                warp_setting = 24, folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')
if meshes_from_file:
    all_meshes, all_origin = load_file_RV_meshes(include_stacks,
                                                 folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')

if generate_actors:
    brain_actor_list_all_stacks = generate_whole_actors(all_meshes,all_origin,structs_to_colors=[],
                                                        include_stacks = include_stacks,
                                                        wireframe=True,opacity= 0.1)

atlas_actors = list()

if stack_fixed == 'Rat_brainstem_atlas':
    atlas_name = 'Rat_brainstem_atlas'

    structs_to_plot = ['5N','7N', '7n', 'Amb', 'LRT', 'IO']

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
        A = actor_mesh(mesh_rel2canon,color= np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255., wireframe=False,opacity = 0.05) 
        atlas_actors.append(A)


#%%
if 0: ##markers
    marker_volume_actors_all_stacks = {}
    all_brain_marker_actors_slice = list()
    for stack_moving in  include_stacks:
        if stack_moving == stack_fixed: continue
        tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m=stack_moving, warp_setting=24, 
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
        
        moving_brain_markers = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack_moving, structure='All'))
        moving_brain_markers_aligned2fixed = np.dot(R, (moving_brain_markers - om - cm).T).T + t + of + cf
        moving_brain_marker_aligned2fixed_actors = [actor_sphere((x,y,z), radius=1, color= (1,0,0)) 
                                                        for x,y,z in moving_brain_markers_aligned2fixed]
            
        marker_volume_actors_all_stacks[stack_moving] = moving_brain_marker_aligned2fixed_actors
       
if 1:
    GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
    tf_parameter_dict = load_alignment_parameters_v2(stack_f=stack_fixed, stack_m='Rat_brainstem_atlas', warp_setting=24, 
                                                 vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
    cf = np.array(tf_parameter_dict['centroid_f'])
    cm = np.array(tf_parameter_dict['centroid_m'])
    of = np.array(tf_parameter_dict['crop_origin_f'])
    om = np.array(tf_parameter_dict['crop_origin_m'])
    params = np.array(tf_parameter_dict['params'])
    Rt = np.reshape(params, (3,4))
    R = Rt[:3,:3]
    t = Rt[:3,3]
    t= np.asarray( [1,1,1])
    
    R =  np.vstack(([.55,0,0],[0,.55,0],[0,0,.55]))
    cf = cf-[0,5,5]
#    R = np.dot(R,np.vstack(([1,0,0],[0,np.cos(10),-np.sin(10)],[0,np.sin(10),np.cos(10)])))
#    R[1,1] = 0.7
#    R[2,2] = 0.6
#    of= [-380,   45,  -84]
    
    mu_new0 = np.dot(R, (GMM['mu'][0] - om - cm).T).T + t + of + cf
    mu_new1 = np.dot(R, (GMM['mu'][1] - om - cm).T).T + t + of + cf
    
    s0 = np.dot(R, ((GMM['mu'][0]+GMM['sigma0']) - om - cm).T).T + t + of + cf - mu_new0
    s1 = np.dot(R, ((GMM['mu'][1]+GMM['sigma1']) - om - cm).T).T + t + of + cf - mu_new1
    
    
    A0,_ = actor_ellipse_vtk(position = mu_new0,radius_mat = 3*s0,color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)
    A1,_ = actor_ellipse_vtk(position = mu_new1,radius_mat = 3*s1,color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)
    
    

vtk_3d_input_list = []
vtk_3d_input_list.extend([a for actor_list in brain_actor_list_all_stacks.values() for a in actor_list.values()])
vtk_3d_input_list.extend([A0,A1])
#vtk_3d_input_list.extend(atlas_actors)
#V = [a for a in marker_volume_actors_all_stacks.values()]
#vtk_3d_input_list.extend(V[0])
vtk_3d_input_list = filter(None, vtk_3d_input_list)

#vtk_3d_input_list.extend(vtk_vol_input_list)
launch_vtk(vtk_3d_input_list)
