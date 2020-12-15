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
#cut_plane_origin = (0,0,140)
#cut_plane_normal = (0,0,1)

#Saggital through lesions
cut_plane_origin = (0,0,50)
cut_plane_normal = (0,0,1)
thickness = 20


include_stacks = {'RV2_L'} #{'RFles02','RFles03','RFles04','RFles05','RFles06','RFles07','RFles08'} 
stack_fixed = 'Rat_brainstem_atlas'

structures_to_display = alignment_structures_sided
'''
brain_atlas = load_original_volume_all_known_structures_v2(stack='Rat_brainstem_atlas', 
                                                                sided=True, 
                                                               include_surround=False, 
                                                               common_shape=False,
                                                               volume_type='annotationAsScore')  
origins = dict()
mesh_dicts = dict()
for name_s in brain_atlas.iterkeys():
    if name_s not in structures_to_display:
        continue
    origins[name_s] = brain_atlas[name_s][1][[0,2,4]]
    mesh_dicts[name_s] = volume_to_polydata(brain_atlas[name_s][0], num_simplify_iter= 3, smooth=True) 
'''
        
if generate_meshes:
    generate_and_save_RV_meshes(include_stacks,stack_fixed = stack_fixed,#'RV16',
                                warp_setting = 24,
                                folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')
all_origin = dict()
all_meshes = dict()

if meshes_from_file:
    all_meshes, all_origin = load_file_RV_meshes(include_stacks,folder = '/home/asya/Documents/Yuncong_code/data/meshes_save/')

#all_meshes['Rat_brainstem_atlas'] = mesh_dicts
#all_origin['Rat_brainstem_atlas'] = origins

if generate_actors:
    brain_actor_list_all_stacks = generate_whole_actors(all_meshes,all_origin,structs_to_colors=[], include_stacks = include_stacks,wireframe=True,opacity= 0.1)



#%%
atlas_actors = list()
if stack_fixed == 'Rat_brainstem_atlas':
    atlas_name = 'Rat_brainstem_atlas'

    structs_to_plot = ['5N','7N', '7n', 'Amb', 'LRT', 'IO']#, 'RN', 'RtTg', 'fr']

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
        A = actor_mesh(mesh_rel2canon,color= np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255., wireframe=False,opacity = 0.2) 
        atlas_actors.append(A)


volume_marker_actors_all_stacks =  generate_RV_marker_actors_for_atlas(include_stacks=include_stacks, radius=2, stacks_to_colors = [])

#%%
vtk_3d_input_list = []
vtk_3d_input_list.extend([a for actor_list in brain_actor_list_all_stacks.values() for a in actor_list.values()])

vtk_3d_input_list.extend([a for a in volume_marker_actors_all_stacks])

vtk_3d_input_list.extend(atlas_actors)
vtk_3d_input_list = filter(None, vtk_3d_input_list)

launch_vtk(vtk_3d_input_list)
