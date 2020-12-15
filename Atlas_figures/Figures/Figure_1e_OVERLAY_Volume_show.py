#%load_ext autoreload
#%autoreload 2

import sys
import os
#import time

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
import matplotlib.patches as mpatches


sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
sys.path.append('/home/asya/Documents/data/utilities')

from utilities2015 import *
from metadata import *
from data_manager import *
from annotation_utilities import *
from registration_utilities import *
#from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from vis3d_utilities_stacy import *
from utils_for_plots import *
import vtk

import xml.etree.ElementTree as ET

from skimage import measure

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import seaborn as sns


thickness_mm = 0.05
thickness = thickness_mm*1000./15.



## Coronal through RF
cut_plane_origin_mm = (-1.,0,0)
cut_plane_normal = (1,0,0)

spacefill = True

cut_plane_origin = np.asarray(cut_plane_origin_mm)*1000./15.

slice_actor_dict_atlas_cor = generate_slice_actors_atlas(cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, thickness = thickness,\
                            doBrainstem = 1,doCortex = 0,opacity= 0.3, structs_to_plot = all_known_structures_sided, constcolor = [],spacefill = spacefill,atlas_name = atlas_name)




## Saggital through PreBotC
cut_plane_origin_mm = (0.,0.,2.15)
cut_plane_normal = (0,0,1)

cut_plane_origin = np.asarray(cut_plane_origin_mm)*1000./15.

slice_actor_dict_atlas_sag = generate_slice_actors_atlas(cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, thickness = thickness,\
                            doBrainstem = 1,doCortex = 0,opacity= 0.3, structs_to_plot = all_known_structures_sided, constcolor = [],spacefill = spacefill,atlas_name = atlas_name)


structs_to_plot = ['Brainstem','7N','IO','LRT','Amb','5N','SpVI','SpVO','SpVC','7n','PreBotC']

structs_to_plot_sided = list()
for n in structs_to_plot:
    if n in singular_structures:
        structs_to_plot_sided.append(n)
    else:
        structs_to_plot_sided.append(convert_to_left_name(n))
        structs_to_plot_sided.append(convert_to_right_name(n))
structure_mesh_actors_rel2canon = list()

for name_s in structs_to_plot_sided:
    fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
    fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=name_s)
    if os.path.isfile(fp_spacefill) and spacefill:
        mesh_rel2canon = load_mesh_stl(fp_spacefill,return_polydata_only=True)      
    else:
        mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
    A = actor_mesh(mesh_rel2canon,np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255., wireframe=False,opacity = 0.05) 
    structure_mesh_actors_rel2canon.append(A)
    
 #   structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas_cor[name_s][1])
 #   structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas_cor[name_s][2])
 #   structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas_cor[name_s][3])
    
    structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas_sag[name_s][1])
    structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas_sag[name_s][2])
    structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas_sag[name_s][3])

#%%
launch_vtk([] \
        + structure_mesh_actors_rel2canon,
               init_angle='horizontal_topDown',
        background_color=(1.0,1.0,1.0),
        interactive=True)#, snapshot_fn='/home/asya/Documents/Yuncong_code/Figures_code/NIRT_snap1')
