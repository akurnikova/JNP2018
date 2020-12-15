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


include_stacks = {'RV13'}
atlas_name = 'Rat_full_atlas'

structs_to_plot = ['Forel','PR','fr']
#structs_to_plot = ['SpVO','SpVI','Pr5']
#structs_to_plot = ['7N', 'LRT']


structs_to_plot_sided = list()
for n in structs_to_plot:
    if n in singular_structures:
        structs_to_plot_sided.append(n)
    else:
        structs_to_plot_sided.append(convert_to_left_name(n))
        structs_to_plot_sided.append(convert_to_right_name(n))
        
structure_mesh_actors_rel2canon = list()

volume_marker_actors_all_stacks =  generate_RV_marker_actors_for_atlas(include_stacks=include_stacks, radius=2, stacks_to_colors = [])
structure_mesh_actors_rel2canon.extend([a for a in volume_marker_actors_all_stacks])

#%%
for name_s in structs_to_plot_sided:
    fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
    mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
    A = actor_mesh(mesh_rel2canon,np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255., wireframe=False,opacity = 0.2) 
    structure_mesh_actors_rel2canon.append(A)

#P = actor_volume_respace(np.float32(density[name_stack]/np.max(density[name_stack])),'probability',spacing = [xstep,ystep,zstep] ,origin =  [xmin, ymin, zmin])
#P = actor_volume_respace(np.float32(density_vol[name_stack]),'score',spacing = [xstep,ystep,zstep] ,origin =  [xmin, ymin, zmin])

launch_vtk(structure_mesh_actors_rel2canon,
           init_angle='sagittal',
    background_color=(1,1,1))
