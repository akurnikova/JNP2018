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
# from data_manager import *
from annotation_utilities import *
# from registration_utilities import *
# from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from byhand_alignment import *
from vis3d_utilities_stacy import *
from neurolucida_to_volume_utilities import *


#for stack in ['RV16','RV19','RV14','RV15']: #['RV4','RV13','RV12']:#
#for stack in ['RV2_R','RV2_L','RV9','RV10']:
for stack in ['RV12']:

#stack = 'RV19' #,'RV14_65hrs','RV13_64hrs','RV19_61hrs'}:
        
    all_names_in_stack = get_list_of_RV_contours(stack)
    
    ##load in contours
    contours, markers,_ = get_sided_contours(stack,downscale = 15)
    vol_bbox_dict = contours_to_volume(contours,stack)
    
    ## 2) Contours to volumes
    if contours.has_key('Brainstem'):
        brainstem_contour_to_volume(contours,stack)
    
    ## 3) Save volumes
    save_vol_bboxes(vol_bbox_dict,stack,downscale = 15)
    save_markers(markers,stack)
    


# Now plot all the stuff     
## Optional from byhand alignment:aligns to origin
#ox,oy,oz = get_new_origin(vol_bbox_dict)
ox,oy,oz = (0,0,0) #get_new_origin(vol_bbox_dict)

#%%# Make polydata actors
polydata_actor_list = []
for name_s in ['7n_L','7N_L','LRT_L']:#['LRT_R','7n_R','7N_R','IO_R']:
    if name_s == 'Brainstem':
        continue
    polydata = volume_to_polydata(vol_bbox_dict[name_s][0], num_simplify_iter=3, smooth=True,)

    xmin, _, ymin, _, zmin, _ = vol_bbox_dict[name_s][1]
    polydata_actor = actor_mesh(polydata, color=(0.,0.,0.),origin=(xmin-ox,ymin-oy,zmin-oz),opacity = 0.4)
    polydata_actor_list.append(polydata_actor)
       
if not len(markers)==0:
    for m_i in range(0,len(markers['All'])):
       nextpt = markers['All'][m_i,:]-(ox,oy,oz)
       marker_actor = actor_sphere(position= nextpt,radius=2,color=(1,0.5,.5))
       polydata_actor_list.append(marker_actor)
        

launch_vtk(polydata_actor_list)

