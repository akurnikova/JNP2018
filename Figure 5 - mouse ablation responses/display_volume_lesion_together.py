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
cut_plane_origin = (0,0,-40) ## 1.4 mm lat
cut_plane_normal = (0,0,1)

#Coronall through lesions
#cut_plane_origin = (-1140,0,0)
#cut_plane_normal = (1,0,0)


#Horizontal through lesions
#cut_plane_origin = (0,160,0)
#cut_plane_normal = (0,1,0)

thickness = 6.66


include_stacks = {'DTA01_L', #+++
                  'DTA02_L', #+++
                  'DTA03_L', #---
                  'DTA04_L', #---
                  'DTA05_L', #---
                  'DTA06_L',#+++
                  'DTA07_L',#+++
                  'DTA09_L',#---
                  'DTA10_L'#+++
            }


type_by_file = {'DTA01':1,'DTA02':2,'DTA03':3,'DTA04':4,'DTA05':5,'DTA06':6,'DTA07':7,'DTA09':9,'DTA10':10,'DTA11':11}
type_by_file = {'DTA01':1,'DTA02':2,'DTA03':-3,'DTA04':-4,'DTA05':5,'DTA06':6,'DTA07':7,'DTA09':9,'DTA10':10,'DTA11':-11}

stack_fixed = 'DTA04_L'

stacks_to_colors = {'DTA01_L':(0.9, 0.6, 0.1),# +++++# orange
                    'DTA02_L': (0.9, 0.6, 0.1),# +++++# orange
                    'DTA03_L': (0.3422, 0.48, 0.6), #----# blue
                    'DTA04_L': (0.3422, 0.48, 0.6), #----# blue
                    'DTA05_L': (0.3422, 0.48, 0.6), #----# blue
                    'DTA06_L': (0.9, 0.6, 0.1),# +++++# orange
                    'DTA07_L': (0.9, 0.6, 0.1),# +++++# orange
                    'DTA09_L': (0.3422, 0.48, 0.6), #----# blue
                    'DTA10_L': (0.9, 0.6, 0.1),# +++++# orange
                    'DTA11_L': (0.3422, 0.48, 0.6), #----# blue
                    }


all_meshes, all_origin, ccentroids, volumes  = generate_lesion_meshes(include_stacks,stack_fixed = stack_fixed,#'RV16',
                                warp_setting = 24)
if generate_actors:
     whole_lesions_actors = generate_whole_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = include_stacks, wireframe=False,opacity= 0.3)

     slice_lesions_actors = generate_slice_actors_lesions(all_meshes,all_origin,stacks_to_colors,include_stacks = include_stacks, wireframe=False,\
                                                                 opacity= 0.3, cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness)

## plot atlas studd # 
#%%
atlas_actors = list()

if 1:
        atlas_meshes, atlas_origin = generate_ms_meshes()
        
        ## make actors
for name_s in atlas_meshes.keys():
        color = np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        color = (0.,0.,0.)
        A = actor_mesh(atlas_meshes[name_s],color=color , wireframe=False,opacity = 0.15) 
        atlas_actors.append(A)
       
slice_actor_dict_atlas = generate_slice_actors_atlas_ms_no_origin(atlas_meshes,cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal,thickness = thickness,\
                                opacity= 0.1, constcolor = (0.5,0.5,0.5))


    
if 1:
    GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
    tf_parameter_dict = load_alignment_parameters_v2(stack_f='DTA04_L', stack_m='Rat_brainstem_atlas', warp_setting=24, 
                                                 vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
    cf = np.array(tf_parameter_dict['centroid_f'])
    cm = np.array(tf_parameter_dict['centroid_m'])
    of = np.array(tf_parameter_dict['crop_origin_f'])
    om = np.array(tf_parameter_dict['crop_origin_m'])
    params = np.array(tf_parameter_dict['params'])
        
    G_movingvol2fixedvol = consolidate(params=params, centroid_m=cm, centroid_f=cf)
            
    Rt = np.reshape(params, (3,4))
    R = Rt[:3,:3]
    t = Rt[:3,3]        
    
    R[1,1] = 0.7
    R[2,2] = 0.6
    of= [-1260,    83,   -91]
    
    mu_new0 = np.dot(R, (GMM['mu'][0] - om - cm).T).T + t + of + cf
    mu_new1 = np.dot(R, (GMM['mu'][1] - om - cm).T).T + t + of + cf
    
    s0 = np.dot(R, ((GMM['mu'][0]+GMM['sigma0']) - om - cm).T).T + t + of + cf - mu_new0
    s1 = np.dot(R, ((GMM['mu'][1]+GMM['sigma1']) - om - cm).T).T + t + of + cf - mu_new1
    
    
    A0,_ = actor_ellipse_vtk(position = mu_new0,radius_mat = 3*s0,color=(.7, 0.4, .7),opacity = 0.3,wireframe = True)
    A1,_ = actor_ellipse_vtk(position = mu_new1,radius_mat = 3*s1,color=(.7, 0.4, .7),opacity = 0.3,wireframe = True)


    Origin = (0,-1300,  -205)# (200, 78, 200)
    yi = np.arange(-1300,-1000,1)
    xi = np.arange(0,305,1)
    zi = np.arange(-205,205,1)
    
    total_vol_P = np.zeros((len(xi),len(yi),len(zi)))
    total_vol_N = np.zeros((len(xi),len(yi),len(zi)))
    
    for stack in  ['DTA01_L','DTA02_L','DTA06_L','DTA07_L','DTA10_L']:
        V = volumes[stack]
        x_offset,y_offset,z_offset = np.round((all_origin[stack][[1,0,2]]-Origin).astype(int))
        Vn = np.zeros(total_vol_P.shape)
        Vn[x_offset:V.shape[0]+x_offset,y_offset:V.shape[1]+y_offset,z_offset:V.shape[2]+z_offset] = V
        
        #if stack in ['RFles06']:#'RFles03','RFles05',
        total_vol_P = total_vol_P+Vn
        #else:
        #    total_vol_N = total_vol_N+Vn

    vtk_3d_input_list = []
    
    #vtk_3d_input_list.extend([a['Lesion'] for a in whole_lesions_actors.values()])
    vtk_3d_input_list.extend(atlas_actors)
    vtk_3d_input_list = filter(None, vtk_3d_input_list)
            
    TV = (total_vol_P>=3).astype(int)
    total_vol = np.swapaxes(TV,1,0)
  #  total_vol2 = total_vol_P>1

    
    Origin_render = (-1300,0,  -205)
    mesh = volume_to_polydata(TV,origin = Origin_render)
    A = actor_mesh(mesh, color = (0.9, 0.6, 0.1), wireframe=False,opacity = 0.5)
    vtk_3d_input_list.append(A)
        
        
    vtk_3d_input_list.append(A)
    vtk_3d_input_list.extend([A0,A1])
    launch_vtk(vtk_3d_input_list)
    launch_vtk(vtk_3d_input_list,interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion_DTA/3D_overlap.png')
    


#%%
GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
tf_parameter_dict = load_alignment_parameters_v2(stack_f='DTA04_L', stack_m='Rat_brainstem_atlas', warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
cf = np.array(tf_parameter_dict['centroid_f'])
cm = np.array(tf_parameter_dict['centroid_m'])
of = np.array(tf_parameter_dict['crop_origin_f'])
om = np.array(tf_parameter_dict['crop_origin_m'])
params = np.array(tf_parameter_dict['params'])
    
G_movingvol2fixedvol = consolidate(params=params, centroid_m=cm, centroid_f=cf)
        
Rt = np.reshape(params, (3,4))
R = Rt[:3,:3]
t = Rt[:3,3]        

R[1,1] = 0.7
R[2,2] = 0.6
of= [-1260,    83,   -91]

mu_new0 = np.dot(R, (GMM['mu'][0] - om - cm).T).T + t + of + cf
mu_new1 = np.dot(R, (GMM['mu'][1] - om - cm).T).T + t + of + cf

s0 = np.dot(R, ((GMM['mu'][0]+GMM['sigma0']) - om - cm).T).T + t + of + cf - mu_new0
s1 = np.dot(R, ((GMM['mu'][1]+GMM['sigma1']) - om - cm).T).T + t + of + cf - mu_new1


A0,mesh0 = actor_ellipse_vtk(position = mu_new0,radius_mat = 3*s0,color=(.7, 0.4, .7),opacity = 0.5,wireframe = True)
A1,mesh1 = actor_ellipse_vtk(position = mu_new1,radius_mat = 3*s1,color=(.7, 0.4, .7),opacity = 0.5,wireframe = True)




vtk_slice_input_list = []
vtk_slice_input_list.extend([A0,A1])
vtk_slice_input_list.extend([a[1] for a in slice_lesions_actors.values()])
vtk_slice_input_list.extend([a[2] for a in slice_lesions_actors.values()])
vtk_slice_input_list.extend([a[3] for a in slice_lesions_actors.values()])

A = make_scaleBar_actor((-1200,05,-125),(0,1,0))
vtk_slice_input_list.append(A)
A = make_scaleBar_actor((-1200,05,-125),(0,0,1))
vtk_slice_input_list.append(A)

for key,a in slice_actor_dict_atlas.iteritems():
    if not key == 'Brainstem':
        vtk_slice_input_list.extend([a[1]])
        vtk_slice_input_list.extend([a[2],a[3]])


if cut_plane_normal == (1,0,0):
#    for stack,c in ccentroids.iteritems():
#        vtk_slice_input_list.extend([actor_sphere(c+[150,0,0],color = stacks_to_colors[stack],radius = 4)])
    launch_vtk(vtk_slice_input_list,init_angle = 'coronal',interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion_DTA/C',snapshot_magnification=1)
    launch_vtk(vtk_slice_input_list,init_angle = 'coronal')

if cut_plane_normal == (0,1,0):
#    for stack,c in ccentroids.iteritems():
#        vtk_slice_input_list.extend([actor_sphere(c+[0,-50,0],color = stacks_to_colors[stack],radius = 4)]) 
    launch_vtk(vtk_slice_input_list,init_angle = 'horizontal_topDown',interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion_DTA/H')
    launch_vtk(vtk_slice_input_list,init_angle = 'horizontal_topDown')
if cut_plane_normal == (0,0,1):
#    for stack,c in ccentroids.iteritems():
#        vtk_slice_input_list.extend([actor_sphere(c+[0,0,+50],color = stacks_to_colors[stack],radius = 4)]) #Sagittal
    launch_vtk(vtk_slice_input_list,init_angle = 'sagittal',interactive=False, snapshot_fn='/home/asya/Documents/data/odor_lesion_DTA/S')
    launch_vtk(vtk_slice_input_list,init_angle = 'sagittal')
