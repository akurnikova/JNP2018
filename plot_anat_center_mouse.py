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

include_stacks = {'GR23_L',
              'GR24_L',
              'GR25_L',
              'GR26_L',
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

#include_stacks = {'GR27_L'}

stack_fixed = 'GR35_L'

cLesion = (1.,0.8,0.)
cLess = (0.8,0.5,0.24)#(0.6,0.3,0.17)#(0.56,0.3,0.05)
cLeast = (0.1,0.1,0.1)
cControl = 	(0.2, 0.4, 0.8)

stack_to_color = {'GluReaCh23':cLeast, ###+-+-
                  'GluReaCh24':cLess, ###+-+-
                  'GluReaCh25':cLesion, ###++++
                  'GluReaCh26':cLeast, ###????
                  'GluReaCh27':cLesion, ###+-+-
                  'GluReaCh28':cControl, ###NO LESION
                  'GluReaCh29':cLess,
                  'GluReaCh30':cLeast, 
                  'GluReaCh31':cLesion, ###++++
                  'GluReaCh32':cControl, ###NO LESION
                  'GluReaCh33':cLeast, 
                  'GluReaCh34':cControl, ###NO LESION
                  'GluReaCh35':cLesion, ###++++
                  'GluReaCh36':cControl, ###NO LESION
                  'GluReaCh37':cLeast, 
                  'GluReaCh38':cLesion, ###++++
                  'GluReaCh39':cControl, ###NO LESION
                  'GluReaCh41':cLeast, 
                  'GluReaCh42':cLesion, ###+-+-
                       }

scaling_factor = 15./1000

volume_actor_list = list()
#%%
atlas_slice_dict = dict()
if 1:
   atlas_meshes = generate_ms_meshes_0origin(stack_fixed)
for name_s in atlas_meshes.keys():
        contour_points,order, P = polydata_cut(atlas_meshes[name_s], origin=(0,0,0),
                                   cut_plane_origin = cut_plane_origin, 
                                   cut_plane_normal = cut_plane_normal)
        atlas_slice_dict[name_s] = contour_points[order]*scaling_factor
        A = actor_mesh(atlas_meshes[name_s],color = (0.,0.,0,),opacity = 0.01)
        volume_actor_list.append(A)


def get_markers_within_shape(mesh,all_markers):
    thickness = 5
    cut_plane_normal = (1,.0,.0)
    markers_in_struct = list()
    for OZ in range(-400,-200,10):
        cut_plane_origin = (OZ,0.,0.)

        plane=vtk.vtkPlane()
    
        ## Load in slice
        contour_points1,order1, _ = polydata_cut(mesh, origin=(0,0,0),
                                   cut_plane_origin = np.asarray(cut_plane_origin)+thickness*np.asarray(cut_plane_normal), 
                                   cut_plane_normal = cut_plane_normal)
        slice_coords_top = contour_points1[order1]
        
        contour_points2,order2, _ = polydata_cut(mesh, origin=(0,0,0),
                                   cut_plane_origin = np.asarray(cut_plane_origin)-thickness*np.asarray(cut_plane_normal), 
                                   cut_plane_normal = cut_plane_normal)
        slice_coords_bot = contour_points2[order2]
        
        polygon_T = Polygon(zip([0,0,0],[0,0,0]))
        polygon_B = Polygon(zip([0,0,0],[0,0,0]))
        
        if len(order1) > 0:    
            X = slice_coords_top[:,1]
            Y = slice_coords_top[:,2]
            polygon_T = Polygon(zip(X,Y))

        if len(order2) > 0:    
            X = slice_coords_bot[:,1]
            Y = slice_coords_bot[:,2]
            polygon_B = Polygon(zip(X,Y))
        
        if len(order1) == 0 and len(order2) == 0:
            continue
        
        ## Count cells
        for i, (x,y,z) in enumerate(all_markers):
            dist = np.abs(plane.DistanceToPlane((x,y,z),cut_plane_normal,cut_plane_origin))
            if dist > thickness: continue
            p = Point(y,z)
            T = p.within(polygon_T)
            B = p.within(polygon_B)
            if T or B:
                markers_in_struct.append([x,y,z])
    return markers_in_struct

#%%
density_meshes = dict()
if 1:
    #fig = plt.figure(figsize=(10,10))
    #ax = fig.add_subplot(111)
    
    density_slice_dict = {}
    for stack in include_stacks:    
    
        P = dict()
        fp = DataManager.get_density_pvol_filepath(stack,0.266)
        P[stack] = pickle.load(open(fp,'rb'))
    
        
        maxcol = np.max(P[stack]['vol'])
        V = P[stack]['vol']/np.sum(P[stack]['vol'])
        sorted_data = np.sort(np.ravel(V))[::-1]
        percentile_d = sorted_data[np.where(np.cumsum(sorted_data)>float(0.05))[0][0]]
        density_vol = V>percentile_d
        density_vol = np.swapaxes(density_vol,0,1)
            
        mesh = volume_to_polydata(density_vol,origin = P[stack]['origin'])
        #A = actor_mesh(mesh,color = stack_to_color['GluReaCh'+stack[2:4]],opacity = 0.3)
        #volume_actor_list.append(A)
        
        contour_points,order, P = polydata_cut(mesh, origin=(0,0,0),
                                       cut_plane_origin = cut_plane_origin, 
                                       cut_plane_normal = cut_plane_normal)
        density_slice_dict[stack] = contour_points[order]*scaling_factor
        density_meshes[stack] = mesh

#%%##
markers_all_stacks = {}
all_brain_marker_actors_slice = list()
for stack in include_stacks:
        if stack == stack_fixed:
            moving_brain_markers_aligned2fixed = bp.unpack_ndarray_file(get_stacy_markers_filepath(stack=stack, structure='All'))
        else:
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
    
        
        markers_all_stacks[stack] = moving_brain_markers_aligned2fixed*scaling_factor
       # volume_actor_list.extend([actor_sphere((x,y,z),radius = 0.5,color = stack_to_color['GluReaCh'+stack[2:4]]) for x,y,z in moving_brain_markers_aligned2fixed])


#%%
if 0:
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    
    for stack in {'GR23_L',
                  'GR24_L',
                  'GR25_L',
                  'GR26_L',
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
                           }:
            if len(density_slice_dict[stack])==0:
                continue
            if cut_plane_normal == (0,0,1):
                X = density_slice_dict[stack][:,0]
                Y = -density_slice_dict[stack][:,1]
            if cut_plane_normal == (1,0,0):
                X = density_slice_dict[stack][:,2]
                Y = -density_slice_dict[stack][:,1]
            if cut_plane_normal == (0,1,0):
                X = density_slice_dict[stack][:,2]
                Y = density_slice_dict[stack][:,0]
        
                
            col = stack_to_color['GluReaCh'+stack[2:4]]
    
            
            ring = LinearRing(zip(X,Y))
            x, y = ring.xy
        
            ax.plot(x, y, color=col, alpha=0.8,
                linewidth=1, solid_capstyle='round', zorder=1)
        
    for name_s in atlas_slice_dict.keys():
            if len(atlas_slice_dict[name_s])==0:
                continue
            if cut_plane_normal == (0,0,1):
                X = atlas_slice_dict[name_s][:,0]
                Y = -atlas_slice_dict[name_s][:,1]
            if cut_plane_normal == (1,0,0):
                X = atlas_slice_dict[name_s][:,2]
                Y = -atlas_slice_dict[name_s][:,1]
            if cut_plane_normal == (0,1,0):
                X = atlas_slice_dict[name_s][:,2]
                Y = atlas_slice_dict[name_s][:,0]
                
            col = (0.3,0.3,0.3)
    
            
            ring = LinearRing(zip(X,Y))
            x, y = ring.xy
        
            ax.plot(x, y, color=col, alpha=1,
                linewidth=1, solid_capstyle='round', zorder=0)
    
    
    for stack in markers_all_stacks.keys():
        CX,CY,CZ  = np.mean(markers_all_stacks[stack],0)
        col = stack_to_color['GluReaCh'+stack[2:4]]
        if cut_plane_normal == (0,0,1):
             ax.plot(CX, -CY, color=col, marker = 'o', markersize = 5)
        if cut_plane_normal == (1,0,0):
             ax.plot(CZ, -CY, color=col, marker = 'o', markersize = 5)
        if cut_plane_normal == (0,1,0):
             ax.plot(CZ, CX, color=col, marker = 'o', markersize = 5)
#%%
#%%
GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
tf_parameter_dict = load_alignment_parameters_v2(stack_f='GR35_L', stack_m='Rat_brainstem_atlas', warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
cf = np.array(tf_parameter_dict['centroid_f'])
cm = np.array(tf_parameter_dict['centroid_m'])
of = np.array(tf_parameter_dict['crop_origin_f'])
om = np.array(tf_parameter_dict['crop_origin_m'])
params = np.array(tf_parameter_dict['params'])
    
Rt = np.reshape(params, (3,4))
R = Rt[:3,:3]
t = Rt[:3,3]        
mu_new0 = np.dot(R, (GMM['mu'][0] - om - cm).T).T + t + of + cf
mu_new1 = np.dot(R, (GMM['mu'][1] - om - cm).T).T + t + of + cf

s0 = np.dot(R, ((GMM['mu'][0]+GMM['sigma0']) - om - cm).T).T + t + of + cf - mu_new0
s1 = np.dot(R, ((GMM['mu'][1]+GMM['sigma1']) - om - cm).T).T + t + of + cf - mu_new1


A0,mesh0 = actor_ellipse_vtk(position = mu_new0,radius_mat = 3*s0,color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)
A1,mesh1 = actor_ellipse_vtk(position = mu_new1,radius_mat = 3*s1,color=(.7, 0.4, .7),opacity = 0.2,wireframe = True)

volume_actor_list.extend([A0,A1])
#fig3d = plt.figure()
#ax = fig3d.add_subplot(111, projection='3d')
list_coords = pd.DataFrame()
for stack in markers_all_stacks.keys():
    mesh = density_meshes[stack]
    markers = get_markers_within_shape(mesh,markers_all_stacks[stack]/scaling_factor)
    markers = markers_all_stacks[stack]/scaling_factor
    markers = markers[markers[:,1] > 105,:]
    markers = markers[markers[:,1] < 160,:]
    CX,CY,CZ  = np.mean(markers,0)
    col = stack_to_color['GluReaCh'+stack[2:4]]
   # if stack in df_sort.index:
   #     col = df_sort.loc[stack]['color'][0:3]
   # else: continue
    #ax.plot([CX], [-CY],[CZ], color=col, marker = 'o', markersize = 5)
    
    print CY
    volume_actor_list.append(actor_sphere((CX,CY,CZ),color = col,radius = 3))
    
    list_coords= list_coords.append(pd.DataFrame([[stack,CX,CY,CZ]], columns = ['ms','cx','cy','cz']),ignore_index = True)
#ax.axis('equal')
list_coords = list_coords.set_index('ms')
launch_vtk(volume_actor_list,init_angle = 'Saggital')
#%%
'''
volume_actor_list2 = []
point_count0 = dict()
point_count1 = dict()
for stack in ['GR35_L']:#include_stacks:
  #  stack = 'GR35_L'
    points = vtk.vtkPoints()
    M = markers_all_stacks[stack]/scaling_factor
    for pt in M:
        points.InsertNextPoint(pt)
        volume_actor_list2.extend([actor_sphere(pt,radius = 2)])
    pointsPolydata = vtk.vtkPolyData()
    pointsPolydata.SetPoints(points)
    
    selectEnclosedPoints0 = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPoints0.SetInputData(pointsPolydata)
    selectEnclosedPoints0.SetSurfaceData(mesh0)
    selectEnclosedPoints0.Update()
    numpoints_inside0 = 0
    
    selectEnclosedPoints1 = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPoints1.SetInputData(pointsPolydata)
    selectEnclosedPoints1.SetSurfaceData(mesh1)
    selectEnclosedPoints1.Update()
    numpoints_inside1 = 0
    
    for i in range(markers_all_stacks[stack].shape[0]):
        numpoints_inside0 += selectEnclosedPoints0.IsInside(i)
        numpoints_inside1 += selectEnclosedPoints1.IsInside(i)
    
    point_count0[stack] = numpoints_inside0
    point_count1[stack] = numpoints_inside1

volume_actor_list2.extend([actor_mesh(mesh0,color = (1,0,1),opacity = 0.3),actor_mesh(mesh0,color = (1,0,1),opacity = 0.3)])
launch_vtk(volume_actor_list2)
'''