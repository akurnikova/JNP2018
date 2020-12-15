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
from registration_utilities_2 import *
#from conversion import *
from vis3d_utilities import *

from volume_display_utilities import *
from vis3d_utilities_stacy import *
import vtk
import seaborn as sns

import xml.etree.ElementTree as ET
sys.path.append('/home/asya/Documents/GMM_test')
from visualization import *
from sklearn import mixture

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

from pomegranate import *

volume_actor_list = list()
d = .15

n_overlap = 4

dI = .08

include_stacks = {'RV13','RV9','RV10','RV14','RV4','RV19'}

atlas_name = 'Rat_brainstem_atlas'
individual_density_meshes = dict()
C = sns.color_palette("Set1", n_colors=8, desat=.5)
stacks_to_color = {'RV19':C[0],'RV14':C[1],'RV10':C[2],'RV4':C[3],'RV13':C[4],'RV9':C[6]}#,'RV9','RV10','RV14','RV4','RV19')
## Load in stacks - this is slow
if 1:
    P_dict = dict()
    for stack in include_stacks:
        fp = DataManager.get_density_pvol_filepath(stack,.133)#.0666)
        P = pickle.load(open(fp,'rb'))
        P_dict[stack] = P

## Create density actors for each stack
if 0:
    for stack in include_stacks:
        P = P_dict[stack]
        color = np.array(stacks_to_color[stack])
        ## Option 1 for density vol
        #density_vol = (P['vol']/np.max(P['vol'])) > 0.7
    
        ## Option 2 for density vol
        maxcol = np.max(P['vol'])
        V = P['vol']/np.sum(P['vol'])
        sorted_data = np.sort(np.ravel(V))[::-1]
        percentile_d = sorted_data[np.where(np.cumsum(sorted_data)>float(dI))[0][0]]
        density_vol = V>percentile_d
        
        imagedata = volume_to_imagedata_respace(np.float32(density_vol), origin=P['origin'], spacing = [P['xstep'],P['ystep'],P['zstep']])
        surface = vtk.vtkMarchingCubes()
        surface.SetInputData(imagedata)
        surface.ComputeNormalsOn();
        surface.SetValue(0, 0.5);
        surface.Update()
            
        mesh = surface.GetOutput()
        individual_density_meshes[stack] = mesh
      #  color = np.array([0.1,0.1,0.1])
        #A = actor_mesh(mesh, color = color, wireframe=False,opacity = 0.3)
        #A.GetProperty().SetColor(color)
        #volume_actor_list.append(A)

## Create OVERLAP volume
if 1:
    Origin = (-450, -350, -250)# (200, 78, 200)
    xi = np.arange(-450,500,3)
    yi = np.arange(-350,200,3)
    zi = np.arange(-250,250,3)
    
    total_vol = np.zeros((len(xi),len(yi),len(zi)))
    
    total_D = np.zeros((len(xi),len(yi),len(zi)))
    
    for stack in include_stacks:
        print stack
        P = P_dict[stack]
        
        maxcol = np.max(P['vol'])
        V = P['vol']/np.sum(P['vol'])
        
        sorted_data = np.sort(np.ravel(V))[::-1]
        percentile_d = sorted_data[np.where(np.cumsum(sorted_data)>float(d))[0][0]] 
        density_vol = V>percentile_d
        
        color = np.array(stacks_to_color[stack])
        
        x_offset = np.round((P['origin'][0]-Origin[0])/3).astype(int)
        y_offset = np.round((P['origin'][1]-Origin[1])/3).astype(int)
        z_offset = np.round((P['origin'][2]-Origin[2])/3).astype(int)
        Vn = np.zeros(total_vol.shape)
        Vn[x_offset:V.shape[0]+x_offset,y_offset:V.shape[1]+y_offset,z_offset:V.shape[2]+z_offset] = density_vol
        total_vol = total_vol+Vn
    
        ## Normalize by max color
        Vn2 = np.zeros(total_vol.shape)
        Vn2[x_offset:V.shape[0]+x_offset,y_offset:V.shape[1]+y_offset,z_offset:V.shape[2]+z_offset] = P['vol']/maxcol
        total_D = total_D+Vn2

    if 1:
        imagedata = volume_to_imagedata_respace(np.float32(total_vol>=n_overlap), origin=Origin, spacing = [3,3,3])
    else:
        D = total_D/np.sum(total_D)
        nz_D = D[np.nonzero(D)]
        sorted_data = np.sort(np.ravel(nz_D))[::-1]
        percentile_D = sorted_data[np.where(np.cumsum(sorted_data)>float(d))[0][0]] 
        imagedata = volume_to_imagedata_respace(np.float32(D>=percentile_D), origin=Origin, spacing = [3,3,3])
 
    surface = vtk.vtkMarchingCubes()
    surface.SetInputData(imagedata)
    surface.ComputeNormalsOn();
    surface.SetValue(0, 0.5);
    surface.Update()
    mesh = surface.GetOutput()
    mesh_combined = mesh
    #  color = np.array([0.1,0.1,0.1])
    if 0: #append overlap surface
        A = actor_mesh(mesh, color = (0.,0.,0.4), wireframe=False,opacity = 0.2)
        volume_actor_list.append(A)

#%%
def get_markers_within_shape(mesh,all_markers):
    thickness = 5
    cut_plane_normal = (1,.0,.0)
    markers_in_struct = list()
    for OZ in range(-200,200,10):
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
            if T:# or B:
                markers_in_struct.append([x,y,z])
    return markers_in_struct

#all_markers_list = load_atlas_marker_actors(stack = 'RV4', warp_setting = 24,markertype=str('All'))
#A = get_markers_within_shape(mesh_combined,all_markers_list)
#print A
    
#%%# Threshold points and fit to GMM
if 1:
    #include_stacks = {'RV4'}
    
    all_markers = {}
    all_markers_list = ()
    subset_points = dict()
    for stack in include_stacks:            
        all_markers[stack] = load_atlas_marker_actors(stack = stack, warp_setting = 24,markertype=str('All'))
        if len(all_markers_list) ==0:
            all_markers_list = all_markers[stack]
        else:
            all_markers_list = np.vstack((all_markers_list,all_markers[stack]))
            
        
        points_in_surface = get_markers_within_shape(mesh_combined,all_markers[stack])
        points_in_surface = np.asarray(points_in_surface)
        
        subset_points[stack] = points_in_surface
        
        
            
        for pt in points_in_surface:
            volume_actor_list.append(actor_sphere(pt,radius = 1,color = stacks_to_color[stack]))
    
        gmm = mixture.GaussianMixture(n_components=2, covariance_type='diag',max_iter = 1000)
        gmm.fit(points_in_surface)
        
        print stack
        print gmm.weights_
        
        GMM_center = gmm.means_[0]        
        GMM_center = gmm.means_[1]
        #volume_actor_list.append(actor_sphere(GMM_center,radius = 2,color = stacks_to_color[stack]))
    
    weights = ()
    subset_markers_list = ()
    for stack in include_stacks:
        Npts =  subset_points[stack].shape[0]
        W = np.ones((Npts,1))*(60./Npts)
        if len(subset_markers_list) ==0:
            subset_markers_list = subset_points[stack]
            weights = W
        else:
            subset_markers_list = np.vstack((subset_markers_list,subset_points[stack]))
            weights = np.vstack((weights,W))
        
    points_in_surface = get_markers_within_shape(mesh_combined,all_markers_list)
    points_in_surface = np.asarray(points_in_surface)
    
    subset_points['All'] = points_in_surface
        
    if 1:## Pomegranate
        model = GeneralMixtureModel.from_samples(MultivariateGaussianDistribution, n_components=2, X=subset_markers_list, weights = weights.flatten())
            
        GMM_center = model.distributions[0].parameters[0]
        volume_actor_list.append(actor_sphere(GMM_center,radius = 4,color = (0,0,0)))
        
        GMM_center = model.distributions[1].parameters[0]
        volume_actor_list.append(actor_sphere(GMM_center,radius = 4,color = (0,0,0)))
        
        
        points = subset_markers_list
        mu = np.vstack((model.distributions[0].parameters[0],model.distributions[1].parameters[0]))
        sigma1 = np.linalg.eig(model.distributions[0].parameters[1])[0]
        sigma2 = np.linalg.eig(model.distributions[1].parameters[1])[0]

        sigma1 = np.sqrt(np.diagonal(model.distributions[0].parameters[1]))
        sigma2 = np.sqrt(np.diagonal(model.distributions[1].parameters[1]))
    
        stdev = np.sqrt(np.vstack((sigma1,sigma2)))
        w = model.weights
        
        A,_ = actor_ellipse_vtk(position = (mu[0][0],mu[0][1],mu[0][2]),radius_mat = (3*sigma1[0],3*sigma1[1],3*sigma1[2]),color=(1., 0., 1.),opacity = 0.3)
        volume_actor_list.append(A)
        #A = actor_ellipse_vtk(position = (mu[0][0],mu[0][1],mu[0][2]),radius_mat = (sigma1[0],sigma1[1],sigma1[2]),color=(1., 0., 1.),opacity = 0.5)
        #volume_actor_list.append(A)
        
        A,_ = actor_ellipse_vtk(position = (mu[1][0],mu[1][1],mu[1][2]),radius_mat = (3*sigma2[0],3*sigma2[1],3*sigma2[2]),color=(1., 0., 1.),opacity = 0.3)
        volume_actor_list.append(A)
        #A = actor_ellipse_vtk(position = (mu[1][0],mu[1][1],mu[1][2]),radius_mat = (sigma2[0],sigma2[1],sigma2[2]),color=(1., 0., 1.),opacity = 0.5)
        #volume_actor_list.append(A)
                
    if 0: ##sklearn
        gmm = mixture.GaussianMixture(n_components=2, covariance_type='diag',max_iter = 1000)
        gmm.fit(points_in_surface)
            
        GMM_center = gmm.means_[0]
        volume_actor_list.append(actor_sphere(GMM_center,radius = 4,color = (1,0,1)))
            
        GMM_center = gmm.means_[1]
        volume_actor_list.append(actor_sphere(GMM_center,radius = 4,color = (1,0,1)))
        
    #Visualization
    #visualize_3d_gmm(points_in_surface, gmm.weights_, gmm.means_.T, np.sqrt(gmm.covariances_).T)
#%%
    
#%%
def get_AIC(realValues,estimatedValues,K):
    n = np.float(len(realValues))
    res = subset_points[stack] - estimatedValues
    RSS = np.sum(res**2)
    K=1.
    AIC=n*log(RSS/n)+2*K
    K = 1.
    return AIC
    #aicc_adj = (2*K*K + 2*K)/(n-K-1)
    
    
    ## TEST THE GMM

gmm1 = mixture.GaussianMixture(n_components=1, covariance_type='diag',max_iter = 100,n_init=1000)
gmm2 = mixture.GaussianMixture(n_components=2, covariance_type='diag',max_iter = 100,n_init=1000)
gmm3 = mixture.GaussianMixture(n_components=3, covariance_type='diag',max_iter = 100,n_init=1000)
#gmm4 = mixture.GaussianMixture(n_components=4, covariance_type='diag',max_iter = 100,n_init=1000)





for stack in ['All']:
    gmm1.fit(subset_points[stack])
    gmm2.fit(subset_points[stack])
    gmm3.fit(subset_points[stack])
#    gmm4.fit(subset_points[stack])
    
    bic = list()
    bic.append(gmm1.bic(subset_points[stack]))
    bic.append(gmm2.bic(subset_points[stack]))
    bic.append(gmm3.bic(subset_points[stack]))
#    bic.append(gmm4.bic(subset_points[stack]))
    
    aic = list()

    
    aic.append(gmm1.aic(subset_points[stack]))
    aic.append(gmm2.aic(subset_points[stack]))
    aic.append(gmm3.aic(subset_points[stack]))
    
    
#    aic.append(gmm4.aic(subset_points[stack]))
    
    plt.plot([1,2,3],bic,color = [1,0,1])
    plt.plot([1,2,3],aic,color = [0,1,1])
    
#fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/BIC_AIC.pdf'
#plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

#%%
## Load stuff for atlas display
if 1:
    structs_to_plot = ['7N','7n','IO','LRT','Amb','5N','SpVI','SpVO','SpVC','Brainstem']#,'PreBotC'],#,'SCInG','SCSuG','IC','ZI','RN','fr','scp','SNR']
    structs_to_plot_sided = list()
    for n in structs_to_plot:
        if n in singular_structures:
            structs_to_plot_sided.append(n)
        else:
            structs_to_plot_sided.append(convert_to_left_name(n))
            structs_to_plot_sided.append(convert_to_right_name(n))
        
    for name_s in structs_to_plot_sided:
        fp = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=name_s)
        fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=name_s)
        mesh_rel2canon = load_mesh_stl(fp,return_polydata_only=True)
         
        
        color = np.array(name_unsided_to_color[convert_to_original_name(name_s)])/255.
        color = np.array([0.1,0.1,0.1])
        A = actor_mesh(mesh_rel2canon, color = color, wireframe=False,opacity = 0.1)
        volume_actor_list.append(A)
        #structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas[name_s][1])
        #structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas[name_s][2])
        #structure_mesh_actors_rel2canon.append(slice_actor_dict_atlas[name_s][3])

volume_actor_list.append(make_scaleBar_actor((5,0,400),(1,0,0)))
volume_actor_list.append(make_scaleBar_actor((5,0,400),(0,1,0)))
volume_actor_list.append(make_scaleBar_actor((5,0,400),(0,0,1)))

launch_vtk(volume_actor_list, interactive = False, snapshot_fn='/home/asya/Documents/Yuncong_code/overlayB2.png', snapshot_magnification=1)
#launch_vtk(volume_actor_list,init_angle = 'sagittal', interactive = False, snapshot_fn='/home/asya/Documents/Yuncong_code/overlayB1.png', snapshot_magnification=1)    

#%%###
print array([ -82.0,   23.7,  115.3])*15/1000 + (-9.5,7.4, 0.)
print array([-130.2,   18.3,  101.3])*15/1000 + (-9.5,7.4, 0.)
