#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 11 10:09:38 2018

@author: yuncong
"""
import sys
import os

import numpy as np
import bloscpack as bp

sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from metadata import *
from data_manager import *

from registration_utilities import *
from annotation_utilities import *
from vis3d_utilities import *
from aligner import *

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

from mpl_toolkits.mplot3d import Axes3D

atlas_name = 'Rat_brainstem_atlas'

structures_to_average = ['5N', '7N', '7n', 'Amb', 'LRT','IO']#, 'RN', 'RtTg', 'fr'] #all_known_structures



structures_to_average_sided = list()
for n in structures_to_average:
    if n in singular_structures:
        structures_to_average_sided.append(n)
    else:
        structures_to_average_sided.append(convert_to_left_name(n))
        structures_to_average_sided.append(convert_to_right_name(n))

#%%#Load instance centroids

fp = DataManager.get_instance_centroids_filepath(atlas_name=atlas_name)
instance_centroids_rel2fixed = load_pickle(fp)

canonical_centroid_rel2fixed_mm = (7.4,  -9.5, 0.)
canonical_centroid_mm = canonical_centroid_rel2fixed_mm


instance_centroids_rel2fixed_mm = {}
for name,oldlist in instance_centroids_rel2fixed.iteritems():
    mmlist = list()
    for val in oldlist:
        mmlist.append(val*15./1000.)
    instance_centroids_rel2fixed_mm[name] = mmlist
    

##Average the centroids
nominal_locations_rel2canon_mm, canonical_centroid_rel2fixed_mm, canonical_normal = \
average_location(instance_centroids_rel2fixed_mm)

nominal_locations_rel2fixed = {k: canonical_centroid_rel2fixed_mm + c for k, c in nominal_locations_rel2canon_mm.iteritems()}


### FOR SAVING
canonical_centroid_rel2fixed = (0., 0., 0.)
instance_centroids_rel2canon = {name_s: np.array(instance_centroids) - canonical_centroid_rel2fixed
                                for name_s, instance_centroids in instance_centroids_rel2fixed.iteritems()}

nominal_locations_rel2canon, canonical_centroid_rel2fixed, canonical_normal = \
average_location(instance_centroids_rel2fixed)

save_pickle(nominal_locations_rel2canon, DataManager.get_structure_mean_positions_filepath(atlas_name=atlas_name))



#%%

instance_centroids_subset = {}
for name in structures_to_average_sided:
    instance_centroids_subset[name] = instance_centroids_rel2fixed_mm[name]

cov_mat_allStructures, radii_allStructures, ellipsoid_matrix_allStructures = \
compute_covar_from_instance_centroids(instance_centroids_subset)


#%%



fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111, projection='3d')

for name_s, centroids in instance_centroids_subset.iteritems():
    if name_s[-1] == 'R':
        continue
    color_plot = np.asarray(name_unsided_to_color[convert_to_original_name(name_s)])/255.


    centroids2 = np.array(centroids) + canonical_centroid_mm

    ax.scatter(centroids2[:,0], centroids2[:,1], centroids2[:,2],
               marker='o', s=100, alpha=.1, color = color_plot)

    c = nominal_locations_rel2fixed[name_s] + canonical_centroid_mm

    ax.scatter(c[0], c[1], c[2],
               color=color_plot, marker='*', s=100)

    # Plot uncerntainty ellipsoids
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii_allStructures[name_s][0] * np.outer(np.cos(u), np.sin(v))
    y = radii_allStructures[name_s][1] * np.outer(np.sin(u), np.sin(v))
    z = radii_allStructures[name_s][2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(u)):
        for j in range(len(v)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], ellipsoid_matrix_allStructures[name_s]) + c

#     ax.plot_surface(x, y, z, color='b')
    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=color_plot, alpha=0.2)

# plot canonical centroid
#ax.scatter(canonical_centroid[0], canonical_centroid[1], canonical_centroid[2], color=(0,0,0), marker='^', s=200)


# Plot mid-sagittal plane
if 0:
    R1 = (nominal_locations_rel2canon_mm['LRT_L'])
    R2 = (nominal_locations_rel2canon_mm['5N_L'])
    canonical_midplane_xx, canonical_midplane_yy = np.meshgrid(range(int(R1[0]), int(R2[0]), 100), range(int(R2[1]), int(R1[1]), 100), indexing='xy')
    canonical_midplane_z = -(canonical_normal[0]*(canonical_midplane_xx-canonical_centroid_mm[0]) + \
    canonical_normal[1]*(canonical_midplane_yy-canonical_centroid_mm[1]) + \
    canonical_normal[2]*(-canonical_centroid_mm[2]))/canonical_normal[2]
    
    ax.plot_surface(canonical_midplane_xx, canonical_midplane_yy, canonical_midplane_z, alpha=.2)

# ax.set_axis_off()
#ax.set_xlim3d([xlim[0], xlim[1]]);
#ax.set_ylim3d([ylim[0], ylim[1]]);
#ax.set_zlim3d([zlim[0], zlim[1]]);
ax.view_init(azim = 270, elev = 0)
# Hide y-axis (https://stackoverflow.com/questions/12391271/matplotlib-turn-off-z-axis-only-in-3-d-plot)
ax.w_yaxis.line.set_lw(0.)
ax.set_yticks([])

ax.set_aspect(1.0)
ax.set_title('Covariance ellipsoids for alignment structures')
plt.xlabel('A-P position (mm)', fontsize=14,labelpad = 30)
plt.ylabel('M-L position (mm)', fontsize=14,labelpad = 0)
plt.legend()
plt.show()

fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/FigS1a_covariance_for_positions.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

#%%


fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111, projection='3d')

for name_s, centroids in instance_centroids_subset.iteritems():
    if name_s[-1] == 'R':
        continue
    color_plot = np.asarray(name_unsided_to_color[convert_to_original_name(name_s)])/255.


    centroids2 = np.array(centroids) + canonical_centroid_mm

    ax.scatter(centroids2[:,0], centroids2[:,1], centroids2[:,2],
               marker='o', s=100, alpha=.1, color = color_plot)

    c = nominal_locations_rel2fixed[name_s] + canonical_centroid_mm

    ax.scatter(c[0], c[1], c[2],
               color=color_plot, marker='*', s=100)

    # Plot uncerntainty ellipsoids
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    x = radii_allStructures[name_s][0] * np.outer(np.cos(u), np.sin(v))
    y = radii_allStructures[name_s][1] * np.outer(np.sin(u), np.sin(v))
    z = radii_allStructures[name_s][2] * np.outer(np.ones_like(u), np.cos(v))
    for i in range(len(u)):
        for j in range(len(v)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], ellipsoid_matrix_allStructures[name_s]) + c

#     ax.plot_surface(x, y, z, color='b')
    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=color_plot, alpha=0.2)

# plot canonical centroid
#ax.scatter(canonical_centroid[0], canonical_centroid[1], canonical_centroid[2], color=(0,0,0), marker='^', s=200)


# Plot mid-sagittal plane
if 0:
    R1 = (nominal_locations_rel2fixed['LRT_L'])
    R2 = (nominal_locations_rel2fixed['RN_L'])
    canonical_midplane_xx, canonical_midplane_yy = np.meshgrid(range(int(R1[0]), int(R2[0]), 100), range(int(R2[1]), int(R1[1]), 100), indexing='xy')
    canonical_midplane_z = -(canonical_normal[0]*(canonical_midplane_xx-canonical_centroid_mm[0]) + \
    canonical_normal[1]*(canonical_midplane_yy-canonical_centroid_mm[1]) + \
    canonical_normal[2]*(-canonical_centroid_mm[2]))/canonical_normal[2]
    
    ax.plot_surface(canonical_midplane_xx, canonical_midplane_yy, canonical_midplane_z, alpha=.2)

# ax.set_axis_off()
#ax.set_xlim3d([xlim[0], xlim[1]]);
#ax.set_ylim3d([ylim[0], ylim[1]]);
#ax.set_zlim3d([zlim[0], zlim[1]]);
ax.view_init(azim = 270, elev = 90)
# Hide y-axis (https://stackoverflow.com/questions/12391271/matplotlib-turn-off-z-axis-only-in-3-d-plot)
ax.w_zaxis.line.set_lw(0.)
ax.set_zticks([])

ax.set_aspect(1.0)
ax.set_title('Covariance ellipsoids for alignment structures')
plt.xlabel('A-P position (mm)', fontsize=14,labelpad = 30)
plt.ylabel('DV position (mm)', fontsize=14,labelpad = 0)
plt.legend()
plt.show()


fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/FigS1a_covariance_for_positionsV.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
