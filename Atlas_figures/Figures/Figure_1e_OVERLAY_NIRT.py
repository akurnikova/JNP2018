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


do_volume_display = 1


thickness_mm = 0.1

## Coronal through nIRT
#cut_plane_origin_mm = (-2.2,0,0)
#cut_plane_normal = (1,0,0)

## Saggital through nIRT
cut_plane_origin_mm = (0.,0.,1.55)
cut_plane_normal = (0,0,1)


cut_plane_origin = np.asarray(cut_plane_origin_mm)*1000./15.
thickness = thickness_mm*1000./15.

include_stacks = {'RV4','RV14','RV13','RV19','RV9','RV10'}

stack_to_color = dict()
for i,stack in enumerate(include_stacks):
    stack_to_color[stack] = sns.color_palette()[i]

atlas_name = 'Rat_brainstem_atlas'
spacefill = True

contour_densities = {0.1} #{0.25,0.5,0.75}

#%%
P = dict()
slice_coords_dict_atlas, temp = load_atlas_slice(atlas_name = 'Rat_brainstem_atlas',cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, spacefill = True)
for stack in include_stacks:
    fp = DataManager.get_density_pvol_filepath(stack,0.13333333333333)
    P[stack] = pickle.load(open(fp,'rb'))

## create legend for the densities
patch_handles = []
for stack in include_stacks:
     patch_handles.append(mpatches.Patch(color=stack_to_color[stack], label='Density %s'%stack))
        
### Plot densities as contours
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)

for stack in include_stacks:
    contours = load_density_contours_v2(P[stack],cut_plane_origin,cut_plane_normal, contour_densities = contour_densities)
    plot_density_contours(ax, contours = contours, contour_densities = contour_densities,basecolor = stack_to_color[stack])

ax.set_yticks([])
ax.set_xticks([-3,-2.5,-2.25,-2,-1.75,-1.50,-1.25,-1.00,-.75,-.50,-.25,0.,.25,.50,.75,1.,1.25,1.5,1.75,2.,3.])
ax.set_aspect(1.0)

plt.legend(handles=patch_handles)

fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/Fig1E_NIRT_density_contours_overlay_cor.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)


### Plot densities as points
fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)

for stack in include_stacks:
    plot_markers_slice(ax,stack,cut_plane_normal,cut_plane_origin,col = stack_to_color[stack],thickness = thickness)

ax.set_yticks([])
ax.set_xticks([-3.,-2.,-1.75,-1.50,-1.25,-1.00,-.75,-.50,-.25,0.,.25,.50,.75,1.,1.25,1.5,1.75,2.,3.])
ax.set_aspect(1.0)

plt.legend(handles=patch_handles)

fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/Fig1E_NIRT_point_overlay_cor.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
