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

import seaborn as sns
do_slice = 1
do_volume_display = 1


##Saggital
cut_plane_origin = (0,0,-90)
cut_plane_normal = (0,0,1)

##Coronal
cut_plane_origin = (-200,0,0)
cut_plane_normal = (1,0,0)

cut_plane_origin = (-290,0,0)
cut_plane_normal = (1,0,0)

#thickness = 1

include_stacks = {'RV4','RV14','RV13','RV19','RV9','RV10'}

stack_to_color = dict()
for i,stack in enumerate(include_stacks):
    stack_to_color[stack] = sns.color_palette()[i]

atlas_name = 'Rat_atlas'
spacefill = True

contour_densities = {0.5} #{0.25,0.5,0.75}

slice_coords_dict_atlas, temp = load_atlas_slice(atlas_name = 'Rat_atlas',cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, spacefill = True)
if 1:
    P = dict()
    for stack in include_stacks:
        fp = DataManager.get_density_pvol_filepath(stack)
        P[stack] = pickle.load(open(fp,'rb'))

fig = pl.figure(figsize=(10,10))
ax = fig.add_subplot(111)
plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)

for stack in include_stacks:
    contours = load_density_contours(P[stack],cut_plane_origin,cut_plane_normal, contour_densities = contour_densities)
    plot_density_contours(ax, contours = contours, contour_densities = contour_densities,basecolor = stack_to_color[stack])

ax.set_yticks([])
ax.set_xticks([-100,0,100])
ax.set_aspect(1.0)

## create legend for the densities
if len(contour_densities)>1:
    basecolor = [1,0,1]
    patch_handles = []
    for d in contour_densities:
        colorscale = (1.2-d/max(contour_densities))
        color = np.asarray(basecolor)*colorscale
        patch_handles.append(mpatches.Patch(color=color, label='Density %s'%d))
    plt.legend(handles=patch_handles)
else:
    patch_handles = []
    for stack in include_stacks:
        patch_handles.append(mpatches.Patch(color=stack_to_color[stack], label='Density %s'%stack))
    plt.legend(handles=patch_handles)


fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Fig1e_retrofacial_1.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)