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


do_volume_display = 1

thickness_mm = 0.15
##Saggital
cut_plane_origin_mm = (0.,0.,1.75)
cut_plane_normal = (0,0,1)

##Coronal
#cut_plane_origin_mm = (-200,0,0)
#cut_plane_normal = (1,0,0)


cut_plane_origin = np.asarray(cut_plane_origin_mm)*1000./15
thickness = thickness_mm*1000/15

stack ='RV14'
include_stacks = {stack}#'RV14'} #
atlas_name = 'Rat_brainstem_atlas'
spacefill = True

contour_densities = {0.01,0.05,0.1,0.25}

#%% 
## Load the stuff
if 1:
    slice_coords_dict_atlas, temp = load_atlas_slice(atlas_name = atlas_name,cut_plane_origin = cut_plane_origin,cut_plane_normal = cut_plane_normal, spacefill = True)
    fp = DataManager.get_density_pvol_filepath(stack,0.13333333333333)
    P = pickle.load(open(fp,'rb'))

#%%
fig = pl.figure(figsize=(20,10))
ax = fig.add_subplot(111)

plot_atlas_slice(ax,slice_coords_dict_atlas,cut_plane_normal,usecolors = False)
load_plot_density(ax,P,cut_plane_origin,cut_plane_normal)

contours = load_density_contours_v2(P,cut_plane_origin,cut_plane_normal, contour_densities = contour_densities)

plot_density_contours(ax, contours = contours, contour_densities = contour_densities,basecolor = [0,0,1])
#plot_density_contours_from_Volume(ax,stack,cut_plane_origin,cut_plane_normal,contour_densities = contour_densities,step = 5)

plot_markers_slice(ax,stack,cut_plane_normal,cut_plane_origin,col = (0.,0.1,0.))


ax.set_yticks([])
ax.set_xticks([-3,-2,-1,0,1,2,3])
ax.set_aspect(1.0)

matplotlib.pyplot.text(-3, 2, 'Lateral = %s'%cut_plane_origin_mm[2])

## create legend for the densities
basecolor = [0,0,1]
patch_handles = []
for d in contour_densities:
    colorscale = d/max(contour_densities)
    color = np.asarray(basecolor)*colorscale
    patch_handles.append(mpatches.Patch(color=color, label='Density %s'%d))
plt.legend(handles=patch_handles)

fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Fig1d_RV14_markers_and_density.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

