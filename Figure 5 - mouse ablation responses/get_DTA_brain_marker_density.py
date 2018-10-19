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
from registration_utilities_2 import *
#from conversion import *
from vis3d_utilities import *
from neurolucida_to_volume_utilities import *

from vis3d_utilities_stacy import *
from volume_display_utilities import *
import cPickle as pickle
import vtk
import scipy.ndimage as snd
import scipy.stats as st
from shapely import geometry
import shapely.affinity
import pandas as pd
from scipy.stats import ttest_ind
import seaborn as sns

include_stacks = {'DTA03_R', 
                  'DTA06_R', 
                  'DTA05_R', 
                  'DTA10_R',
                  'DTA04_R',
                  'DTA01_L',
                  'DTA05_L',#
                  'DTA10_L', #
                  'DTA02_L', 
                  'DTA04_L', 
                  'DTA07_L', 
                  'DTA03_L', 
                  'DTA06_L',  
                  'DTA09_L',
                  'DTA11_L'
                  }
stack_fixed = 'DTA04_L'


stack_to_type = {'DTA03_R':'R', 
                  'DTA06_R':'R', 
                  'DTA05_R':'R', 
                  'DTA10_R':'R',
                  'DTA04_R':'R',
                  'DTA01_L':'Les',
                  'DTA02_L':'Les',
                  'DTA03_L':'No',
                  'DTA04_L':'No',
                  'DTA05_L':'Les',
                  'DTA07_L':'Les', 
                  'DTA06_L':'Les',                    
                  'DTA09_L':'No',
                  'DTA10_L':'Les',
                  'DTA11_L':'No'}

stack_to_col = {'DTA03_R':'b', 
                  'DTA06_R':'b', 
                  'DTA05_R':'b', 
                  'DTA10_R':'b',
                  'DTA04_R':'b',
                  'DTA01_L':'m',
                  'DTA02_L':'m',
                  'DTA03_L':'g',
                  'DTA04_L':'g',
                  'DTA05_L':'m',
                  'DTA07_L':'m', 
                  'DTA06_L':'m',                    
                  'DTA09_L':'g',
                  'DTA10_L':'m',
                  'DTA11_L':'g'}



structs_to_colors_slice = {'7n_R':(0.7,0.7,0), '7n_L':(0.7,0.7,0),'7N_R':(0.5,0,0),'7N_L':(0.5,0,0),'5N_R':(0,0,0.5),'5N_L':(0,0,0.5),'LRt_R':(0,0,0.8),'LRt_L':(0,0,0.8),'Amb_L':(1,.5,.5),'Amb_R':(1,.5,.5),'Brainstem':(0.7,0.7,0.7),'Lesion':(1.,0.,0.)}

aligned_markers, raw_markers, aligned_lesion_markers, lesion_markers = load_aligned_marker_actors(include_stacks=include_stacks, stack_fixed = stack_fixed,warp_setting = 24, markertype = str('All'))

#%%
if 0:
    stack1 = 'DTA05_L'
    plt.plot(-aligned_markers[stack1][:,0],-aligned_markers[stack1][:,1],'b.')
    stack2 = 'DTA07_L'
    plt.plot(-aligned_markers[stack2][:,0],-aligned_markers[stack2][:,1],'c.')
    plt.plot(-aligned_lesion_markers[stack2][:,0],-aligned_lesion_markers[stack2][:,1],'m.')


GMM = pickle.load(open('/home/asya/Documents/Yuncong_code/GMM_fit_results.pickle','rb'))
tf_parameter_dict = load_alignment_parameters_v2(stack_f='DTA04_L', stack_m='Rat_brainstem_atlas', warp_setting=24, 
                                             vol_type_f='annotationAsScore', vol_type_m='annotationAsScore', downscale=15)
cf = np.array(tf_parameter_dict['centroid_f'])
cm = np.array(tf_parameter_dict['centroid_m'])
of = np.array(tf_parameter_dict['crop_origin_f'])
om = np.array(tf_parameter_dict['crop_origin_m'])
params = np.array(tf_parameter_dict['params'])
    
t= np.asarray( [1,1,1])
R =  np.vstack(([.55,0,0],[0,.55,0],[0,0,.55]))
cf = cf-[-17,15,5]

mu_new0 = np.dot(R, (GMM['mu'][0] - om - cm).T).T + t + of + cf
mu_new1 = np.dot(R, (GMM['mu'][1] - om - cm).T).T + t + of + cf

s0 = np.dot(R, ((GMM['mu'][0]+GMM['sigma0']) - om - cm).T).T + t + of + cf - mu_new0
s1 = np.dot(R, ((GMM['mu'][1]+GMM['sigma1']) - om - cm).T).T + t + of + cf - mu_new1

ellipse = ((mu_new0[0], mu_new0[1]),(2*s0[0], 3*s0[1]),0)

# Let create a circle of radius 1 around center point:
circ = shapely.geometry.Point(ellipse[0]).buffer(1)

# Let create the ellipse along x and y:
ell  = shapely.affinity.scale(circ, int(ellipse[1][0]), int(ellipse[1][1]))


#%%
##
## Density behind facial
#20:40, 60:80
A1 = -1137
A2 = -1150
B1 = 150
B2 = 163

poly1 = ell#geometry.Polygon([[A1,B1],[A2,B1],[A2,B2],[A1,B2]])

#poly1a = geometry.Polygon([[A1+5,B1-5],[A2-10,B1-5],[A2-10,B2+5],[A1+5,B2+5]])

poly2 = geometry.Polygon([[-1110,100],[-1160,130],[-1170,200],[-1110,200]])
poly2 = poly2-poly1

all_counts = pd.DataFrame()
for stack in include_stacks:
    count1 = 0
    count2 = 0
    
    countLesion = 0
    M = aligned_markers[stack][:,[0,1]]
    for x,y in M:
        p = Point(x,y)
        if p.within(poly1):
            count1+=1
        if p.within(poly2):
            count2+=1
    
    ML = aligned_lesion_markers[stack]
    if len(ML)>0:
        ML = ML[:,[0,1]]
        for x,y in ML:
            p = Point(x,y)
            if p.within(poly1):
                countLesion+=1
            
    all_counts = all_counts.append(pd.DataFrame([[stack_to_type[stack],count1,count2,len(M),countLesion]],columns = ['stack','1','2','tot','Lesion']))

all_counts['ratio'] = (all_counts['1'])/(all_counts['2'])#+all_counts['Lesion'])

fig = plt.figure()
sns.stripplot(x = 'stack', y = 'ratio', data = all_counts)
#sns.boxplot(x = 'stack', y = 'ratio', data = all_counts)


MN1 = np.mean(all_counts[all_counts['stack']=='Les'].ratio)
SERR1 = (np.std(all_counts[all_counts['stack']=='Les'].ratio))/sqrt(sum(all_counts['stack']=='Les'))
plt.plot(1,MN1,'ro')
plt.plot([1,1],[MN1+SERR1,MN1-SERR1],'r-')

MN1 = np.mean(all_counts[all_counts['stack']=='No'].ratio)
SERR1 = (np.std(all_counts[all_counts['stack']=='No'].ratio))/sqrt(sum(all_counts['stack']=='No'))
plt.plot(2,MN1,'co')
plt.plot([2,2],[MN1+SERR1,MN1-SERR1],'c-')


MN3 = np.mean(all_counts[all_counts['stack']=='R'].ratio)
SERR3 = (np.std(all_counts[all_counts['stack']=='R'].ratio))/sqrt(sum(all_counts['stack']=='R'))
plt.plot(0,MN3,'bo')
plt.plot([0,0],[MN3+SERR3,MN3-SERR3],'b-')


print ttest_ind(all_counts[all_counts['stack']=='Les'].ratio, all_counts[all_counts['stack']=='R'].ratio, equal_var=False)
print ttest_ind(all_counts[all_counts['stack']=='Les'].ratio, all_counts[all_counts['stack']=='No'].ratio, equal_var=False)
print ttest_ind(all_counts[all_counts['stack']=='No'].ratio, all_counts[all_counts['stack']=='R'].ratio, equal_var=False)


fig_title = '/home/asya/Documents/data/odor_lesion_DTA/cell_density_averages.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)

#%%##




contours, markers = get_stacy_contours('DTA04_L', downscale = 15)
fig2 = plt.figure()
ax = plt.subplot(111)


x, y = ell.exterior.xy
ax.plot(x, y, color=(0,1,1), alpha=1,
            linewidth=1, zorder=1)

x, y = poly1.exterior.xy
ax.plot(x, y, color=(1,0,1), alpha=1,
            linewidth=1, zorder=1)
x, y = poly2.exterior.xy
ax.plot(x, y, color=(1,0,1), alpha=1,
            linewidth=1, zorder=1)


for c in contours['7N_L']:
    if c[0][2] == -45.:
        x = c[:,0]
        y = c[:,1]
        ax.plot(x, y, color=(1,0,0), alpha=1,
            linewidth=1, zorder=1)
for c in contours['IO_L']:
    if c[0][2] == -45.:
        x = c[:,0]
        y = c[:,1]
        ax.plot(x, y, color=(0,1,1), alpha=1,
            linewidth=1, zorder=1)
for c in contours['7n_L']:
    if c[0][2] == -45.:
        x = c[:,0]
        y = c[:,1]
        ax.plot(x, y, color=(1,1,0), alpha=1,
            linewidth=1, zorder=1)
for c in contours['LRT_L']:
    if c[0][2] == -45.:
        x = c[:,0]
        y = c[:,1]
        ax.plot(x, y, color=(0,0,1), alpha=1,
            linewidth=1, zorder=1)
plt.axis('equal')

fig_title = '/home/asya/Documents/data/odor_lesion_DTA/contour_area_show.pdf'
plt.savefig(fig_title, format='pdf',dpi=fig.dpi)