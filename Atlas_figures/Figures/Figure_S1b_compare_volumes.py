#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 10:06:48 2018

@author: asya
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

import pandas as pd
import seaborn as sns
import matplotlib.lines as mlines

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

atlas_name = 'Rat_brainstem_atlas'

from volume_display_utilities import *

def get_volume_measurements(structures):
    
    instance_volumes = {}
    averaged_volumes = {}
    adjusted_volumes = {}
    
    scale_factor = (15./1000.)**3
    for struct in structures:
        
        ## get the instance volumes
        instance_vol = list()
        for index in range(12):
            fp = DataManager.get_instance_mesh_filepath(atlas_name=atlas_name, structure=struct, index=index)
            if not os.path.isfile(fp):
                continue
            verts,faces = load_mesh_stl(fp)
            V_calc = vtk.vtkMassProperties()
            V_calc.SetInputData(mesh_to_polydata(verts,faces))
            instance_vol.append(V_calc.GetVolume()*scale_factor)
        instance_volumes[struct] = instance_vol
        
        fp_av = DataManager.get_mesh_filepath(stack_m=atlas_name, structure=convert_to_left_name(struct))
        if not os.path.isfile(fp_av):
            print 'not found %s'%struct
            continue
        
        verts,faces = load_mesh_stl(fp_av)
        V_calc = vtk.vtkMassProperties()
        V_calc.SetInputData(mesh_to_polydata(verts,faces))
        averaged_volumes[struct] = V_calc.GetVolume()*scale_factor
        
        '''
        fp_spacefill = DataManager.get_mesh_filepath_spacefill(stack_m=atlas_name, structure=convert_to_left_name(struct))
        if not os.path.isfile(fp_spacefill):
            continue
        
        verts,faces = load_mesh_stl(fp_spacefill)
        V_calc = vtk.vtkMassProperties()
        V_calc.SetInputData(mesh_to_polydata(verts,faces))
        adjusted_volumes[struct] = V_calc.GetVolume()*scale_factor
        '''
    df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in instance_volumes.iteritems()]))
    df_av = (pd.DataFrame(averaged_volumes, index=['averaged']))
    df_adj = df_av
    #df_adj = (pd.DataFrame(adjusted_volumes, index=['adjusted']))
    
    return df, df_av, df_adj

def plot_volume_measurements(df, df_av, df_adj,title = 'Structures',save_eps = False):
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    sns.stripplot(ax=ax,data=df[0:11],color = sns.xkcd_rgb["greyish"],alpha = 0.3,zorder = 0,size = 10)
    sns.pointplot(ax=ax,data=df_av, markers="d",color = 'k', join=False,zorder = 1)
    sns.pointplot(ax=ax,data=df_adj, markers="d",color = sns.xkcd_rgb["rose red"], join=False,zorder = 2)
    ax.set(xlabel='Structure name', ylabel='Volume (mm^3)')
    ax.set(title=title)
    
    handles = []
    handles.append(mlines.Line2D([], [], linestyle='', color=sns.xkcd_rgb["greyish"], marker="o", label='Instance'))
    handles.append(mlines.Line2D([], [], linestyle='', color='k', marker="d", label='Average'))
    handles.append(mlines.Line2D([], [], linestyle='', color=sns.xkcd_rgb["rose red"], marker="d", label='Adjusted average'))
    ax.legend(handles=handles)
        
    if save_eps:
        fig_title = '/home/asya/Documents/Yuncong_code/Figures_code/Atlas_figs/FigS1_'+title+'.pdf'
        plt.savefig(fig_title, format='pdf',dpi=fig.dpi)
 

save_eps = True
#%%
df, df_av, df_adj = get_volume_measurements(alignment_structures)
plot_volume_measurements(df, df_av, df_adj,title = 'Alignment structures',save_eps = True)

#%%
df, df_av, df_adj = get_volume_measurements(['3N', '4N','Dk', 'Pr5','SpVI','SpVO','SPVmu','SCInG','ZI','AP', 'SNR'])
plot_volume_measurements(df, df_av, df_adj,title = 'Other outlined structures',save_eps = save_eps)

#%%
df, df_av, df_adj = get_volume_measurements(['IRt','PCRt','Gi','PnO','PnC','MdV','MdD','SubC','DMSp5','intertrigeminal','GiA','GiV'])
plot_volume_measurements(df, df_av, df_adj,title='Select reticular nuclei',save_eps = save_eps)