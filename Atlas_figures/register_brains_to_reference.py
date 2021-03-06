#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 14:22:28 2018

@author: asya
"""

import sys
import os
import time
import matplotlib.pyplot as plt
import numpy as np
sys.path.append('/home/asya/Documents/Yuncong_code/utilities')
from utilities2015 import *
from registration_utilities import *
from annotation_utilities import *
from metadata import *
from data_manager import *
from aligner import *
from registration_utilities_v0 import save_alignment_results


from registration_utilities_2 import *

GR_stacks = {'GR23_L', 'GR26_L', 'GR30_L', 'GR35_L', 'GR40_L',
'GR24_L', 'GR27_L', 'GR31_L', 'GR37_L', 'GR41_L',
'GR25_L', 'GR29_L', 'GR33_L', 'GR38_L', 'GR42_L'}

## Set the stack moving, and fixed
stack_fixed = 'DTA04_L'
stacks_moving = {'Rat_brainstem_atlas'} # RV19,RV4, RV14, RV13

warp_setting = 24# 24 is affine transformation

def run_aligner(stack_moving,stack_fixed,warp_setting, showplots = 'off', save = 'on'):
    alignment_name_dict=dict(stack_m=stack_moving, 
                                  stack_f=stack_fixed,
                                  warp_setting=warp_setting,
                              vol_type_f='annotationAsScore',
                              vol_type_m='annotationAsScore')
    
    
    aligner_parameters = generate_aligner_parameters(warp_setting=warp_setting, 
                                                     alignment_name_dict=alignment_name_dict)
    
    volume_fixed = aligner_parameters['volume_fixed']
    volume_moving = aligner_parameters['volume_moving']
    
    
    aligner = Aligner4(volume_fixed, volume_moving, labelIndexMap_m2f=aligner_parameters['label_mapping_m2f'])
    
    aligner.set_centroid(centroid_m='structure_centroid', centroid_f='structure_centroid', 
                         indices_m=[aligner_parameters['structure_to_label_moving']['7N_L']])
    
    aligner.set_label_weights(label_weights=aligner_parameters['label_weights_m'])

    
    gradients_f = compute_gradient(volume_fixed, smooth_first=True)
    aligner.load_gradient(gradients=gradients_f) # 120s-170 = 2 mins; all 28, 220s
    
    trial_num = 1
    
    T_all_trials = []
    scores_all_trials = []
    traj_all_trials = []
    
    for _ in range(trial_num):
    
        try:
            T, scores = aligner.optimize(tf_type=aligner_parameters['transform_type'], 
                                         max_iter_num=100,
                                         history_len=20, 
                                         terminate_thresh_rot=0.0001,
                                         terminate_thresh_trans=0.01,
                                         grad_computation_sample_number=aligner_parameters['grad_computation_sample_number'],
                                         lr1=1, lr2=.1,
    #                                     init_T=grid_search_T, 
                                          affine_scaling_limits=(.4, 0.6)
                                        )
            T_all_trials.append(T)
            scores_all_trials.append(scores)
            traj_all_trials.append(aligner.Ts)
            
        except Exception as e:
            sys.stderr.write('%s\n' % e)
    
    ## select the best trial here
    Ts = np.array(aligner.Ts)
    best_trial = np.argsort([np.max(scores) for scores in scores_all_trials])[-1]
    
    # best_trial = 1
    T = T_all_trials[best_trial]
    scores = scores_all_trials[best_trial]
    print 'Best trial:', best_trial
    print max(scores), scores[-1]
    
    crop_origin_m = aligner_parameters['volume_moving_bbox'][[0,2,4]]
    print crop_origin_m
    
    crop_origin_f = aligner_parameters['volume_fixed_bbox'][[0,2,4]]
    print crop_origin_f
    
    if showplots == 'on':

        plt.plot(Ts[:, [0,1,2,4,5,6,8,9,10]]);
        plt.title('rotational params');
        plt.xlabel('Iteration');
        plt.show();
        
        plt.plot(Ts[:, [3,7,11]]);
        plt.title('translation params');
        plt.xlabel('Iteration');
        plt.show();
    
        print T.reshape((3,4))
        plt.figure();
        plt.plot(scores);
        plt.show();
    
    ## save the best trial here
    if save == 'on':
        save_alignment_results(T_all_trials[best_trial], 
                       aligner.centroid_m, aligner.centroid_f, 
                       crop_origin_m, crop_origin_f,
                       score_traj=scores_all_trials[best_trial],
                       parameter_traj=traj_all_trials[best_trial],
                      alignment_name_dict=alignment_name_dict)
    
 #   return aligner,T_all_trials, scores_all_trials,traj_all_trials

for stack_moving in stacks_moving:
    run_aligner(stack_moving,stack_fixed,warp_setting, showplots = 'off', save = 'on')
