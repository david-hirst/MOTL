#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 13:48:28 2022

@author: david
"""

## This is for factorizaing the simulated based on outputs of tuning K number of factors  

## There is an option to select whether learning, target or total dataset is facotrized
## additionally there is an option to select whether data is scaled or not

## Load libraries

import os
import json
import numpy as np
# import sys
from mofapy2.run.entry_point import entry_point

## set up paths
# set base directory as the current working directory which should be the script's parent folder
BasePath = os.getcwd()
# Set the name of the parent folder that will contain the simulated data for this configuration
Config = "SimConfig1"
ConfigDir = os.path.join(BasePath,Config)
SimsDir = os.path.join(ConfigDir,"SimData") 


## set options for factorization
YDataType = 'Lrn' # options are Trg or Lrn or All
scale_views = False  #should views be scaled by standard deviation? True or False
center_groups = True #this is the default could try not centering though
seed = 1234567
FVarTH = 0 #features should have variance higher than this
FVarRankTH = 0 #select filters with rank higher than this - the higher the variance the higher the rank
TrainingIter = 50000 # max iterations for training, shouldn't be necessary but just in case

## facotrization output directory
FctrznsDir = os.path.join(ConfigDir,'Fctrzn_' + YDataType + '_scale_' + str(scale_views) + '_center_' + str(center_groups) + '_ActiveK')
if not os.path.exists(FctrznsDir):
    os.mkdir(FctrznsDir)


## import simulation meta data
sim_meta_in = os.path.join(ConfigDir,"sim_meta.json")
with open(sim_meta_in) as infile:
    sim_meta = json.load(infile)
    
M = len(sim_meta['D'])
YDists = sim_meta['YDists']
Sims = sim_meta['Sims']

## import vector of number of active factors per simulated dataset
## if using the combined dataset then assuming that using number of factors 
## is found by tuning the learning set previously

if YDataType == 'All':
    factorsBySim = np.loadtxt(os.path.join(ConfigDir, 'ActiveKLrn.csv'), delimiter=',')
else:
    factorsBySim = np.loadtxt(os.path.join(ConfigDir, 'ActiveK'+YDataType+'.csv'), delimiter=',')
    
factorsBySim = factorsBySim.astype(int)

###### Set factoriztion options

G = 1 # no groups so hard coded as 1 in order to input data into mofa model

likelihoods = ['gaussian' for m in range(M)]
for m in range(M):
    if YDists[m] == 'Poisson':
        likelihoods[m] = 'poisson'
    elif YDists[m] == 'Bernoulli':
        likelihoods[m] = 'bernoulli'

        

#######
#### Loop through simulated data sets to factorize
#######

for s in range(Sims):

    ## location of simulated data  
    sim_in_dir = os.path.join(SimsDir,"sim" + str(s))
      
    # get datset groupings
    IfLrn = np.loadtxt(os.path.join(sim_in_dir, 'IfLrn.csv'), delimiter=',')
       
    # empty lists for data and feature variance indicators
    
    YData= [[None for g in range(G)] for m in range(M)]
    
    FVar= [[None for g in range(G)] for m in range(M)]
    FVarRk= [[None for g in range(G)] for m in range(M)]
    FVarKp= [[None for g in range(G)] for m in range(M)]
    
    for m in range(M):
        for g in range(G):
            
            YData[m][g] = np.loadtxt(os.path.join(sim_in_dir, 'Y' + str(m) + '.csv'), delimiter=',')
            
            FVar[m][g] = np.loadtxt(os.path.join(sim_in_dir, 'FVar' + YDataType +'Y' + str(m) + '.csv'), delimiter=',')
            FVarRk[m][g] = np.loadtxt(os.path.join(sim_in_dir, 'FVar' + YDataType +'RkY' + str(m) + '.csv'), delimiter=',')                    
            FVarKp[m][g] = (FVar[m][g] > FVarTH) & (FVarRk[m][g] > FVarRankTH)
            
            if YDataType == 'Lrn':
                YData[m][g] = YData[m][g][IfLrn==1,:]                         
            elif YDataType == 'Trg':
                YData[m][g] = YData[m][g][IfLrn==0,:]
                
            YData[m][g] = YData[m][g][:,FVarKp[m][g]]
            
            if YDists[m] == "Beta":
                YData[m][g] = np.log2((YData[m][g] + 0.001)/(1-YData[m][g] + 0.001))
                
    # initialse entry point
    ent = entry_point()
    
     ## DEFINE data options
    ent.set_data_options(
    scale_groups = False, 
    scale_views = scale_views,
    center_groups = center_groups
     )
    
    ## input data as list of matrices
    ent.set_data_matrix(YData, 
                        likelihoods = likelihoods)
    
    
    ## how many factors, min of known active factors and number of samples
    if ent.dimensionalities["N"] < factorsBySim[s]:
        KtoUse = ent.dimensionalities["N"]
    else:
        KtoUse = factorsBySim[s]
    
    
    ## set model options
    ent.set_model_options(factors=KtoUse)
    
    ## SET TRAINING OPTIONS
    # if nothing here for dropR2 then factors will not be dropped
    # set the same seed each time for reproducible results
    # otherwise a random seed is initialised by the mofa model
    ent.set_train_options(iter = TrainingIter, 
                          seed = seed)
    
    # Build the model 
    ent.build()
    
    # Run the model
    ent.run()
    
    # save the results
    
    # create directory for the simulated dataset
    
    FctrznOutDir = os.path.join(FctrznsDir,'Fctrzn' + str(s))
    if not os.path.exists(FctrznOutDir):
        os.mkdir(FctrznOutDir)
        

    # Extract the following which the python function doesn't export:
    # Expected W squared - ordered like the datasets extracted with the save model function
    # Expected log Tau 
    # Additionally the expected Z matrix - unordered -  to use for checking factor ordering
    # and the factor order - as per the function used by MOFA based on variance explained
    
    expAll = ent.model.getExpectations(only_first_moments=False)
    
    np.savetxt(os.path.join(FctrznOutDir, 'ZChk_raw' + '.csv'), expAll['Z']['E'], delimiter=',')
    
    r2 = ent.model.calculate_variance_explained()
    order_factors = np.argsort( np.array(r2).sum(axis=(0,1), where= ~np.isnan(np.array(r2))) )[::-1]
    
    np.savetxt(os.path.join(FctrznOutDir, 'order_factors.csv'), order_factors, delimiter=',')
    np.savetxt(os.path.join(FctrznOutDir, 'r2.csv'), r2[0], delimiter=',')
    
    for m in range(M):
        
        np.savetxt(os.path.join(FctrznOutDir, 'WSq' + str(m) + '.csv'), 
                    expAll['W'][m]['E2'][:,order_factors], delimiter=',')
        
        
        if likelihoods[m] != 'poisson': 
            np.savetxt(os.path.join(FctrznOutDir, 'TauLn' + str(m) + '.csv'), 
                    expAll['Tau'][m]['lnE'][0,:], delimiter=',')
    
    # use the premade function to save results as an hdf5 file - needs to happen after the calculate variance explained function above
    
    ModelOutFile = os.path.join(FctrznOutDir,'Model'  + '.hdf5')
    ent.save(ModelOutFile, 
             # save_data=True,
             expectations="all")
        
    ## export record of features kept
    
    for m in range(M):
        for g in range(G):
            np.savetxt(os.path.join(FctrznOutDir, 'FVarKp_m' + str(m) + 'g_' + str(g) +'.csv'), 
                       FVarKp[m][g], delimiter=',')
    


## metadata to use further down the line

FctrMeta = {
    'YDataType': YDataType,
    'scale_views': scale_views, 
    'center_groups': center_groups,
    'seed': seed,
    'G': G,
    'likelihoods': likelihoods,
    'FVarTH': FVarTH,
    'FVarRankTH': FVarRankTH,
    'TrainingIter': TrainingIter
    }

FctrMeta_out = os.path.join(FctrznsDir,"FctrMeta.json")

with open(FctrMeta_out, "w") as outfile:
    json.dump(FctrMeta,outfile)



 
