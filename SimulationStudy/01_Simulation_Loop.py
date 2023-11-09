#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 15:28:19 2022

@author: david
"""

# This code imports the function for simulating a multi-omics data set
# based on the specified number of simulation iterations
# it simulates multi-omics datasets using the SimMultiOmics function from the imported Multi_Omics_Simulator module
# for each simulation instance a folder is created and the simulated data is saved in csv files
# a single json file is saved - it contains metadat for the simulation configuration

import sys
import os
import numpy as np
import json

# set the random seed

seed = 1234567
np.random.seed(seed)

# set base directory as the current working directory which should be the script's parent folder
BasePath = os.getcwd()

#  import the module containing the simulation function - expected to be in same folder as this script
sys.path.append(os.path.join(BasePath))
import Multi_Omics_Simulator as sm

#####
###### Set the output folders
#####

# Set the name of the parent folder that will contain the simulated data for this configuration
Config = "SimConfig1"
ConfigDir = os.path.join(BasePath,Config)

# folder name to hold data for each simulation as well as name of output meadata file - THESE DONT NEED TO BE CHANGED   
SimsDir = os.path.join(ConfigDir,"SimData")   
sim_meta_out = os.path.join(ConfigDir,"sim_meta.json")


#####
##### set the simulation parameters
#####

Sims = 3

CLrn = 20 
CLrnNPool = [10,20,30] 
CTrg = 2 
CTrgN = 5
K = 20 
ZMuBase = 5
ZSdRange = [1,1]
ZMuShift = [-2,2] 
ZShiftP = 1/4

D = [500, 500, 500] 
YDists = ["Poisson", "Normal", "Bernoulli"]


WThetaRange = [0.15, 0.25]
WMuPois = 5 
WSdRangePois = [0.5, 1.5]
WMuNorm = 0
WSdRangeNorm = [0.5, 1.5]
WMuBeta = 0 
WSdRangeBeta = [0.25, 0.75]
WMuBern = 0 
WSdRangeBern = [0.10, 0.20]
YSdRangeNorm = [0.25, 0.75] 
YaRangeBeta = [50, 150]

#########
#### simulate datasets and save to computer
#########

for s in range(Sims):  
        
    # simulate the data using the imported function
    SimData = sm.SimMultiOmics(
        CLrn = CLrn, CLrnNPool = CLrnNPool, CTrg = CTrg, CTrgN = CTrgN, 
        K = K, ZMuBase = ZMuBase, ZSdRange = ZSdRange, ZMuShift = ZMuShift, ZShiftP = ZShiftP, 
        D = D, YDists = YDists, 
        WThetaRange = WThetaRange, 
        WMuPois = WMuPois, WSdRangePois = WSdRangePois, 
        WMuNorm = WMuNorm, WSdRangeNorm = WSdRangeNorm, 
        WMuBeta = WMuBeta, WSdRangeBeta = WSdRangeBeta, 
        WMuBern = WMuBern, WSdRangeBern = WSdRangeBern,
        YSdRangeNorm = YSdRangeNorm, YaRangeBeta = YaRangeBeta
        )
    
    # create directories to store results in
    
    # these root directories only created during first sim
    
    if not os.path.exists(ConfigDir):
        os.mkdir(ConfigDir)
        
    if not os.path.exists(SimsDir):
        os.mkdir(SimsDir)
        
    # a results folder for each sim
    
    sim_out_dir = os.path.join(SimsDir,"sim" + str(s))

    if not os.path.exists(sim_out_dir):
        os.mkdir(sim_out_dir)
    
    # save simulated data to file
        
    np.savetxt(os.path.join(sim_out_dir, 'CIds.csv'), SimData['CIds'], delimiter=',')
    np.savetxt(os.path.join(sim_out_dir, 'IfLrn.csv'), SimData['IfLrn'], delimiter=',')
    np.savetxt(os.path.join(sim_out_dir, 'Nc.csv'), SimData['Nc'], delimiter=',')
    np.savetxt(os.path.join(sim_out_dir, 'ZMuC.csv'), SimData['ZMuC'], delimiter=',')
    
    for m in range(len(D)):
        np.savetxt(os.path.join(sim_out_dir, 'W' + str(m) + '.csv'), SimData['W'][m], delimiter=',')
        np.savetxt(os.path.join(sim_out_dir, 'Y' + str(m) + '.csv'), SimData['Y'][m], delimiter=',')
        
    np.savetxt(os.path.join(sim_out_dir, 'Z.csv'), SimData['Z'], delimiter=',')

#### save the meta data

sim_meta = {
    'Sims': Sims,
    'CLrn': CLrn, 'CLrnNPool': CLrnNPool, 'CTrg': CTrg, 'CTrgN': CTrgN, 
    'K': K, 'ZMuBase': ZMuBase, 'ZSdRange': ZSdRange, 'ZMuShift': ZMuShift, 'ZShiftP': ZShiftP, 
    'D': D, 'YDists': YDists, 
    'WThetaRange': WThetaRange, 
    'WMuPois': WMuPois, 'WSdRangePois': WSdRangePois, 
    'WMuNorm': WMuNorm, 'WSdRangeNorm': WSdRangeNorm, 
    'WMuBeta': WMuBeta, 'WSdRangeBeta': WSdRangeBeta,
    'WMuBern': WMuBern, 'WSdRangeBern': WSdRangeBern,
    'YSdRangeNorm': YSdRangeNorm, 'YaRangeBeta': YaRangeBeta,
    'seed': seed
    }


with open(sim_meta_out, "w") as outfile:
    json.dump(sim_meta,outfile)
    
    
print("simulation loop finished, data saved in " + ConfigDir)    
