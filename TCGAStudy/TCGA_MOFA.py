#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:02:59 2023

@author: davidh
"""

## This is for factorizaing a TCGA multiomics set built from mRNS, miRNS, DNAme and SNV
## Either a full Trg set or a Lrn set is factorizaed with this

## Load libraries
from time import time
import os
import json
import numpy as np
from mofapy2.run.entry_point import entry_point

## LOAD FUNCTIONS
from MOFA_functions import *

## PARAMETERS

## input
Prjct = "Trg_LAML_SKCM_Full"
TopD = '5000D'  ## number of features maintained during variance filtering XD
Prior_K = 100 # how many factors to start with - 100 for Lrn, trying 30/60 for Trg Full, 100 for Trg multi

omics = ['mRNA','miRNA', 'DNAme']
likelihoods = ['gaussian', 'gaussian', 'gaussian']

M = len(omics)
G = 1 ## groups - always set to 1

data_options = {
    "scale_views": False,
    "scale_groups": False,
    "center_groups": True
}

## drop_factor_threshold: threshold for percentage of variance explained for dropping factors during training
## threshold of 0.001 for Lrn and either 0.01, 0.005 or 0.001 for Trg Full
training_options = {
    "TrainingIter": 10000, ## max training iterations
    "freqELBOChk": 5, ## how often to check the elbo, this is the default for R so using this
    "drop_factor_threshold": 0.001,
    "seed": 1234567
}

## DIRECTORIES
## set up paths and directories
BasePath = os.getcwd()
InputDir = os.path.join(BasePath, Prjct + '_' + TopD)
## set output
OutputDir = os.path.join(InputDir,'Fctrzn_' + str(Prior_K) + 'K_' + str(training_options["drop_factor_threshold"])[2:]+'TH')
if not os.path.exists(OutputDir):
    os.mkdir(OutputDir)

script_start_time = time()

## import meta data
expdat_meta_in = os.path.join(InputDir,"expdat_meta.json")
with open(expdat_meta_in) as infile:
    expdat_meta = json.load(infile)

## IMPORT OMICS DATA
data = import_data(G = G, M = M, omics = omics, InputDir = InputDir, metadata = expdat_meta)

## RUN MOFA
ent = runMOFA(data, data_options, likelihoods, omics, Prior_K, training_options)

script_end_time = time()

## SAVE MODEL AND RESULTS
saveResultsModel(ent, OutputDir, training_options, data_options, G, script_start_time, script_end_time)


#################################################################################################################
#################################################################################################################

#### factorize target datasets (subsets of the full/reference datasets)


## PARAMETERS
## input
Prjct = "Trg_LAML_SKCM_SS5"
TopD = '5000D'  ## number of features maintained during variance filtering XD
Prior_K = 10 # how many factors to start with - number of samples for target subset datasets

omics = ['mRNA','miRNA', 'DNAme']
likelihoods = ['gaussian', 'gaussian', 'gaussian']

M = len(omics)
G = 1 ## groups - always set to 1

data_options = {
    "scale_views": False,
    "scale_groups": False,
    "center_groups": True
}

## drop_factor_threshold: threshold for percentage of variance explained for dropping factors during training
## threshold of 0.001 for Lrn and either 0.01, 0.005 or 0.001 for Trg Full
training_options = {
    "TrainingIter": 10000, ## max training iterations
    "freqELBOChk": 5, ## how often to check the elbo, this is the default for R so using this
    "drop_factor_threshold": 0.01,
    "seed": 1234567
}

## DIRECTORIES
## set up paths and directories
BasePath = os.getcwd()
InputDir = os.path.join(BasePath, Prjct + "_" + TopD)

## import meta data
expdat_meta_base_in = os.path.join(InputDir,"expdat_meta_SS.json")
with open(expdat_meta_base_in) as infile:
    expdat_meta_base = json.load(infile)

## FOR EACH SUBSET
SS_count = expdat_meta_base['SS_count']

for ss in range(SS_count):

    script_start_time = time()

    SS = 'SS_' + str(ss + 1)
    print(SS)

    InputDirSS = os.path.join(InputDir, SS)

    # OutputDir = os.path.join(InputDir,'Fctrzn')
    OutputDir = os.path.join(InputDirSS, 'Fctrzn_' + str(Prior_K) + 'K')
    if not os.path.exists(OutputDir):
        os.mkdir(OutputDir)

    ## import meta data
    expdat_meta_in = os.path.join(InputDirSS, "expdat_meta.json")
    with open(expdat_meta_in) as infile:
        expdat_meta = json.load(infile)

    ## IMPORT OMICS DATA
    data = import_data(G=G, M=M, omics=omics, InputDir=InputDirSS, metadata=expdat_meta)

    ## RUN MOFA
    ent = runMOFA(data, data_options, likelihoods, omics, Prior_K, training_options)

    script_end_time = time()

    ## SAVE MODEL AND RESULTS
    saveResultsModel(ent, OutputDir, training_options, data_options, G, script_start_time, script_end_time)
