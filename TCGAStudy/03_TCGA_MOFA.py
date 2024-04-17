#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 11:02:59 2023

@author: davidh
"""

#############
## FACTORIZE TCGA MULTI-OMICS LEARNING, REFERENCE, OR TARGET DATASETS WITH MOFA
#############

## This script is divided into 2 main sections
## Sections are run independently, depending on the desired task

## 01 LEARNING OR REFERENCE DATASET FACTORIZATION WITH MOFA
## 02 TARGET DATASET FACTORIZATION WITH MOFA

## NOTE: 
## What we refer to in the MOTL paper as a REFERENCE datasets we originally named FULL TARGET datasets
## Hence the names of objects, created by our code, that correspond to REFERENCE datasets typically include 'trg_XXX_full'
## What we refer to in the MOTL paper as a TARGET dataset we originally named TARGET SUSBETS
## Hence the names of objects, created by our code, that correspond to TARGET datasets typically include 'trg_XXX_ss'

## Load libraries
from time import time
import os
import json

## LOAD FUNCTIONS
from MOFA_functions import *

#############################################################################
## 01 LEARNING OR REFERENCE DATASET FACTORIZATION WITH MOFA
#############################################################################

## This section is for factorizaing a TCGA multi-omics set built from any combination of mRNA, miRNA, DNAme and SNV data
## Either a LEARNIN dataset or a REFERENCE (FULL TARGET) dataset can be factorized using this section

## PARAMETERS

## input
Prjct = "Trg_LAML_SKCM_Full" # either 'Lrn' for LEARNING dataset or 'Trg_Project1_Project2_Full' for REFERENCE dataset
TopD = '5000D'  ## number of features maintained during variance filtering XD
Prior_K = 100 # how many factors to start with - 100 for LEARNING, 60 for single-project REFERENCE, 100 for multi-project REFERENCE

# omics = ['mRNA','miRNA', 'DNAme', 'SNV']
# likelihoods = ['gaussian', 'gaussian', 'gaussian', 'bernoulli']

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
## threshold of 0.001 for LEARNING and either 0.01 for REFERENCE
training_options = {
    "TrainingIter": 10000, ## max training iterations
    "freqELBOChk": 5, ## how often to check the elbo, this is the default for R so using this
    "drop_factor_threshold": 0.01,
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
saveResultsModel(ent, OutputDir, training_options, data_options, G, likelihoods, script_start_time, script_end_time)


#############################################################################
## 02 TARGET DATASET FACTORIZATION WITH MOFA
#############################################################################

## PARAMETERS
## input
Prjct = "Trg_LAML_SKCM_SS5" # format = 'Trg_Project1_Project2_SSX' where X is the number of samples per project in a TARGET dataset
TopD = '5000D'  ## number of features maintained during variance filtering XD
Prior_K = 10 # how many factors to start with - total number of samples for TARGET datasets

# omics = ['mRNA','miRNA', 'DNAme', 'SNV']
# likelihoods = ['gaussian', 'gaussian', 'gaussian', 'bernoulli']

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
## threshold of 0.01 for TARGET datasets
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

## FACTORIZE EACH TARGET DATASET
SS_count = expdat_meta_base['SS_count'] #number of TARGET datasets (subsets of the REFERENCE dataset)

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
    saveResultsModel(ent, OutputDir, training_options, data_options, G, likelihoods, script_start_time, script_end_time)
