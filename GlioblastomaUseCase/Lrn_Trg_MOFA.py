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


## FUNCTIONS

def import_data(G, M, omics, InputDir, metadata):
    print("Read omics")
    YData = [[None for g in range(G)] for m in range(M)]

    for m in range(M):
        for g in range(G):
            print(omics[m])
            YData[m][g] = np.genfromtxt(os.path.join(InputDir, omics[m] + '.csv'), delimiter=',', missing_values='').T

    ## SAMPLE AND FEATURE NAMES RESHAPING
    N = [YData[0][p].shape[0] for p in range(len(YData[0]))]
    D = [YData[m][0].shape[1] for m in range(len(YData))]

    smpls = [[metadata['smpls'][n] for n in range(N[g])] for g in range(G)]
    ftrs = [[metadata['ftrs_' + omics[m]][d] for d in range(D[m])] for m in range(M)]

    MaxK = min(min(N), min(D))

    data = {
        "YData": YData,
        "smpls": smpls,
        "ftrs": ftrs,
        "MaxK": MaxK
    }

    return(data)

def runMOFA(data, data_options, likelihoods, omics, Prior_K, training_options):

    # initialise entry point
    ent = entry_point()

    ## set data options - using the defaults which centers but doesn't scale
    ent.set_data_options(scale_views=data_options["scale_views"],
                         scale_groups=data_options["scale_groups"],
                         center_groups=data_options["center_groups"])

    ## load the data into the MOFA object
    ent.set_data_matrix(data["YData"],
                        likelihoods=likelihoods,
                        views_names=omics,
                        samples_names=data["smpls"],
                        features_names=data["ftrs"])

    ## how many factors to start with -  min of prior ceiling estimate from PCA and number of samples or min features
    if data["MaxK"] < Prior_K:
        KStart = data["MaxK"]
    else:
        KStart = Prior_K

    ## set model options

    ent.set_model_options(factors=KStart)

    ## SET TRAINING OPTIONS
    # im using R default for frequency of checking elbo (5 vs 1 for python) and dropping factors (1 for both)

    ent.set_train_options(iter=training_options["TrainingIter"],
                          freqELBO=training_options["freqELBOChk"],
                          dropR2=training_options["drop_factor_threshold"],
                          seed=training_options["seed"])

    # Build the model
    ent.build()

    # Run the model
    ent.run()

    return(ent)

def saveResultsModel(ent, OutputDir, training_options, data_options, G, script_start_time, script_end_time):

    ## MODEL
    ModelOutFile = os.path.join(OutputDir, 'Model' + '.hdf5')
    ent.save(ModelOutFile, save_data=False, expectations="all")

    ## METADATA
    FctrMeta = {
        "seed": training_options["seed"],
        "TrainingIter": training_options["TrainingIter"],
        "freqELBOChk": training_options["freqELBOChk"],
        "drop_factor_threshold": training_options["drop_factor_threshold"],
        "scale_views": data_options["scale_views"],
        "scale_groups": data_options["scale_groups"],
        "center_groups": data_options["center_groups"],
        "G": G,
        "likelihoods": likelihoods,
        "script_start_time": script_start_time,
        "script_end_time": script_end_time
    }
    FctrMeta_out = os.path.join(OutputDir, "FctrMeta.json")

    with open(FctrMeta_out, "w") as outfile:
        json.dump(FctrMeta, outfile)

## PARAMETERS

## input
Prjct = "Trg_Normal_PN"
TopD = '5000D'  ## number of features maintained during variance filtering XD
Prior_K = 7 # how many factors to start with - 100 for Lrn, trying 30/60 for Trg Full, 100 for Trg multi

omics = ['mRNA','DNAme']
likelihoods = ['gaussian', 'gaussian']

M = len(omics)
G = 1 ## groups - always set to 1

data_options = {
    "scale_views": False,
    "scale_groups": False,
    "center_groups": True
}

## drop_factor_threshold: threshold for percentage of variance explained for dropping factors during training
## threshold of 0.001 for Lrn and 0.01 for Trg Full
training_options = {
    "TrainingIter": 10000, ## max training iterations
    "freqELBOChk": 5, ## how often to check the elbo, this is the default for R so using this
    "drop_factor_threshold": 0.01, ## am trying both 001 and 01 for usecase Lrn
    "seed": 1234567
}

## DIRECTORIES
## set up paths and directories
BasePath = os.getcwd()
InputDir = os.path.join(BasePath, Prjct + '_' + TopD,'MOFA')
# InputDir = os.path.join(BasePath, Prjct + '_' + TopD)
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
