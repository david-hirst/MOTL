#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues April 9 2024

@author: Morgane Térézol
"""

## Load require libraries

import os
import numpy as np
from mofapy2.run.entry_point import entry_point
import json

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
