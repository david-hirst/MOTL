#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 11:05:07 2022

@author: david
"""

# This is a function for simulating multi-omics data from latent factors 
# loosely based on the generative model from the MOFA paper
# various paramaters and hyperparamters can be changed
# A score matrix and M weight matrices are generated with guassian distirbutions
# Binary matrices are generated for each omics - elementwise multiplication with guassian distributed weight matrices produces final weight matrices
# The product of each weight matrix and the score matrix produces a matrix of parameter values from which the observed data is simulated based on specified likelihood


import numpy as np
import sys

## define function for simulating a multi-omics dataset
   
def SimMultiOmics(CLrn = 20, CLrnNPool = [10,25,40], CTrg = 2, CTrgN = 7,
                  K = 30, ZMuBase = 5, ZSdRange = [1,1], ZMuShift = [-2,2], ZShiftP = 1/4,
                  D = [1000,1000,1000,1000], YDists = ["Poisson", "Normal", "Beta", "Bernoulli"], 
                  WThetaRange = [0.1, 0.2],
                  WMuPois = 5, WSdRangePois = [0.5, 1.5],
                  WMuNorm = 0, WSdRangeNorm = [0.5, 1.5],
                  WMuBeta = 0, WSdRangeBeta = [0.25, 0.75], 
                  WMuBern = 0, WSdRangeBern = [0.5, 1.5],
                  YSdRangeNorm = [0.25, 0.75], YaRangeBeta = [50, 150]):
    
    ### Explanation of parameters
    # CLrn: number of learning dataset sample groups
    # CLrnNPool: pool of sizes for learning dataset sample groups
    # CTrg: number of target dataset sample groups
    # CTrgN: size of each target dataset sample group
    # K: number of latent factors
    # ZMuBase: Base mean parameter for factor sample scores
    # ZSdRange: range from which to choose standard deviation for sample scores for each factor
    # ZMuShift: potential shifts in mean parameter for sample factor scores
    # ZShiftP: probability that score mean parameter shifts for a particular sample group and factor
    # D: number of features for each omic
    # YDists: observed data likelihoods for each omic
    # WThetaRange: parameter for generating bernoulli random variables to induce sparsity in weight matrices
    # WMuPois, WMuNorm, etc : the mean parameter for simulating weight matrices depending on the corresponding observed data likelihood
    # WSdRangePois, WSdRangeNorm, etc : the standard deviation parameter for simulating weight matrices depending on the corresponding observed data likelihood
    # YSdRangeNorm: range from which to select featurewise standard deviations for generating guassian observed data
    # YaRangeBeta: range from which to select featurewise alpha parameter values when generating beta distributed observed data
    

    ########## Check feature counts and distributions
    
    # number of views
    M = len(D)
    
    if M != len(YDists) :
        print("The vectors D and YDists must be the same length")
        sys.exit()
     
    for m in range(M): 
        if YDists[m]!="Poisson" and YDists[m]!="Normal" and YDists[m]!="Beta" and YDists[m]!="Bernoulli":
            print("Distributions must be either Poisson, Normal, Beta or Bernoulli")
            sys.exit()
            
             
    ########## Create Sample clusters (groups)   
     
    # combined learning and target sets
    C = CLrn + CTrg
    Nc = np.append(np.random.choice(a = CLrnNPool, size = CLrn),np.repeat(CTrgN,CTrg))
    # assign sample groupings
    CIds = np.repeat(a = range(C), repeats = Nc)
    IfLrn = CIds <= max(range(CLrn))
    
    
    ######### Simulate score values
    
    # pool base and shift values
    ZMuPool = np.append(ZMuBase,ZMuBase+np.array(ZMuShift))
    # sampling probabilities for the pooled values
    ZSampP = np.append(1-ZShiftP, np.repeat(ZShiftP/len(ZMuShift),len(ZMuShift)))
    
    # empty lists for parameter values and sampled scores
    ZMuC = []
    ZMuS = []
    Z = []
    
    # cluster means
    for k in range(K):
        ZMuC.append(np.random.choice(a = ZMuPool, size = C, p = ZSampP))
        
    # sample means based on corresponding cluster value
    for k in range(K):
        ZMuS.append(np.repeat(a = ZMuC[k], repeats = Nc))
        
    # generate scores for each sample
    for k in range(K):
        Z.append(np.random.normal(loc = ZMuS[k], 
                                  scale = np.random.uniform(low=ZSdRange[0], high=ZSdRange[1], size=1), 
                                  size=len(ZMuS[k])))
          
    ZMuC = np.vstack(ZMuC).transpose()   
    Z = np.vstack(Z).transpose()
    
    ######### Simulate weight matrices
    
    W = []        
    
    for m in range(M):
        if YDists[m] == "Poisson":
            W.append(
                np.column_stack(
                    [np.random.normal(loc = WMuPois, 
                                      scale = np.random.uniform(low=WSdRangePois[0], high=WSdRangePois[1], size=1), 
                                      size = D[m]) * 
                     np.random.binomial(n = 1, 
                                        p = np.random.uniform(low=WThetaRange[0], high=WThetaRange[1], size=1), 
                                        size = D[m]) 
                     for k in range(K)]
                    )
                )
        elif YDists[m] == "Normal":
            W.append(
                np.column_stack(
                    [np.random.normal(loc = WMuNorm, 
                                      scale = np.random.uniform(low=WSdRangeNorm[0], high=WSdRangeNorm[1], size=1), 
                                      size = D[m]) * 
                     np.random.binomial(n = 1, 
                                        p = np.random.uniform(low=WThetaRange[0], high=WThetaRange[1], size=1), 
                                        size = D[m]) 
                     for k in range(K)]
                    )
                )
        elif YDists[m] == "Beta":
            W.append(
                np.column_stack(
                    [np.random.normal(loc = WMuBeta, 
                                      scale = np.random.uniform(low=WSdRangeBeta[0], high=WSdRangeBeta[1], size=1), 
                                      size = D[m]) * 
                     np.random.binomial(n = 1, 
                                        p = np.random.uniform(low=WThetaRange[0], high=WThetaRange[1], size=1), 
                                        size = D[m]) 
                     for k in range(K)]
                    )
                )
        elif YDists[m] == "Bernoulli":
            W.append(
                np.column_stack(
                    [np.random.normal(loc = WMuBern, 
                                      scale = np.random.uniform(low=WSdRangeBern[0], high=WSdRangeBern[1], size=1), 
                                      size = D[m]) * 
                     np.random.binomial(n = 1, 
                                        p = np.random.uniform(low=WThetaRange[0], high=WThetaRange[1], size=1), 
                                        size = D[m]) 
                     for k in range(K)]
                    )
                )
        else: 
            print("Distributions must be either Poisson, Normal, Beta or Bernoulli")
            sys.exit()
            
    ##### Simulate Data
        
    ## Calculate signal
    
    S = []
    for m in range(M):
        S.append(Z @ W[m].T)
        
    ## Use the signal to simulate the data
    
    Y = []
    
    for m in range(M):
        if YDists[m] == "Poisson":
            Y.append(
                np.column_stack(
                    [np.random.poisson(np.log(1+np.exp(S[m][:,d])))
                     for d in range(D[m])]
                    )            
                )
        elif YDists[m] == "Normal":
            Y.append(
                np.column_stack(
                    [np.random.normal(loc = S[m][:,d],
                                      scale = np.random.uniform(low=YSdRangeNorm[0], high=YSdRangeNorm[1], size=1)
                                      ) for d in range(D[m])]
                    )
                )
        elif YDists[m] == "Beta":
            Yalphas = np.random.uniform(low=YaRangeBeta[0], high=YaRangeBeta[1], size=D[m])
            Y.append(
                np.column_stack(
                    [np.random.beta(a = Yalphas[d], 
                                    b = Yalphas[d] / np.power(1+np.exp(-S[m][:,d]),-1) - Yalphas[d] + 0.1
                                    ) for d in range(D[m])]
                    )
                )
        elif YDists[m] == "Bernoulli":
            Y.append(
                np.column_stack(
                    [np.random.binomial(n=1,
                                        p = np.power(1+np.exp(-S[m][:,d]),-1)
                                        ) for d in range(D[m])]
                    )            
                )
        else: 
            print("Distributions must be either Poisson, Normal, Beta or Bernoulli")
            sys.exit()
            
    return {'Nc': Nc, 'CIds': CIds, 'IfLrn': IfLrn, 'ZMuC': ZMuC, 'W': W, 'Y': Y, 'Z': Z}

    ## this returns a dictionary containing the following: 
    # Nc: The number of samples in each group of samples. The learning dataset groups are first
    # CIds: The group membership for each sample
    # IfLrn: Whether a not each sample is in the learning dataset
    # ZMuC: Matrix containing sample score means by sample group (rows) and factor (column)
    # W: list of weight matrices, in same order as the likelihoods specified for observed data. Each matrix has features weights by feature (rows) and factor (columns)
    # Y: list of omics data matrices, in same order as the likelihoods specified for observed data. Each matrix has an observed data point by sample (rows) and feature (column)
    
        