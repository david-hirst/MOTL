
################
#### Import factorization and input data and learn intercepts
###############

## MOFA doesn't learn intercepts but it could be good to have them
## this script can be used to learn intercepts via MLE

library(ggplot2)
library(MOFA2)
library(rhdf5)
library(rjson)

source("TL_VI_functions.R")

Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

### Set directories

ExpDataDir = 'Lrn_5000D'
FctrznDir = file.path(ExpDataDir,'Fctrzn_100K')

## import meta data and model

expdat_meta = readRDS(file.path(ExpDataDir,'expdat_meta.rds'))

InputModel = file.path(FctrznDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

intercepts_calculation(expdat_meta = expdat_meta, 
                       Fctrzn = Fctrzn,
                       FctrznDir = FctrznDir, 
                       ExpDataDir = ExpDataDir)

