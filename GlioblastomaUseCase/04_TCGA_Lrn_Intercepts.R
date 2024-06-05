
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
FctrznDir = file.path(ExpDataDir,'Fctrzn_100K_001TH')

## import meta data and model

expdat_meta = readRDS(file.path(ExpDataDir,'expdat_meta.rds'))

InputModel = file.path(FctrznDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

# loop through the views and estimate the intercept
# for gaussian data its just the mean
# for other data will try mle with a naive estimator as backup

intercepts_calculation(expdat_meta = expdat_meta, 
                       Fctrzn = Fctrzn,
                       FctrznDir = FctrznDir, 
                       ExpDataDir = ExpDataDir)


###
### some testing
###

# EstimatedIntercepts = readRDS(file.path(FctrznDir,"EstimatedIntercepts.rds"))
# Intercepts = EstimatedIntercepts$Intercepts
# InterceptsMethod = EstimatedIntercepts$InterceptsMethod
# InterceptsNaive = EstimatedIntercepts$InterceptsNaive
# 
# 
# hist(Intercepts$mRNA)
# sum(is.na(Intercepts$mRNA))
# 
# hist(Intercepts$miRNA)
# sum(is.na(Intercepts$miRNA))
# 
# hist(Intercepts$DNAme)
# sum(is.na(Intercepts$DNAme))
# 
# table(InterceptsMethod$SNV)
# sum(is.na(Intercepts$SNV))
# 
# rnd_idx = sample(1:length(InterceptsNaive$SNV),1000,replace=FALSE)
# plot(x=InterceptsNaive$SNV[rnd_idx],y=Intercepts$SNV[rnd_idx])
# abline(a=0,b=1)
