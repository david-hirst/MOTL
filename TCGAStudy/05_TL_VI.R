## MT - 20240227

#############
## APPLY MOTL TO TCGA MULTI-OMICS TARGET DATASETS
#############

# This script can be used to apply MOTL to the TARGET datasets created from a particular combination of TCGA projects
# The result will be a factorization of each TARGET dataset

## IMPORT FUNCTIONS
source("TCGA_preprocessedData_functions.R")
source("TL_VI_functions.R")

## LIBRARIES
library("MOFA2")
library("rhdf5")
library("SummarizedExperiment")
library("dplyr")

## DATA PARAMETERS
Seed = 1234567
TopD = 5000
LrnK = 100
viewsTrg = c('mRNA','miRNA','DNAme') # The omics in the TARGET dataset
Prjcts = c("TCGA-LAML", "TCGA-SKCM") # The TCGA projects from which the TARGET datasets are created
SS_size = 5 # The subset size, which is the number of samples per project in each TARGET dataset
YTrgFull <- list() # placeholder for a multi-omics dataset for all samples that will be subsetted to create TARGET datasets

## FOLDER PARAMETERS
DataDir = "DownloadedData_unfltrd/" # location of downloaded TCGA data

LrnDir = paste0('Lrn_',TopD,'D') # When the LEARNING dataset data files are save
LrnFctrnDir = file.path(LrnDir,paste0('Fctrzn_',LrnK,'K')) # where the LEARNING dataset MOFA factorization is saved

TrgPrjcts <- paste0(substr(Prjcts,6,nchar(Prjcts)), collapse='_')
TrgBrcd <- paste0('Trg_', TrgPrjcts, "_Full_",TopD,'D/expdat_meta.rds')
TrgDirSS <- paste0('Trg_', TrgPrjcts, '_SS',SS_size,'_',TopD,'D')

## TL PARAMETERS
## LrnSimple = TRUE ## if TRUE then E[W^2] and E[LnTau] are calculated form E[W] and E[Tau] otherwise imported
CenterTrg = FALSE # Center the TARGET datasets with own featurewise means or use estimated intercepts from the LEARNING dataset factorization?

if (CenterTrg){
  TLDirName = 'TL_VI_C'
} else {
  TLDirName = 'TL_VI'
}

##############################################################################
# ------------------------------------------- IMPORT LEARNING DATASET DATA----
##############################################################################

script_start_time = Sys.time()

# IMPORT FACTORIZATION RESULTS

## METADATA
expdat_meta_Lrn = readRDS(file.path(LrnDir,"expdat_meta.rds"))

## FACTORIZATION MODEL OBJECT
InputModel = file.path(LrnFctrnDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

viewsLrn = Fctrzn@data_options$views
likelihoodsLrn = Fctrzn@model_options$likelihoods
MLrn = Fctrzn@dimensions$M

## LOAD IN THE EXPECTATIONS FOR TAU AND REORDER FOR CONSISTENCY
Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)

## E[log(tau)] CALCULATION (ONLY RELEVANT FOR GAUSSIAN DATA)
Fctrzn@expectations[["TauLn"]] = sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnDir)

## load or calculate E[W^2] values
Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir)

## import the W intercepts if not centering, create vectors containing 0s if centering
Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg, Fctrzn, LrnFctrnDir)

##############################################################################
# -----------------------DATA FOR ALL SAMPLES FROM TARGET DATASET PROJECTS----
##############################################################################

## harmonize the vectors of views and likelihoods
views = viewsLrn[is.element(viewsLrn,viewsTrg)]
likelihoods = likelihoodsLrn[views]

### Import omics data for all TARGET dataset samples
YTrgFull_brds <- readRDS(TrgBrcd)
Prjcts <- unique(YTrgFull_brds$brcds_mRNA$prjct)

## Prepare YTrgFull, which will be subsetted to create TARGET datasets
YTrgFull$mRNA <- mRNAdataPreparation(Prjcts = Prjcts, brcds_mRNA = YTrgFull_brds$brcds_mRNA, InputDir = DataDir)
YTrgFull$miRNA <- miRNAdataPreparation(Prjcts = Prjcts, brcds_miRNA = YTrgFull_brds$brcds_miRNA, InputDir = DataDir)
YTrgFull$DNAme <- DNAmedataPreparation(Prjcts = Prjcts, brcds_DNAme = YTrgFull_brds$brcds_DNAme, InputDir = DataDir)
# YTrgFull$SNV <- SNVdataPreparation(Prjcts = Prjcts, brcds_SNV = YTrgFull_brds$brcds_SNV, InputDir = DataDir)

##############################################################################
# --------------------CREATE AND PREPROCESS TARGET DATASETS AND APPLY MOTL----
##############################################################################

## PARAMETERS
minFactors = 10 ## floor when dropping factors - we used the number of samples
StartDropFactor = 1 # after which iteration to start dropping factors
FreqDropFactor = 1 ## how often to drop factors
StartELBO = 1 #which iteration to start checking ELBO on, excl initiation iteration
FreqELBO = 5 #how often to assess the ELBO
DropFactorTH = 0.01 # factor with lowest max variance, that is less than this, is dropped
MaxIterations = 3
MinIterations = 2 # as per MOFA with defaults in python, >= 2 (hardcoded floor), exl initial setup
ConvergenceIts = 2 # as per MOFA python defaults
ConvergenceTH = 0.0005 # default for MOFA - fast option

## import predefined subset information
## these are the samples already allocated to each TARGET dataset 
brcds_SS = readRDS(file.path(TrgDirSS,'brcds_SS.rds'))

## loop to subset YTrgFull to create, pre-process and apply MOTL to TARGET datsets
loop_start_time = Sys.time()

## FOR EACH TARGET DATASET
for(ss in 1:brcds_SS$SS_count){
  
  ss_start_time = Sys.time()
  SS = paste0('SS_',ss)
  print(SS)
  
  ## directories to export data to
  SSOutDir = file.path(TrgDirSS, SS)
  TL_SSOutDir = file.path(SSOutDir, TLDirName)
  if(!dir.exists(TL_SSOutDir)){
    dir.create(TL_SSOutDir, showWarnings = FALSE, recursive = TRUE)
  }
  
  ## extract sample ids for the TARGET dataset
  smpls = brcds_SS$smpls_SS[[SS]]
  
  ## select a subset of YTrgFull to create the TARGET dataset
  YTrgSS = TCGATargetDataPreparation(views = views, YTrgFull = YTrgFull, brcds_SS = brcds_SS, 
                                     SS = SS, Fctrzn = Fctrzn, smpls = smpls, expdat_meta_Lrn = expdat_meta_Lrn)
  
  ## INIT PARAMETERS
  TL_param = initTransferLearningParamaters(views = views, YTrg = YTrgSS, Fctrzn = Fctrzn, 
                                            likelihoods = likelihoods, expdat_meta_Lrn = expdat_meta_Lrn)
  
  ## TRANSFER LEARNING WITH MOTL
  TL_output <- transferLearning_function(TL_param = TL_param, MaxIterations= MaxIterations, MinIterations = MinIterations, minFactors = minFactors, 
                            StartDropFactor = StartDropFactor, FreqDropFactor = FreqDropFactor, StartELBO = StartELBO, 
                            FreqELBO = FreqELBO, DropFactorTH = DropFactorTH, ConvergenceIts = ConvergenceIts, ConvergenceTH = ConvergenceTH, 
                            CenterTrg = CenterTrg, ss_start_time = ss_start_time, outputDir = TL_SSOutDir)
  
  ## save overall meta data
  script_end_time = Sys.time()
  
  TL_meta_data = list(
    'Seed' = Seed,
    'CenterTrg' = CenterTrg,
    'StartDropFactor' = StartDropFactor,
    'FreqDropFactor' = FreqDropFactor,
    'DropFactorTH' = DropFactorTH,
    'minFactors' = minFactors,
    'MinIterations' = MinIterations,
    'MaxIterations' = MaxIterations,
    'StartELBO' = StartELBO,
    'FreqELBO' = FreqELBO,
    'ConvergenceTH' = ConvergenceTH,
    'ConvergenceIts' = ConvergenceIts,
    # 'PoisRateCstnt' = 0.0001,
    # 'LrnSimple' = TRUE,
    'script_start_time' = script_start_time,
    'loop_start_time' = loop_start_time,
    'script_end_time' = script_end_time
  )
  saveRDS(TL_meta_data, file.path(TrgDirSS,paste0(TLDirName,'_meta_data.rds')))
  
}

# -------------------------------------------------------------
