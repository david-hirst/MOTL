## MT - 20240227

## SET WORK DIRECTORY
# setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/MOTL/TCGAStudy/TCGA_analysisTests/")
setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/MOLTI_TL-VI/MOTL/TCGAStudy/TCGA_analysisTests/")

## SET FUNCTIONS
source("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/MOLTI_TL-VI/MOTL/TCGAStudy/TCGA_preprocessedData_functions.R")
source("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/MOLTI_TL-VI/MOTL/TCGAStudy/TL_VI_functions.R")

## LIBRARIES
library("MOFA2")
library("rhdf5")
library("SummarizedExperiment")
library(dplyr)

## DATA PARAMETERS
Seed = 1234567
TopD = 5000
LrnK = 100
# viewsTrg = c('mRNA','miRNA','DNAme', 'SNV')
viewsTrg = c('mRNA','miRNA','DNAme')
# Prjcts = c("TCGA-LAML", "TCGA-PAAD", "TCGA-SKCM")
# Prjcts = "TCGA-PAAD"
Prjcts = c("TCGA-LAML", "TCGA-SKCM")
SS_size = 5
YTrgFull <- list()

## FOLDER PARAMETERS
DataDir = "DownloadedData_unfltrd/"

LrnDir = paste0('Lrn_',TopD,'D')
LrnFctrnDir = file.path(LrnDir,paste0('Fctrzn_',LrnK,'K'))

TrgPrjcts <- paste0(substr(Prjcts,6,nchar(Prjcts)), collapse='_')
TrgBrcd <- paste0('Trg_', TrgPrjcts, "_Full_",TopD,'D/expdat_meta.rds')
TrgDirSS <- paste0('Trg_', TrgPrjcts, '_SS',SS_size,'_',TopD,'D')
# TrgDir <- paste0('Trg_', TrgPrjcts, '_',TopD,'D')

DataDir = "../TCGA_database/DownloadedData_unfltrd/"
LrnDir <- "../TCGA_database/MOFA_TL/Lrn_5000D/"
LrnFctrnDir <- "../TCGA_database/MOFA_TL/Lrn_5000D/Fctrzn_100K/"
TrgBrcd <- "Trg_PAAD_5000D/expdat_meta.rds"
TrgDirSS <- "Trg_PAAD_SS5_5000D/"


## TL PARAMETERS
## LrnSimple = TRUE ## if TRUE then E[W^2] and E[LnTau] are calculated form E[W] and E[Tau] otherwise imported
CenterTrg = FALSE # Center Trg with own means or use estimated intercepts?

if (CenterTrg){
  TLDirName = 'TL_VI_C'
} else {
  TLDirName = 'TL_VI'
}

# ------------------------------------------- LEARNING SET ----

script_start_time = Sys.time()

# IMPORT FACTORIZATION RESULTS

## METADATA
expdat_meta_Lrn = readRDS(file.path(LrnDir,"expdat_meta.rds"))

## FATORIZATION MODEL OBJECT
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

# -------------------------------------------------------------

# --------------------------------------------- TARGET SET ----

## harmonize the vectors of views and likelihoods
views = viewsLrn[is.element(viewsLrn,viewsTrg)]
likelihoods = likelihoodsLrn[views]

### Import omics data for all Trg set samples
YTrgFull_brds <- readRDS(TrgBrcd)
Prjcts <- unique(YTrgFull_brds$brcds_mRNA$prjct)

## Prepare YTrgFull
YTrgFull$mRNA <- mRNAdataPreparation(Prjcts = Prjcts, brcds_mRNA = YTrgFull_brds$brcds_mRNA, InputDir = DataDir)
YTrgFull$miRNA <- miRNAdataPreparation(Prjcts = Prjcts, brcds_miRNA = YTrgFull_brds$brcds_miRNA, InputDir = DataDir)
YTrgFull$DNAme <- DNAmedataPreparation(Prjcts = Prjcts, brcds_DNAme = YTrgFull_brds$brcds_DNAme, InputDir = DataDir)
# YTrgFull$SNV <- SNVdataPreparation(Prjcts = Prjcts, brcds_SNV = YTrgFull_brds$brcds_SNV, InputDir = DataDir)

# -------------------------------------------------------------

# -------------------------------------- TRANSFER LEARNING ----

## PARAMETERS
minFactors = 15 ## floor when dropping factors - number of samples in evlauations
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
brcds_SS = readRDS(file.path(TrgDirSS,'brcds_SS.rds'))

## loop to subset, preprocess, perform TL and export
loop_start_time = Sys.time()

## FOR EACH TARGET SUBSET
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
  
  ## extract sample ids for the subset
  smpls = brcds_SS$smpls_SS[[SS]]
  
  ## Prepare YTrg subset data 
  YTrgSS = TCGATargetDataPreparation(views = views, YTrgFull = YTrgFull, brcds_SS = brcds_SS, 
                                     SS = SS, Fctrzn = Fctrzn, smpls = smpls, expdat_meta_Lrn = expdat_meta_Lrn)
  
  ## INIT PARAMETERS
  TL_param = initTransferLearningParamaters(views = views, YTrg = YTrgSS, Fctrzn = Fctrzn, 
                                            likelihoods = likelihoods, expdat_meta_Lrn = expdat_meta_Lrn)
  
  ## TRANSFER LEARNING
  TL_output <- transferLearning_function(TL_param = TL_param, MaxIterations= MaxIterations, MinIterations = MinIterations, minFactors = minFactors, 
                            StartDropFactor = StartDropFactor, FreqDropFactor = FreqDropFactor, StartELBO = StartELBO, 
                            FreqELBO = FreqELBO, DropFactorTH = DropFactorTH, ConvergenceIts = ConvergenceIts, ConvergenceTH = ConvergenceTH, 
                            CenterTrg = CenterTrg, outputDir = TL_SSOutDir)
  
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
    'PoisRateCstnt' = TL_output$PoisRateCstnt,
    'LrnSimple' = TRUE,
    'script_start_time' = script_start_time,
    'loop_start_time' = loop_start_time,
    'script_end_time' = script_end_time
  )
  saveRDS(TL_meta_data, file.path(TrgDirSS,paste0(TLDirName,'_meta_data.rds')))
  
}

# -------------------------------------------------------------