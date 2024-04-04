# MOTL: multi-omics matrix factorization with transfer learning

## MOTL code repository contents

This repository contains code used to recreate analyses described in the MOTL paper, as well as functions that can be downloaded to apply MOTL to a target dataset in order to perform matrix factorization with transfer learning.

The repository is organised into three folders:

[SimulationStudy](https://github.com/david-hirst/MOTL/tree/main/SimulationStudy): This folder contains scripts for carrying out the simulation study. More details in the [00_SimStudy_ReadMe.md](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/00_SimStudy_ReadMe.md) file

[TCGAStudy](https://github.com/david-hirst/MOTL/tree/main/TCGAStudy): This folder contains scripts for carrying out the TCGA study, as well as functions that can be downloaded to apply MOTL to a target dataset of interest. More details in the [00_TCGAstudy_ReadMe.md](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/00_TCGAstudy_ReadMe.md) file.

[GlioblastomaUseCase](https://github.com/david-hirst/MOTL/tree/main/GlioblastomaUseCase): This folder contains scripts for the use case illustration of MOTL with glioblastoma data. MOre details in the [00_GlioblastomaUseCase_ReadMe.md](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/00_GlioblastomaUseCase_ReadMe.md) file.

## How to use MOTL on your own data 

MOTL is run using R. 

### Downloads
You will firstly need to download the following files of functions from this github site:
[TCGA_preprocessedData.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_preprocessedData.R)
[TL_VI_functions.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TL_VI_functions.R)

Next download the learning dataset factorization files from [this zenodo repository](https://zenodo.org/records/10848217) and unzip.

### load libraries and functions
```
library(SummarizedExperiment)
library(DESeq2)
library(sesame)
library(dplyr)
library(MOFA2)
library(rhdf5)

source('TCGA_preprocessedData_functions.R')
source('TL_VI_functions.R')
```
### Pre-processing of learning dataset factorization
Location of learning data downloaded from zenodo
```
LrnDir = 'LrnData'
LrnFctrnDir = file.path(LrnDir,'Lrn_5000D_Fctrzn_100K_001TH')
```
Specify whether to center the target dataste during processing or to leave uncentered and use enstmated learning dataset intercepts
```
CenterTrg = FALSE
```
Import the learning dataset metadata, factorization, and initialise values
```
expdat_meta_Lrn = readRDS(file.path(LrnDir,"expdat_meta.rds"))
InputModel = file.path(LrnFctrnDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

viewsLrn = Fctrzn@data_options$views
likelihoodsLrn = Fctrzn@model_options$likelihoods
MLrn = Fctrzn@dimensions$M

Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)
Fctrzn@expectations[["TauLn"]] = sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnDir)
Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir)
Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg, Fctrzn, LrnFctrnDir)
```
### Pre-processing of target dataset factorization
Record the time that the pre-processing starts
```
ss_start_time = Sys.time()
```

Start with some or all of the following omics matrix, all with features in rows and samples in columns. The names of these matrices are not important. The set of column names should be the same for each omics, although the order is not important as they will be ordered automatically. The order of the features is not important, however the type of id used to name the features should be consistant with the TCGA learning dataset, as outlined below.

**expdat_mRNA**: a matrix of mRNA raw counts, genes in rows, samples in columns. Row names should be Ensemble ids without the version suffix; for example *ENSG00000000005*.

**expdat_miRNA**: a matrix of miRNA raw counts, miRNAs in rows, samples in columns. Row names should be miRNA names as per miRBase (v21); for example *hsa-mir-1-1*.

**expdat_DNAme**: a matrix of DNA methylation M-values, cpgs in rows, samples in columns. Row names should be cpg probe ids from either the 450 or epic illumina array; for example *cg09364122*.

**expdat_SNV**: a binary matrix of SNV mutation absence / presence, genes in rows, samples in columns. Row names should be Hugo symbols; for example *AKAP13*.

In the case of expdat_mRNA, a version suffix will need to be added to be consisant with the naming in the TCGA learning data factorization.
```
expdat_mRNA = mRNA_addVersion(expdat = expdat_mRNA, Lrndat = Fctrzn@expectations$W$mRNA)
```

Create a list of the matrices, using the following naming convention for each omics.
```
YTrg_list = list(
  mRNA = expdat_mRNA,
  miRNA = expdat_miRNA,
  DNAme = expdat_DNAme,
  SNV = expdat_SNV
)
```
Initialise values for transfer learning
```
smpls = colnames(YTrg_list[[1]])
viewsTrg = names(YTrg_list)
views = viewsLrn[is.element(viewsLrn,viewsTrg)]
likelihoods = likelihoodsLrn[views]

YTrg_list = TargetDataPreparation(views = views, YTrg_list = YTrg_list, 
                                  Fctrzn = Fctrzn, smpls = smpls, expdat_meta_Lrn = expdat_meta_Lrn,
                                  normalization = 'Lrn', transformation = TRUE)

TL_param = initTransferLearningParamaters(YTrg = YTrg_list, views = views, 
                                          expdat_meta_Lrn = expdat_meta_Lrn, 
                                          Fctrzn = Fctrzn, likelihoods = likelihoods
                                          )
```

### Transfer learning factorization with MOTL

Specify output folder and paramaters for MOTL
```
TL_OutDir = 'MOTL_Fctrzn'
if(!dir.exists(TL_OutDir)){
  dir.create(TL_OutDir, showWarnings = FALSE, recursive = TRUE)
}

minFactors = 6 ## floor when dropping factors
StartDropFactor = 1 # after which iteration to start dropping factors
FreqDropFactor = 1 # how often to drop factors
StartELBO = 1 # which iteration to start checking ELBO on, excl initiation iteration
FreqELBO = 5 # how often to assess the ELBO
DropFactorTH = 0.01 # factor with lowest max variance, that is less than this, is dropped
MaxIterations = 10000
MinIterations = 2 # 
ConvergenceIts = 2 # numbe rof consectutive checks in a row for which the change in elbo is below the threshold
ConvergenceTH = 0.0005 # change in elbo threshold
```
run MOTL to infer and save the factorization as an rds file
```
TL_data = transferLearning_function(TL_param = TL_param, MaxIterations = MaxIterations, MinIterations =  MinIterations, 
                                    minFactors = minFactors, StartDropFactor = StartDropFactor, FreqDropFactor = FreqDropFactor, 
                                    StartELBO = StartELBO, FreqELBO = FreqELBO, DropFactorTH = DropFactorTH, 
                                    ConvergenceIts = ConvergenceIts, ConvergenceTH = ConvergenceTH, 
                                    CenterTrg = CenterTrg, outputDir = TL_OutDir)
```
Extract the Z and W matrices from the MOTL factrization. Z is the inferred score matrix for the target dataset: rows and samples, and columns are factors. Each W has features in the columns and factors in the rows. Factor names correspond to the name form the learning dataset factorization
```
ZMu = TL_data$ZMu
W_mRNA = TL_data$Fctrzn_Lrn_W$mRNA
```




