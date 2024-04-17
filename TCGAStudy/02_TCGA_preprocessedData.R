## MT - 20231110

#############
## CREATE AND PRE-PROCESS TCGA LEARNING, REFERENCE, AND TARGET MULTI-OMICS DATASETS
#############

## This script is divided into 3 main sections
## Sections are run independently, depending on the desired task

## 01 LEARNING DATASET
## This section has the code for creating and pre-processing a TCGA multi-omics learning dataset
## The outputs include metadata and a csv file for each omics matrix
## The outputs are used as input for a MOFA factorization of the learning dataset

## 02 REFERENCE AND TARGET DATASETS FROM ONE PROJECT
## This section has the code for creating and pre-processing single-project TCGA multi-omics reference, and target, datasets
## The TCGA project and the required omics are specified 
## Then sample barcodes are selected - 1 profile per sample submittor, only samples with full multi-omics profiles
## The reference dataset is created from the barcodes, and then pre-processed. The outputs can be used as input to MOFA factorization
## The barcodes are then subsetted randomly. Each subset of barcodes is used to create a target dataset.
## Target datasets are pre-processed, and the output can be used as input to MOFA factorization.
## The full set of barcodes, and the barcode subsets are saved to be used in the MOTL application and for evaluation purposes


## 03 REFERENCE AND TARGET DATASETS FROM MULTIPLE PROJECTS
## Same as previous section, but the reference dataset, and each target dataset, contain samples from multiple TCGA projects


## LIBRARIES
library(TCGAbiolinks)
library(SummarizedExperiment)
library(SEtools)
library(sesame)
library(dplyr)

## FUNCTIONS
source("TCGA_preprocessedData_functions.R")

## ------------------------------------------------------------------------------
## ------------------------------------------------------------ 01 LEARNING DATASET ----

print(Sys.time())
print("Learning set")

## ---------------- SET ENVIRONMENT ----

## SET PARAMETERS 
PrjctExcl = c("TCGA-PAAD","TCGA-LAML","TCGA-SKCM","TCGA-GBM") # projects to exclude
InputDir = "DownloadedData_unfltrd" # where the downloaded TCGA data was saved
BarcodesDir = 'Lrn_barcodes' # where the barcodes (sample ids) for the learning dataset will be saved
GeoMeans = "Lrn"
TopD = 5000 ## how many features to retain

## SET SEED
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

## PROJECTS TO INCLUDE IN THE LEARNING DATASET
Prjcts = selectProjects(InputDir = InputDir, PrjctExcl = PrjctExcl)
## -------------------------------------

## ----------------- SELECT SAMPLES ----

print(Sys.time())
print("Select samples")

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
if(!dir.exists(BarcodesDir)){
  dir.create(BarcodesDir, showWarnings = FALSE, recursive = TRUE)
}

## FOR EACH PROJECT
brcds_list <- sapply(X = Prjcts, FUN = function(Prjct){

  print("[*]-------------------------------------")
  print(paste("Project", Prjct))

  # import rds files
  print("Read omics files")
  expdat_mRNA = readRDS(file.path(InputDir,paste0(Prjct,"_mRNA.rds")))
  expdat_miRNA = readRDS(file.path(InputDir,paste0(Prjct,"_miRNA.rds")))
  expdat_DNAme = readRDS(file.path(InputDir,paste0(Prjct,"_DNAme.rds")))
  expdat_SNV = readRDS(file.path(InputDir,paste0(Prjct,"_SNV.rds")))

  # tidy up and filter omics data (parsing file)
  print("Parse omics data")
  expdat_mRNA_fltr = mRNAFileParsing(expdat_mRNA)
  expdat_miRNA_fltr = miRNAFileParsing(expdat_miRNA)
  expdat_DNAme_fltr = DNAmeFileParsing(expdat_DNAme)
  expdat_SNV_fltr = SNVFileParsing(expdat_SNV)

  # samples with shared omics
  print("Select shared samples between omics")
  omics_list <- list("mRNA" = expdat_mRNA_fltr, "miRNA" = expdat_miRNA_fltr, "DNAme" = expdat_DNAme_fltr, "SNV" = expdat_SNV_fltr)
  smpls_shared <- extractSampleNamesShared(omics_list)

  # for each submittor select the min barcode from those related to samples with shared omics
  print("Select the min barcode for each submittor")
  brcds_list_tmp = sapply(X = omics_list, FUN = function(expdat, smpls_shared, Prjct){
    brcds = selectMinBarcodeSampleNames(expdat = expdat, smpls_shared = smpls_shared, Prjct = Prjct)
    return(brcds)
  }, smpls_shared, Prjct, simplify = FALSE)

  print("[*]-------------------------------------")

  return(brcds_list_tmp)
})

## Merge barcodes
brcds_mRNA <- do.call(rbind, brcds_list["mRNA",])
brcds_miRNA <- do.call(rbind, brcds_list["miRNA",])
brcds_DNAme <- do.call(rbind, brcds_list["DNAme",])
brcds_SNV <- do.call(rbind, brcds_list["SNV",])

## Save barcodes
saveRDS(brcds_mRNA,file.path(BarcodesDir,'brcds_mRNA.rds'))
saveRDS(brcds_miRNA,file.path(BarcodesDir,'brcds_miRNA.rds'))
saveRDS(brcds_DNAme,file.path(BarcodesDir,'brcds_DNAme.rds'))
saveRDS(brcds_SNV,file.path(BarcodesDir,'brcds_SNV.rds'))

## -------------------------------------

## ------- CREATION OF LEARNING SET ----

print(Sys.time())
print("Learning set creation")

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
OutDir = paste0('Lrn_',TopD,"D")
if(!dir.exists(OutDir)){
  dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
}

## IMPORT BARCODES FROM PREVIOUS STEP
brcds_mRNA = readRDS(file.path(BarcodesDir,'brcds_mRNA.rds'))
brcds_miRNA = readRDS(file.path(BarcodesDir,'brcds_miRNA.rds'))
brcds_DNAme = readRDS(file.path(BarcodesDir,'brcds_DNAme.rds'))
brcds_SNV = readRDS(file.path(BarcodesDir,'brcds_SNV.rds'))

## TO COMMENT
## subset based on instersection of projects - in case only subset is desired for testing purposes
# Prjcts = Prjcts[c(1:5)]
# Prjcts = unique(base::intersect(Prjcts,brcds_mRNA$prjct))

brcds_mRNA = brcds_mRNA[is.element(brcds_mRNA$prjct,Prjcts),]
brcds_miRNA = brcds_miRNA[is.element(brcds_miRNA$prjct,Prjcts),]
brcds_DNAme = brcds_DNAme[is.element(brcds_DNAme$prjct,Prjcts),]
brcds_SNV = brcds_SNV[is.element(brcds_SNV$prjct,Prjcts),]

# save samples for shuffling and keeping track of names later on
smpls = substr(brcds_mRNA$brcds,1,16)
smpls = sample(smpls, length(smpls), replace = FALSE)

## ---------- mRNA data preprocessing

print(Sys.time())
print("mRNA")

## Read and merge data files
expdat_mRNA_merged <- mRNAdataPreparation(Prjcts = Prjcts, brcds_mRNA = brcds_mRNA, InputDir = InputDir)

## Filtering
expdat_mRNA_fltr <- mRNAFiltering(expdat_mRNA = expdat_mRNA_merged)

## Normalization
expdat_mRNA_norm <- countsNormalization(expdat = expdat_mRNA_fltr, GeoMeans = GeoMeans)
expdat_mRNA_geoMeans <- expdat_mRNA_norm$GeoMeans

## Transformation
expdat_mRNA_log <- countsTransformation(expdat_count = expdat_mRNA_norm$counts, TopD = TopD)

## Save data
expdat_mRNA_sorted <- orderAndSaveData(expdat = expdat_mRNA_log, smpls = smpls, OutDir = OutDir, fileName = "mRNA.csv")

## CLEAN SPACE
rm(expdat_mRNA_merged, expdat_mRNA_fltr, expdat_mRNA_norm, expdat_mRNA_log)
gc()

## ---------- miRNA data preprocessing 

print(Sys.time())
print("miRNA")

## Read and merge data files
expdat_miRNA_merged <- miRNAdataPreparation(Prjcts = Prjcts, brcds_miRNA = brcds_miRNA, InputDir = InputDir)

## Filtering
expdat_miRNA_fltr <- miRNAFiltering(expdat_miRNA = expdat_miRNA_merged)

## Normalization
expdat_miRNA_norm <- countsNormalization(expdat = expdat_miRNA_fltr, GeoMeans = GeoMeans)
expdat_miRNA_geoMeans <- expdat_miRNA_norm$GeoMeans

## Transformation
expdat_miRNA_log <- countsTransformation(expdat_count = expdat_miRNA_norm$counts, TopD = TopD)

## Save data
expdat_miRNA_sorted <- orderAndSaveData(expdat = expdat_miRNA_log, smpls = smpls, OutDir = OutDir, fileName = "miRNA.csv")

## CLEAN SPACE
rm(expdat_miRNA_merged, expdat_miRNA_fltr, expdat_miRNA_norm, expdat_miRNA_log)
gc()

## ---------- DNAme data preprocessing

print(Sys.time())
print("DNAme")

## Read and merge data files
expdat_DNAme_merged <- DNAmedataPreparation(Prjcts = Prjcts, brcds_DNAme = brcds_DNAme, InputDir = InputDir)

## Filtering
expdat_DNAme_fltr <- DNAmeFiltering(expdat_DNAme = expdat_DNAme_merged, TopD = TopD)

## Save data
expdat_DNAme_sorted <- orderAndSaveData(expdat = expdat_DNAme_fltr, smpls = smpls, OutDir = OutDir, fileName = "DNAme.csv")

## CLEAN SPACE
rm(expdat_DNAme_merged, expdat_DNAme_fltr)
gc()

## ---------- SNV data preprocessing

print(Sys.time())
print("SNV")

## Read and merge data files
expdat_SNV_merged <- SNVdataPreparation(Prjcts = Prjcts, brcds_SNV = brcds_SNV, InputDir = InputDir)

## Filtering
expdat_SNV_fltr <- SNVFiltering(expdat_SNV = expdat_SNV_merged, TopD = TopD)

## Save data
expdat_SNV_sorted <- orderAndSaveData(expdat = expdat_SNV_fltr, smpls = smpls, OutDir = OutDir, fileName = "SNV.csv")

## CLEAN SPACE
rm(expdat_SNV_merged, expdat_SNV_fltr)
gc()

## -------------------------------------

## ------------------ SAVE METADATA ----

## SAVE METADATA
expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted, "SNV" = expdat_SNV_sorted)
brcds_list <- list("brcds_mRNA" = brcds_mRNA, "brcds_miRNA" = brcds_miRNA, "brcds_DNAme" = brcds_DNAme, "brcds_SNV" = brcds_SNV)
GeoMeans_list <- list("GeoMeans_mRNA" = expdat_mRNA_geoMeans, "GeoMeans_miRNA" = expdat_miRNA_geoMeans)

saveMetadata(OutDir = OutDir, Seed = Seed, smpls = smpls, expdat_list = expdat_list, brcds_list = brcds_list, Prjcts = Prjcts, GeoMeans_list = GeoMeans_list)

## -------------------------------------

## ------------------------------------------------------------------------------

## ------------------------------------------------------------------------------
## ------------------------------------------------ 02 REFERENCE AND TARGET DATASETS FROM ONE PROJECT ----

## This section has 2 parts.
## Firstly a dataset is created for the reference dataset
## The reference dataset is the full target dataset for a selected project.
## This means it contains multi-omics data for all samples.
## We use factorizations of the reference dataset to derive groundtruth factors 
## Then target datasets are created. 
## Each target dataset is a multi-omics dataset for a subset of the samples used to create the corresponding reference dataset.

print(Sys.time())
print("Reference set - one project")

## ---------------- SET ENVIRONMENT ----

## SET PARAMETERS
Prjct = "TCGA-LAML" # the TCGA project to use for creating reference and target datasets
InputDir = "DownloadedData_unfltrd" # location of downloaded TCGA data
GeoMeans = "Trg"
TopD = 5000 ## how many features to retain

## SET SEED
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

## -------------------------------------

## ----------------- SELECT SAMPLES ----

print(Sys.time())
print("Select samples")

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
BarcodesDir <- paste0('Trg_',substr(Prjct,6,nchar(Prjct)), "_Full_barcodes")
if(!dir.exists(BarcodesDir)){
  dir.create(BarcodesDir, showWarnings = FALSE, recursive = TRUE)
}

# import rds files
expdat_mRNA = readRDS(file.path(InputDir,paste0(Prjct,"_mRNA.rds")))
expdat_miRNA = readRDS(file.path(InputDir,paste0(Prjct,"_miRNA.rds")))
expdat_DNAme = readRDS(file.path(InputDir,paste0(Prjct,"_DNAme.rds")))

# tidy up and filter omics data (parsing file)
expdat_mRNA_fltr = mRNAFileParsing(expdat_mRNA)
expdat_miRNA_fltr = miRNAFileParsing(expdat_miRNA)
expdat_DNAme_fltr = DNAmeFileParsing(expdat_DNAme)

# samples with shared omics
omics_list <- list("mRNA" = expdat_mRNA_fltr, "miRNA" = expdat_miRNA_fltr, "DNAme" = expdat_DNAme_fltr)
smpls_shared <- extractSampleNamesShared(omics_list)

# for each submittor select the min barcode from those related to samples with shared omics
brcds_list = sapply(X = omics_list, FUN = function(expdat, smpls_shared, Prjct){
  brcds = selectMinBarcodeSampleNames(expdat = expdat, smpls_shared = smpls_shared, Prjct = Prjct)
  return(brcds)
}, smpls_shared, Prjct, simplify = FALSE)

## Merge barcodes
brcds_mRNA <- brcds_list$mRNA
brcds_miRNA <- brcds_list$miRNA
brcds_DNAme <- brcds_list$DNAme

## Save barcodes
saveRDS(brcds_mRNA,file.path(BarcodesDir,'brcds_mRNA.rds'))
saveRDS(brcds_miRNA,file.path(BarcodesDir,'brcds_miRNA.rds'))
saveRDS(brcds_DNAme,file.path(BarcodesDir,'brcds_DNAme.rds'))

## -------------------------------------

## --------- CREATION OF REFERENCE DATASET (THE FULL TARGET DATASET BEFORE SUBSETTING) ----

print(Sys.time())
print("Reference set creation")

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
OutDir = paste0('Trg_',substr(Prjct,6,nchar(Prjct)),'_Full_', TopD,"D")
if(!dir.exists(OutDir)){
  dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
}

## IMPORT BARCODES FROM PREVIOUS STEP
brcds_mRNA = readRDS(file.path(BarcodesDir,'brcds_mRNA.rds'))
brcds_miRNA = readRDS(file.path(BarcodesDir,'brcds_miRNA.rds'))
brcds_DNAme = readRDS(file.path(BarcodesDir,'brcds_DNAme.rds'))

# save samples for shuffling and keeping track of names later on
smpls = substr(brcds_mRNA$brcds,1,16)
smpls = sample(smpls, length(smpls), replace = FALSE)

## ---------- mRNA data preprocessing

print(Sys.time())
print("mRNA")

## Read and merge data files
expdat_mRNA_merged <- mRNAdataPreparation(Prjcts = Prjct, brcds_mRNA = brcds_mRNA, InputDir = InputDir)

## Filtering
expdat_mRNA_fltr <- mRNAFiltering(expdat_mRNA = expdat_mRNA_merged)

## Normalization
expdat_mRNA_norm <- countsNormalization(expdat = expdat_mRNA_fltr, GeoMeans = GeoMeans)

## Transformation
expdat_mRNA_log <- countsTransformation(expdat_count = expdat_mRNA_norm$counts, TopD = TopD)

## Save data
expdat_mRNA_sorted <- orderAndSaveData(expdat = expdat_mRNA_log, smpls = smpls, OutDir = OutDir, fileName = "mRNA.csv")

## PCA to get estimates of necessary factor count
ElbowK_mRNA = PCAtools::pca(mat = expdat_mRNA_sorted[rowSums(is.na(expdat_mRNA_sorted))==0,], scale = TRUE)

## get elbow point and % of variance explained
PCVarPrcnt_mRNA = ElbowK_mRNA$variance/sum(ElbowK_mRNA$variance)
ElbowK_mRNA = PCAtools::findElbowPoint(ElbowK_mRNA$variance)

## CLEAN SPACE
rm(expdat_mRNA_merged, expdat_mRNA_fltr, expdat_mRNA_norm, expdat_mRNA_log)
gc()

## ---------- miRNA data preprocessing

print(Sys.time())
print("miRNA")

## Read and merge data files
expdat_miRNA_merged <- miRNAdataPreparation(Prjcts = Prjct, brcds_miRNA = brcds_miRNA, InputDir = InputDir)

## Filtering
expdat_miRNA_fltr <- miRNAFiltering(expdat_miRNA = expdat_miRNA_merged)

## Normalization
expdat_miRNA_norm <- countsNormalization(expdat = expdat_miRNA_fltr, GeoMeans = GeoMeans)

## Transformation
expdat_miRNA_log <- countsTransformation(expdat_count = expdat_miRNA_norm$counts, TopD = TopD)

## Save data
expdat_miRNA_sorted <- orderAndSaveData(expdat = expdat_miRNA_log, smpls = smpls, OutDir = OutDir, fileName = "miRNA.csv")

## PCA to get estimates of necessary factor count
ElbowK_miRNA = PCAtools::pca(mat = expdat_miRNA_sorted[rowSums(is.na(expdat_miRNA_sorted))==0,], scale = TRUE)

## % of variance explained and elbow point
PCVarPrcnt_miRNA = ElbowK_miRNA$variance/sum(ElbowK_miRNA$variance)
ElbowK_miRNA = PCAtools::findElbowPoint(ElbowK_miRNA$variance)

## CLEAN SPACE
rm(expdat_miRNA_merged, expdat_miRNA_fltr, expdat_miRNA_norm, expdat_miRNA_log)
gc()

## ---------- DNAme data preprocessing

print(Sys.time())
print("DNAme")

## Read and merge data files
expdat_DNAme_merged <- DNAmedataPreparation(Prjcts = Prjct, brcds_DNAme = brcds_DNAme, InputDir = InputDir)

## Filtering
expdat_DNAme_fltr <- DNAmeFiltering(expdat_DNAme = expdat_DNAme_merged, TopD = TopD)

## Save data
expdat_DNAme_sorted <- orderAndSaveData(expdat = expdat_DNAme_fltr, smpls = smpls, OutDir = OutDir, fileName = "DNAme.csv")

## PCA to get estimates of necessary factor count
ElbowK_DNAme = PCAtools::pca(mat = expdat_DNAme_sorted[rowSums(is.na(expdat_DNAme_sorted))==0,], scale = TRUE)

## % of variance explained and elbow point
PCVarPrcnt_DNAme = ElbowK_DNAme$variance/sum(ElbowK_DNAme$variance)
ElbowK_DNAme = PCAtools::findElbowPoint(ElbowK_DNAme$variance)

## CLEAN SPACE
rm(expdat_DNAme_merged, expdat_DNAme_fltr)
gc()

## -------------------------------------

## ------------------ SAVE METADATA ----

## SAVE METADATA
expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted)
brcds_list <- list("brcds_mRNA" = brcds_mRNA, "brcds_miRNA" = brcds_miRNA, "brcds_DNAme" = brcds_DNAme)
PCA_list <- list("PCVarPrcnt_mRNA" = PCVarPrcnt_mRNA, "PCVarPrcnt_miRNA" = PCVarPrcnt_miRNA, "PCVarPrcnt_DNAme" = PCVarPrcnt_DNAme)
Elbow_list <- list("ElbowK_mRNA" = ElbowK_mRNA, "ElbowK_miRNA" = ElbowK_miRNA, "ElbowK_DNAme" = ElbowK_DNAme)

saveMetadata(OutDir = OutDir, Seed = Seed, smpls = smpls, expdat_list = expdat_list, brcds_list = brcds_list, PCA_list = PCA_list, Elbow_list = Elbow_list, Prjcts = Prjct)

## -------------------------------------

## ---------------- TARGET DATASETS (SUBSETS OF THE REFERENCE DATASET) ----

print(Sys.time())
print("Target datasets")

## PARAMETERS FOR IMPORTING DATA
OutDir = paste0('Trg_',substr(Prjct,6,nchar(Prjct)),'_Full_', TopD,"D") # output directory used for reference (full target) dataset
expdat_meta = readRDS(file.path(OutDir,'expdat_meta.rds'))
Prjcts = expdat_meta$Prjcts
TopD = expdat_meta$TopD

## SET PARAMETER
SS_size = 5 ## how many samples per project in each sub set

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
OutDir_SS = paste0('Trg_',paste0(substr(Prjcts,6,nchar(Prjcts)),collapse='_'),'_SS',SS_size,'_',TopD,'D')
if(!dir.exists(OutDir_SS)){
  dir.create(OutDir_SS, showWarnings = FALSE, recursive = TRUE)
}

## ----------------- SELECT SAMPLES ---
smpls_raw = expdat_meta$smpls
SS_count = floor(length(smpls_raw)/SS_size)
smpls_SS = vector("list")

## FOR EACH SUBSET, SELECT SAMPLES RANDOMLY
for (i in 1:SS_count){
  smpls_SS[[i]] = sample(smpls_raw, SS_size, replace=FALSE)
  names(smpls_SS)[i] <- paste("SS", i, sep = "_")
  smpls_raw = smpls_raw[!is.element(smpls_raw, smpls_SS[[i]])]
}

## SELECT CORRESPONDING BARCODES
brcds_mRNA_SS <- lapply(smpls_SS, selectBarcodeSampleNames, expdat_meta$brcds_mRNA)
brcds_miRNA_SS <- lapply(smpls_SS, selectBarcodeSampleNames, expdat_meta$brcds_miRNA)
brcds_DNAme_SS <- lapply(smpls_SS, selectBarcodeSampleNames, expdat_meta$brcds_DNAme)

## SAVE BARCODES FOR EACH SUBSET
brcds_SS = list(
  Seed = Seed,
  SS_size = SS_size,
  SS_count = SS_count,
  smpls_SS = smpls_SS,
  brcds_mRNA_SS = brcds_mRNA_SS,
  brcds_miRNA_SS = brcds_miRNA_SS,
  brcds_DNAme_SS = brcds_DNAme_SS
)
saveRDS(brcds_SS, file.path(OutDir_SS,'brcds_SS.rds'))

## ------- CREATION OF SUBSETS DATA ---

## IMPORT METADATA FROM PREVIOUS STEP 
brcds_SS = readRDS(file.path(OutDir_SS,'brcds_SS.rds'))

InputDir = file.path("DownloadedData_unfltrd")

tmp <- lapply(1:brcds_SS$SS_count, createSubsetOneProject, brcds_SS = brcds_SS, OutDir = OutDir_SS, GeoMeans = GeoMeans, InputDir = InputDir)

## ------------------ SAVE METADATA ---
expdat_meta = list(
  Seed = Seed,
  if_SNV = FALSE,
  SS_count = brcds_SS$SS_count
)

saveRDS(expdat_meta, file.path(OutDir_SS,'expdat_meta_SS.rds'))

expdat_meta.json = rjson::toJSON(expdat_meta)
write(expdat_meta.json,file.path(OutDir_SS,'expdat_meta_SS.json'))
## -------------------------------------

## ------------------------------------------------------------------------------

## ------------------------------------------------------------------------------
## ---------------------------------------------- 03 REFERENCE AND TARGET DATASETS FROM MULTIPLE PROJECTS ----

## This section has 2 parts.
## Firstly a dataset is created for the reference dataset
## The reference dataset is the full target dataset for the selected projects.
## This means it contains multi-omics data for all samples.
## We use factorizations of the reference dataset to derive groundtruth factors 
## Then target datasets are created. 
## Each target dataset is a multi-omics dataset for a subset of the samples used to create the corresponding reference dataset.

print(Sys.time())
print("Target set - multi project")

## ---------------- SET ENVIRONMENT ----

## SET PARAMETERS
Prjcts = c("TCGA-LAML", "TCGA-SKCM")
InputDir = "DownloadedData_unfltrd"
TopD = 5000 ## how many features to retain
GeoMeans = "Trg"

## SET SEED
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

## -------------------------------------

## ----------------- SELECT SAMPLES ----

print(Sys.time())
print("Select samples")

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
BarcodesDir <- paste0('Trg_', paste0(substr(Prjcts,6,nchar(Prjcts)),collapse='_'), '_Full_barcodes')
if(!dir.exists(BarcodesDir)){
  dir.create(BarcodesDir, showWarnings = FALSE, recursive = TRUE)
}

## FOR EACH PROJECT
brcds_list <- sapply(X = Prjcts, FUN = function(Prjct){
  
  print("[*]-------------------------------------")
  print(paste("Project", Prjct))
  
  # import rds files
  print("Read omics files")
  expdat_mRNA = readRDS(file.path(InputDir,paste0(Prjct,"_mRNA.rds")))
  expdat_miRNA = readRDS(file.path(InputDir,paste0(Prjct,"_miRNA.rds")))
  expdat_DNAme = readRDS(file.path(InputDir,paste0(Prjct,"_DNAme.rds")))
  
  # tidy up and filter omics data (parsing file)
  print("Parse omics data")
  expdat_mRNA_fltr = mRNAFileParsing(expdat_mRNA)
  expdat_miRNA_fltr = miRNAFileParsing(expdat_miRNA)
  expdat_DNAme_fltr = DNAmeFileParsing(expdat_DNAme)
  
  # samples with shared omics
  print("Select shared samples between omics")
  omics_list <- list("mRNA" = expdat_mRNA_fltr, "miRNA" = expdat_miRNA_fltr, "DNAme" = expdat_DNAme_fltr)
  smpls_shared <- extractSampleNamesShared(omics_list)
  
  # for each submittor select the min barcode from those related to samples with shared omics
  print("Select the min barcode for each submittor")
  brcds_list_tmp = sapply(X = omics_list, FUN = function(expdat, smpls_shared, Prjct){
    brcds = selectMinBarcodeSampleNames(expdat = expdat, smpls_shared = smpls_shared, Prjct = Prjct)
    return(brcds)
  }, smpls_shared, Prjct, simplify = FALSE)
  
  print("[*]-------------------------------------")
  
  return(brcds_list_tmp)
})

## Merge barcodes
brcds_mRNA <- do.call(rbind, brcds_list["mRNA",])
brcds_miRNA <- do.call(rbind, brcds_list["miRNA",])
brcds_DNAme <- do.call(rbind, brcds_list["DNAme",])

## Save barcodes
saveRDS(brcds_mRNA,file.path(BarcodesDir,'brcds_mRNA.rds'))
saveRDS(brcds_miRNA,file.path(BarcodesDir,'brcds_miRNA.rds'))
saveRDS(brcds_DNAme,file.path(BarcodesDir,'brcds_DNAme.rds'))

## -------------------------------------

## --------- CREATION OF REFERENCE DATASET (THE FULL TARGET DATASET BEFORE SUBSETTING) ----

print(Sys.time())
print("Reference set creation")

## CREATE A DIRECTORY TO SAVE THE FINAL FILES
OutDir = paste0('Trg_', paste0(substr(Prjcts,6,nchar(Prjcts)),collapse='_'), '_Full_',TopD,"D")
if(!dir.exists(OutDir)){
  dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
}

## IMPORT BARCODES FROM PREVIOUS STEP
brcds_mRNA = readRDS(file.path(BarcodesDir,'brcds_mRNA.rds'))
brcds_miRNA = readRDS(file.path(BarcodesDir,'brcds_miRNA.rds'))
brcds_DNAme = readRDS(file.path(BarcodesDir,'brcds_DNAme.rds'))

# save samples for shuffling and keeping track of names later on
smpls = substr(brcds_mRNA$brcds,1,16)
smpls = sample(smpls, length(smpls), replace = FALSE)

## ---------- mRNA data preprocessing

print(Sys.time())
print("mRNA")

## Read and merge files
expdat_mRNA_merged <- mRNAdataPreparation(Prjcts = Prjcts, brcds_mRNA = brcds_mRNA, InputDir = InputDir)

## Filtering
expdat_mRNA_fltr <- mRNAFiltering(expdat_mRNA = expdat_mRNA_merged)

## Normalization
expdat_mRNA_norm <- countsNormalization(expdat = expdat_mRNA_fltr, GeoMeans = GeoMeans)

## Transformation
expdat_mRNA_log <- countsTransformation(expdat_count = expdat_mRNA_norm$counts, TopD = TopD)

## Save data
expdat_mRNA_sorted <- orderAndSaveData(expdat = expdat_mRNA_log, smpls = smpls, OutDir = OutDir, fileName = "mRNA.csv")

## PCA to get estimates of necessary factor count
ElbowK_mRNA = PCAtools::pca(mat = expdat_mRNA_sorted[rowSums(is.na(expdat_mRNA_sorted))==0,], scale = TRUE)

## get elbow point and % of variance explained
PCVarPrcnt_mRNA = ElbowK_mRNA$variance/sum(ElbowK_mRNA$variance)
ElbowK_mRNA = PCAtools::findElbowPoint(ElbowK_mRNA$variance)

## CLEAN SPACE
rm(expdat_mRNA_merged, expdat_mRNA_fltr, expdat_mRNA_norm, expdat_mRNA_log)
gc()

## ---------- miRNA data preprocessing 

print(Sys.time())
print("miRNA")

## Read and merge files
expdat_miRNA_merged <- miRNAdataPreparation(Prjcts = Prjcts, brcds_miRNA = brcds_miRNA, InputDir = InputDir)

## Filtering
expdat_miRNA_fltr <- miRNAFiltering(expdat_miRNA = expdat_miRNA_merged)

## Normalization
expdat_miRNA_norm <- countsNormalization(expdat = expdat_miRNA_fltr, GeoMeans = GeoMeans)

## Transformation
expdat_miRNA_log <- countsTransformation(expdat_count = expdat_miRNA_norm$counts, TopD = TopD)

## Save data
expdat_miRNA_sorted <- orderAndSaveData(expdat = expdat_miRNA_log, smpls = smpls, OutDir = OutDir, fileName = "miRNA.csv")

## PCA to get estimates of necessary factor count
ElbowK_miRNA = PCAtools::pca(mat = expdat_miRNA_sorted[rowSums(is.na(expdat_miRNA_sorted))==0,], scale = TRUE)

## % of variance explained and elbow point
PCVarPrcnt_miRNA = ElbowK_miRNA$variance/sum(ElbowK_miRNA$variance)
ElbowK_miRNA = PCAtools::findElbowPoint(ElbowK_miRNA$variance)

## CLEAN SPACE
rm(expdat_miRNA_merged, expdat_miRNA_fltr, expdat_miRNA_norm, expdat_miRNA_log)
gc()

## ---------- DNAme data preprocessing

print(Sys.time())
print("DNAme")

## Read and merge files
expdat_DNAme_merged <- DNAmedataPreparation(Prjcts = Prjcts, brcds_DNAme = brcds_DNAme, InputDir = InputDir)

## Filtering
expdat_DNAme_fltr <- DNAmeFiltering(expdat_DNAme = expdat_DNAme_merged, TopD = TopD)

## Save data
expdat_DNAme_sorted <- orderAndSaveData(expdat = expdat_DNAme_fltr, smpls = smpls, OutDir = OutDir, fileName = "DNAme.csv")

## PCA to get estimates of necessary factor count
ElbowK_DNAme = PCAtools::pca(mat = expdat_DNAme_sorted[rowSums(is.na(expdat_DNAme_sorted))==0,], scale = TRUE)

## % of variance explained and elbow point
PCVarPrcnt_DNAme = ElbowK_DNAme$variance/sum(ElbowK_DNAme$variance)
ElbowK_DNAme = PCAtools::findElbowPoint(ElbowK_DNAme$variance)

## CLEAN SPACE
rm(expdat_DNAme_merged, expdat_DNAme_fltr)
gc()

## -------------------------------------

## ------------------ SAVE METADATA ----

expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted)
brcds_list <- list("brcds_mRNA" = brcds_mRNA, "brcds_miRNA" = brcds_miRNA, "brcds_DNAme" = brcds_DNAme)
PCA_list <- list("PCVarPrcnt_mRNA" = PCVarPrcnt_mRNA, "PCVarPrcnt_miRNA" = PCVarPrcnt_miRNA, "PCVarPrcnt_DNAme" = PCVarPrcnt_DNAme)
Elbow_list <- list("ElbowK_mRNA" = ElbowK_mRNA, "ElbowK_miRNA" = ElbowK_miRNA, "ElbowK_DNAme" = ElbowK_DNAme)

saveMetadata(OutDir = OutDir, Seed = Seed, smpls = smpls, expdat_list = expdat_list, brcds_list = brcds_list, PCA_list = PCA_list, Elbow_list = Elbow_list, Prjcts = Prjcts)

## -------------------------------------

## ---------------- TARGET DATASETS (SUBSETS OF THE REFERENCE DATASET) ----

print(Sys.time())
print("Target datasets")

SS_size <- 5
# output directory used for reference (full target) dataset
OutDir = paste0('Trg_', paste0(substr(Prjcts,6,nchar(Prjcts)),collapse='_'), '_Full_',TopD,"D")

## PARAMETERS IMPORTATION
## CREATE A DIRECTORY TO SAVE THE FINAL FILES
SSOutDir = paste0('Trg_', paste0(substr(Prjcts,6,nchar(Prjcts)),collapse='_'), '_SS',SS_size, '_', TopD,"D")
if(!dir.exists(SSOutDir)){
  dir.create(SSOutDir, showWarnings = FALSE, recursive = TRUE)
}

## PARAMETERS IMPORTATION
expdat_meta = readRDS(file.path(OutDir,'expdat_meta.rds'))

## ----------------- SELECT SAMPLES ---

Projects = unique(expdat_meta$brcds_mRNA$prjct)

smpls_raw_list <- sapply(Projects, function(prjct, brcds){
  smpls <- brcds[brcds$prjct==prjct,"brcds"]
  smpls <- substr(smpls,1,16)
  return(smpls)
},expdat_meta$brcds_mRNA)

SS_count <- sapply(smpls_raw_list, function(smpls_prjct, SS_size){
  SS_count = floor(length(smpls_prjct)/SS_size)
  return(SS_count)
}, SS_size)

SS_count_min <- min(SS_count)

smpls_SS = vector("list")

## FOR EACH SUBSET, SELECT SAMPLES RANDOMLY
for(i in 1:SS_count_min){
  ## RANDOM SELECTION OF SAMPLE NAMES
  smpls_SS_list <- lapply(smpls_raw_list, function(smpls_prjct, SS_size){
    smpls_SS_tmp <- sample(smpls_prjct, SS_size, replace=FALSE)
    return(smpls_SS_tmp)
  }, SS_size)
  ## CREATE ONE LIST
  smpls_SS[[i]] <- unlist(smpls_SS_list, use.names = FALSE)
  names(smpls_SS)[i] <- paste("SS", i, sep = "_")
  ## REMOVE SELECTED SAMPLES
  smpls_raw_list <- sapply(Projects, function(prjct, smpls_raw_list, smpls_SS_list){
    smpls <- smpls_raw_list[[prjct]]
    smpls_raw_tmp <- smpls[!is.element(smpls, smpls_SS_list[[prjct]])]
    return(smpls_raw_tmp)
  }, smpls_raw_list, smpls_SS_list)
}

## SELECT CORRESPONDING BARCODES
brcds_mRNA_SS <- lapply(smpls_SS, selectBarcodeSampleNames, expdat_meta$brcds_mRNA)
brcds_miRNA_SS <- lapply(smpls_SS, selectBarcodeSampleNames, expdat_meta$brcds_miRNA)
brcds_DNAme_SS <- lapply(smpls_SS, selectBarcodeSampleNames, expdat_meta$brcds_DNAme)

## SAVE BARCODES FOR EACH SUBSET
brcds_SS = list(
  Seed = Seed,
  SS_size = SS_size,
  SS_count = SS_count_min,
  smpls_SS = smpls_SS,
  brcds_mRNA_SS = brcds_mRNA_SS,
  brcds_miRNA_SS = brcds_miRNA_SS,
  brcds_DNAme_SS = brcds_DNAme_SS
)

saveRDS(brcds_SS, file.path(SSOutDir,'brcds_SS.rds'))

gc()

## ------- CREATION OF SUBSETS DATA ---

## IMPORT METADATA FROM PREVIOUS STEP
brcds_SS = readRDS(file.path(SSOutDir,'brcds_SS.rds'))

tmp <- lapply(1:brcds_SS$SS_count, createSubsetMultiProjects, brcds_SS = brcds_SS, OutDir = SSOutDir, GeoMeans = GeoMeans, InputDir = InputDir)

## ------------------ SAVE METADATA ---
expdat_meta = list(
  Seed = Seed,
  if_SNV = FALSE,
  SS_count = brcds_SS$SS_count
)

saveRDS(expdat_meta, file.path(SSOutDir,'expdat_meta_SS.rds'))

expdat_meta.json = rjson::toJSON(expdat_meta)
write(expdat_meta.json,file.path(SSOutDir,'expdat_meta_SS.json'))
## -------------------------------------

## ------------------------------------------------------------------------------
