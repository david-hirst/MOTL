## downloaded from github 20240222

## LIBRARIES
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)
library(dplyr)

## FUNCTIONS
source('TCGA_preprocessedData_functions.R')

## ------------------------------------------------------------------------------
## ------------------------------------------------------------ LEARNING SET ----

print(Sys.time())
print("Learning set")

## ---------------- SET ENVIRONMENT ----

## SET PARAMETERS 
PrjctExcl = c("TCGA-GBM") # projects to exclude
InputDir = file.path('..','..','DownloadedData_unfltrd')
BarcodesDir = 'Lrn_barcodes'
GeoMeans = "Lrn"
TopD = 5000 ## how many features to retain

## SET SEED
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

## EXCLUDE ANY CANCER TYPES
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
