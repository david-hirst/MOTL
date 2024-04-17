## MT - 20231108

## CALL FUNCTION TO DOWNLOAD TCGA DATA FROM GDC PORTAL

## Running this script downloads mRNA, miRNA, DNAme and SNV omics data for TCGA samples.
## It does this for each sample, from each TCGA project, using the TCGAbiolinks R package.
## For each TCGA project and omics, files are downloaded separately for each sample.
## The data is consolidated and an .rds is saved for each project and omics.
## Each .rds file is saved in the directory specified below.

## LIBRARIES
library(TCGAbiolinks)

## FUNCTIONS
source("TCGA_preprocessedData_functions.R")

## SET PARAMETERS
SaveDirectory = "DownloadedData_unfltrd"

## MAIN
## DOWNLOAD 4 OMICS DATA FILES FROM ALL TCGA PROJECT
downloadAllTCGAProjectData(SaveDirectory = SaveDirectory)
