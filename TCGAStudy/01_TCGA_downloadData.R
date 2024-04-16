## MT - 20231108

## FUNCTIONS TO DOWNLOAD TCGA DATA FROM GDC PORTAL

## LIBRARIES
library(TCGAbiolinks)

## FUNCTIONS
source("TCGA_preprocessedData_functions.R")

## SET PARAMETERS
SaveDirectory = "DownloadedData_unfltrd"

## MAIN
## DOWNLOAD 4 OMICS DATA FILES FROM ALL TCGA PROJECT
downloadAllTCGAProjectData(SaveDirectory = SaveDirectory)
