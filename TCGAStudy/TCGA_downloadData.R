## MT - 20231108

## FUNCTIONS TO DOWNLOAD TCGA DATA FROM GDC PORTAL

## SET ENVIRONMENT
setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/MOTL/TCGAStudy/TCGA_analysisTests/")

## LIBRARIES
library(TCGAbiolinks)

## FUNCTIONS
downloadTCGAProjectData <- function(Prjct, SaveDirectory){
  #' Download four omics data from one TCGA project 
  #'
  #' For a TCGA project name, this function:
  #' 1. Create queries for mRNA, miRNA, DNAme and SNV data
  #' 2. Download the corresponding raw data files
  #' 3. Read these raw data files
  #' 4. Concatenate files and save them into .rds files
  #' 5. Remove the raw data files
  #' 
  #' @param Prjct TCGA project name
  #' @param SaveDirectory Directory name to save omics data files
  
  ## 
  print(paste("Download data of", Prjct))
  
  ## Set parameters
  fpc = 20 # files per chunk setting
  
  ## set up download queries
  
  query_mRNA = GDCquery(
    project = Prjct, 
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification"
  )
  
  query_miRNA = GDCquery(
    project = Prjct,  
    data.category = "Transcriptome Profiling", 
    data.type = "miRNA Expression Quantification"
  )
  
  query_DNAme = GDCquery(
    project = Prjct,
    data.category = "DNA Methylation", 
    data.type = "Methylation Beta Value",
    platform = "Illumina Human Methylation 450"
  )
  
  query_SNV = GDCquery(
    project = Prjct, 
    data.category = "Simple Nucleotide Variation", 
    access = "open", 
    data.type = "Masked Somatic Mutation"
  )
  
  ## result summaries of the queries
  
  query_mRNA_res = getResults(query_mRNA)
  query_miRNA_res = getResults(query_miRNA)
  query_DNAme_res = getResults(query_DNAme)
  query_SNV_res = getResults(query_SNV)
  
  ## files counts for mRNA and DNAme
  
  files_mRNA = length(query_mRNA_res$file_name)
  files_miRNA = length(query_miRNA_res$file_name)
  files_DNAme = length(query_DNAme_res$file_name)
  files_SNV = length(query_SNV_res$file_name)
  
  ## files per chunk check - can't have any chunks of just 1 file
  
  fpc_mRNA = ((files_mRNA - floor(files_mRNA/fpc)*fpc)==1)+0
  fpc_miRNA = ((files_miRNA - floor(files_miRNA/fpc)*fpc)==1)+0
  fpc_DNAme = ((files_DNAme - floor(files_DNAme/fpc)*fpc)==1)+0
  fpc_SNV = ((files_SNV - floor(files_SNV/fpc)*fpc)==1)+0
  
  ## download the files
  
  GDCdownload(query_mRNA, files.per.chunk = fpc+fpc_mRNA)
  GDCdownload(query_miRNA, files.per.chunk = fpc+fpc_miRNA)
  GDCdownload(query_DNAme, files.per.chunk = fpc+fpc_DNAme)
  GDCdownload(query_SNV, files.per.chunk = fpc+fpc_SNV)
  
  ## load the files into R and aggregate for each omics
  
  expdat_mRNA = GDCprepare(query = query_mRNA)
  expdat_miRNA = GDCprepare(query = query_miRNA)
  expdat_DNAme = GDCprepare(query = query_DNAme)
  expdat_SNV = GDCprepare(query = query_SNV)
  
  ## save these as rds files
  
  saveRDS(expdat_mRNA, file.path(SaveDirectory,paste0(Prjct,"_mRNA.rds")))
  saveRDS(expdat_miRNA, file.path(SaveDirectory,paste0(Prjct,"_miRNA.rds")))
  saveRDS(expdat_DNAme, file.path(SaveDirectory,paste0(Prjct,"_DNAme.rds")))
  saveRDS(expdat_SNV, file.path(SaveDirectory,paste0(Prjct,"_SNV.rds")))
  
  ## delete the downloaded raw data files
  
  unlink(file.path("GDCdata",Prjct),recursive=TRUE)
}

downloadAllTCGAProjectData <- function(SaveDirectory){
  #' Download four omics data from all TCGA projects
  #'
  #' @param SaveDirectory Directory name to save data files
  
  ## get the vector of TCGA projects
  ProjectsAll = TCGAbiolinks:::getGDCprojects()$project_id
  ProjectsTCGA = ProjectsAll[substr(ProjectsAll,1,5)=="TCGA-"]
  
  ## create a directory to save the final files
  if(!dir.exists(SaveDirectory)){
    dir.create(SaveDirectory, showWarnings = FALSE, recursive = TRUE)
  }
  
  ## MT: add condition
  ## check which projects have files
  SavedFiles = list.files(SaveDirectory)
  if(length(SavedFiles) != 0){
    SavedFiles = SavedFiles[substr(SavedFiles,1,5)=="TCGA-"]
    SavedFiles = strsplit(SavedFiles,"_")
    SavedFiles = as.data.frame(do.call("rbind",SavedFiles))
    ## only get files for projects that haven't been downloaded yet
    SavedFiles = aggregate(V2~V1,FUN=length,data=SavedFiles)
    SavedFiles = SavedFiles$V1[SavedFiles$V2==4]
    ProjectsTCGA = ProjectsTCGA[!is.element(ProjectsTCGA,SavedFiles)]
  }
  
  ## MT: change loop with apply (better, faster) + create function downloadTCGAProjectData(Prjct, SaveDirectory)
  lapply(X = ProjectsTCGA, FUN = function(Prjct, SaveDirectory){
    downloadTCGAProjectData(Prjct = Prjct, SaveDirectory = SaveDirectory)
  }, SaveDirectory)
}

## SET PARAMETERS
SaveDirectory = "DownloadedData_unfltrd"

## MAIN
## DOWNLOAD 4 OMICS DATA FILES FROM ALL TCGA PROJECT
downloadAllTCGAProjectData(SaveDirectory = SaveDirectory)
