##############
#### preprocess target data before direct factorization with mofa
##############

library(SummarizedExperiment)
library(sesame)
library(dplyr)
library(rjson)
library(biomaRt)

source('TCGA_preprocessedData_functions.R')

## SET SEED
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

##
### parameters and directories
##

TopD = 5000

InDir = file.path('RawData')
# HGCCgroups = c("CL", "MS", "PN")
HGCCgroups = c("Normal", "PN")
OutDir = file.path(paste0('Trg_',paste(HGCCgroups,collapse='_'),'_',TopD,"D"),"MOFA")
if(!dir.exists(OutDir)){
  dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
}

##
### sample metadata
##

smpls = read.table(
  file.path(InDir,"phenodata_RNAseq.tsv"),
  header = TRUE
)

# tidy filter and shuffle sample names
smpls$Sample = sub("[.]","-",smpls$Sample)
smpls = smpls$Sample[is.element(smpls$HGCC,HGCCgroups)]
smpls = sample(smpls, length(smpls), replace = FALSE)

#####
# mRNA
####

### import

expdat_mRNA = read.table(
  file.path(InDir,"1_genes_vst_filtered_DESeq2_normalized_HGCC.tsv"),
  header = TRUE
)
expdat_mRNA = as.matrix(expdat_mRNA)

### tidy, filter and order based on samples
colnames(expdat_mRNA) = sub("[.]","-",colnames(expdat_mRNA))
expdat_mRNA = expdat_mRNA[,smpls]

### filter features

expdat_mRNA = expdat_mRNA[!regexpr('PAR_Y',rownames(expdat_mRNA), fixed = TRUE)>0,]
rownames(expdat_mRNA) = substr(rownames(expdat_mRNA),1,regexpr('[.]',rownames(expdat_mRNA))-1)

biomart_mRNA = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_mRNA <- getBM(attributes = c("chromosome_name", "ensembl_gene_id"), 
                      filters = "ensembl_gene_id", 
                      values = rownames(expdat_mRNA), bmHeader = T, mart = biomart_mRNA)

FtrsKeep = (rowMeans(is.na(expdat_mRNA))<0.2)  &
  (rowMeans(expdat_mRNA==0, na.rm = TRUE)<0.9) &
  (rowVars(expdat_mRNA, na.rm = TRUE)>0) &
  (!is.element(rownames(expdat_mRNA),
               biomart_mRNA$`Gene stable ID`[biomart_mRNA$`Chromosome/scaffold name`=='Y']
  )) &
  (is.element(rownames(expdat_mRNA), biomart_mRNA$`Gene stable ID`))

FtrsKeep[is.na(FtrsKeep)]=FALSE
expdat_mRNA = expdat_mRNA[FtrsKeep,]

## topD features

FtrsKeep = base::rank(-rowVars(expdat_mRNA, na.rm = TRUE), ties.method = "first") <= TopD
expdat_mRNA = expdat_mRNA[FtrsKeep,]

## save csv

write.table(expdat_mRNA, file.path(OutDir, 'mRNA.csv'), 
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = FALSE)

#####
# DNAme
####

### import and put into correct format

## read each line in and split
expdat_DNAme = strsplit(
  readLines(file.path(InDir,"5mC_bVals_HGCC.txt"), n = -1), 
  split = "\t"
)

# keep sample names and tidy
DNAme_smpls = expdat_DNAme[[1]]
DNAme_smpls = sub("^2","X2",DNAme_smpls)

## only keep lines with specified number of elements
correct_length = 15
expdat_DNAme_chk = lapply(expdat_DNAme, function(x) length(x) == correct_length)
expdat_DNAme = expdat_DNAme[unlist(expdat_DNAme_chk)]
expdat_DNAme = do.call(rbind, expdat_DNAme)

# make sure only rows with probe ids are kept
expdat_DNAme = expdat_DNAme[is.element(substr(expdat_DNAme[,1],1,2),c('cg','ch')),]

# move ids to row names
rownames(expdat_DNAme) = expdat_DNAme[,1]
expdat_DNAme = expdat_DNAme[,-1]
colnames(expdat_DNAme) = DNAme_smpls

# convert to numeric
expdat_DNAme = matrix(
  as.numeric(expdat_DNAme),
  ncol = ncol(expdat_DNAme),
  dimnames = list(
    rownames(expdat_DNAme),
    colnames(expdat_DNAme)
  )
)

# convert to M-Values
expdat_DNAme = sesame::BetaValueToMValue(expdat_DNAme)

# filter and order based on samples
expdat_DNAme = expdat_DNAme[,smpls]

# filter features - they already did typical (including sex chromosomes based) DNAme filtering 
FtrsKeep = (rowMeans(is.na(expdat_DNAme))<0.2) &
  (rowVars(expdat_DNAme, na.rm = TRUE)>0)
FtrsKeep[is.na(FtrsKeep)]=FALSE
expdat_DNAme = expdat_DNAme[FtrsKeep,]

FtrsKeep = base::rank(-rowVars(expdat_DNAme, na.rm = TRUE), ties.method = "first") <= TopD
expdat_DNAme = expdat_DNAme[FtrsKeep,] 

## save csv

write.table(expdat_DNAme, file.path(OutDir, 'DNAme.csv'), 
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = FALSE)

#######
### export metadata
######

expdat_list <- list("mRNA" = expdat_mRNA, "DNAme" = expdat_DNAme)

saveMetadata(OutDir = OutDir, 
             Seed = Seed, 
             smpls = smpls, 
             expdat_list = expdat_list, 
             brcds_list = NULL, 
             Prjcts = NULL, 
             topD = TopD,
             GeoMeans_list = NULL)

