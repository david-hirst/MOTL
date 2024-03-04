## FUNCTIONS

## --------------------------------------------------- OTHERS ----

selectProjects <- function(InputDir, PrjctExcl){
  #' Select project names 
  #'
  #' @param InputDir Directory name where unfiltered data are saved
  #' @param PrjctExcl Project names to exclude
  #' @returns Project names kept for the analysis as a list
  
  ## Remove excluded projects
  Prjcts = list.files(InputDir)
  Prjcts = strsplit(Prjcts,"_")
  Prjcts = as.data.frame(do.call("rbind",Prjcts))
  Prjcts = unique(Prjcts$V1)
  Prjcts = Prjcts[!is.element(Prjcts, PrjctExcl)]
  return(Prjcts)
}

## --------------------------------------------------------------

## ----------------------------------------- SELECT BARCODES ----

extractSampleNamesShared <- function(omics_list){
  #' Extract sample names shared between omics
  #' 
  #' @param omics_list list of omics data
  #' @returns list of samples that are in all omics
  
  smpls_list <- lapply(omics_list, function(omics_el){
    smpls = substr(colnames(omics_el),1,16)
    return(smpls)
  })
  smpls_shared = Reduce(base::intersect, smpls_list)
  
  return(smpls_shared)
}

selectMinBarcodeSampleNames <- function(expdat, smpls_shared, Prjct){
  #' Select minimum barcode for each sample name
  #' 
  #' @param expdat omic data
  #' @param smpls_shared list of samples that are in all omics
  #' @param Prjct TCGA project name
  #' @return list of barcodes selected, shared between all omics
  
  ## Extract and parse column names (sample names)
  smpls_names <- substr(colnames(expdat),1,16)
  ## Select the lower barcode for each sample (one replicate)
  brcds_expdat = data.frame(brcds = colnames(expdat)[is.element(smpls_names, smpls_shared)]) %>%
    dplyr::mutate(submittor = substr(brcds,1,12)) %>%
    dplyr::group_by(submittor) %>%
    dplyr::slice_min(brcds, n=1, with_ties = FALSE) %>%
    as.data.frame()
  ## Add project name to the table
  brcds_expdat$prjct = rep(Prjct,dim(brcds_expdat)[1])
  return(brcds_expdat)
}

## ----- READ FILES

mRNAFileParsing <- function(expdat_mRNA){
  #' Parse and filter out raw mRNA file for one project
  #' 
  #' Keep samples if na mean is < 0.2
  #' Keep samples if variance is > 0 
  #' 
  #' @param expdat_mRNA raw mRNA data 
  #' @return parsed mRNA data
  
  expdat_mRNA_tmp = assay(expdat_mRNA,1)
  SmplsKeep_mRNA = (colMeans(is.na(expdat_mRNA_tmp))<0.2) & (colVars(expdat_mRNA_tmp, na.rm = TRUE)>0)
  SmplsKeep_mRNA[is.na(SmplsKeep_mRNA)]=FALSE # I don't understand this line
  expdat_mRNA_tmp = expdat_mRNA_tmp[,SmplsKeep_mRNA]
  return(expdat_mRNA_tmp)
}

miRNAFileParsing <- function(expdat_miRNA){
  #' Parse and filter out raw miRNA file for one project
  #' 
  #' Keep samples if na meam is < 0.2
  #' Keep samples if varaiance is > 0
  #' 
  #' @param expdat_miRNA raw miRNA data 
  #' @return parsed miRNA data
  
  rownames(expdat_miRNA) = expdat_miRNA$miRNA_ID
  expdat_miRNA_tmp = expdat_miRNA[,substr(colnames(expdat_miRNA),1,10)=="read_count"]
  colnames(expdat_miRNA_tmp) = substr(colnames(expdat_miRNA_tmp),12,nchar(colnames(expdat_miRNA_tmp)))
  expdat_miRNA_tmp = as.matrix(expdat_miRNA_tmp)
  mode(expdat_miRNA_tmp) = "integer"
  SmplsKeep_miRNA = (colMeans(is.na(expdat_miRNA_tmp))<0.2) & (colVars(expdat_miRNA_tmp, na.rm = TRUE)>0)
  SmplsKeep_miRNA[is.na(SmplsKeep_miRNA)]=FALSE
  expdat_miRNA_tmp = expdat_miRNA_tmp[,SmplsKeep_miRNA]
  return(expdat_miRNA_tmp)
}

DNAmeFileParsing <- function(expdat_DNAme){
  #' Parse and transformed raw DNAme file for one project 
  #' 
  #' Transform Beta values into M values
  #' Keep samples if na meam is < 0.2
  #' Keep samples if varaiance is > 0 
  #' 
  #' @param expdat_DNAme raw DNAme data 
  #' @return parsed DNAme data
  
  expdat_DNAme_tmp = BetaValueToMValue(assay(expdat_DNAme,1))
  SmplsKeep_DNAme = (colMeans(is.na(expdat_DNAme_tmp))<0.2) & (colVars(expdat_DNAme_tmp, na.rm = TRUE)>0)
  SmplsKeep_DNAme[is.na(SmplsKeep_DNAme)]=FALSE
  expdat_DNAme_tmp = expdat_DNAme_tmp[,SmplsKeep_DNAme]
  return(expdat_DNAme_tmp)
}

SNVFileParsing <- function(expdat_SNV){
  #' Parse raw SNV file for one project
  #' 
  #' @param expdat_SNV raw SNV data 
  #' @return parsed SNV data
  
  ## Extract some columns from raw data
  expdat_SNV_parsed = expdat_SNV[
    is.element(expdat_SNV$Variant_Classification,
               c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                 "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                 "Splice_Site", "Translation_Start_Site")
    ) & !is.na(expdat_SNV$Hugo_Symbol),]
  
  # PrjctSffx = substr(Prjct,6,nchar(Prjct))
  
  # ## Create a data frame with the previous table
  # expdat_SNV_df = data.frame(
  #   Tumor_Sample_Barcode = expdat_SNV_tmp$Tumor_Sample_Barcode, 
  #   Hugo_Symbol = expdat_SNV_tmp$Hugo_Symbol
  #   #prjct = rep(PrjctSffx,dim(expdat_SNV_tmp)[1])
  # )
  
  ## convert to wide format matrix
  expdat_SNV_matrix = expdat_SNV_parsed %>%
    dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    dplyr::summarize(n = n()) %>%
    as.data.frame()
  
  ## Reshape SNV data to have samples in columns and features in rows
  expdat_SNV_reshape = reshape2::acast(expdat_SNV_matrix, Hugo_Symbol ~ Tumor_Sample_Barcode)
  
  # expdat_SNV_reshape = reshape(expdat_SNV_matrix, idvar = "Hugo_Symbol", timevar = "Tumor_Sample_Barcode", direction = "wide")
  # head(expdat_SNV_reshape[c(1:5)])
  
  return(expdat_SNV_reshape)
}

## --------------------------------------------------------------

## --------------------------------- OMIC DATA PREPROCESSING ----

## ----- READ AND PREPROCESS FILE

mRNAPreprocessing <- function(Prjct, brcds_mRNA, InputDir){
  #' Read and preprocess mRNA data from one TCGA project
  #'
  #' For one project in TCGA, read the mRNA data file.
  #' Then, create a SummarizedExperiment object with raw data.
  #' Select the samples that are in the brcds_mRNA object.
  #'
  #' @param Prjct project name
  #' @param brcds_mRNA data.frame that contains barcodes of mRNA data
  #' @param InputDir Directory input name where unfiltered data are saved
  #' @@returns tidy mRNA data with the selected barcode sample names (given in 
  #' the brcds_mRNA object). This return is a SummarizedExperiment object.
  
  PrjctSffx = substr(Prjct,6,nchar(Prjct))
  
  print(PrjctSffx)
  
  # import rds file
  expdat_mRNA_tmp = readRDS(file.path(InputDir,paste0(Prjct,"_mRNA.rds")))
  
  # tidy up mRNA
  expdat_mRNA_tmp = SummarizedExperiment(
    assays = list(counts = assay(expdat_mRNA_tmp,1)),
    rowRanges = expdat_mRNA_tmp@rowRanges,
    colData = data.frame(
      barcode = colnames(expdat_mRNA_tmp),
      batch_brcd_drvd = substr(colnames(expdat_mRNA_tmp),22,28),
      prjct = rep(PrjctSffx,dim(expdat_mRNA_tmp)[2])
    ))
  
  # samplewise filtering
  SmplsKeep_mRNA = is.element(colnames(assay(expdat_mRNA_tmp)), brcds_mRNA$brcds)
  expdat_mRNA_tmp = expdat_mRNA_tmp[,SmplsKeep_mRNA]
  
  gc()
  
  return(expdat_mRNA_tmp)
}

miRNAPreprocessing <- function(Prjct, brcds_miRNA, InputDir){
  #' Read and preprocess miRNA data from one TCGA project
  #'
  #' For one project in TCGA, read the miRNA data file.
  #' Then, create a SummarizedExperiment object with raw data.
  #' Select the samples that are in the brcds_miRNA object.
  #'
  #' @param Prjct project name
  #' @param brcds_miRNA data.frame that contains barcodes of miRNA data
  #' @param InputDir Directory input name where unfiltered data are saved
  #' @@returns tidy miRNA data with the selected barcode sample names (given in 
  #' the brcds_miRNA object). This return is a SummarizedExperiment object.
  
  PrjctSffx = substr(Prjct,6,nchar(Prjct))
  
  print(PrjctSffx)
  
  # import rds file
  expdat_miRNA_tmp = readRDS(file.path(InputDir,paste0(Prjct,"_miRNA.rds")))
  
  # tidy up miRNA
  rownames(expdat_miRNA_tmp) = expdat_miRNA_tmp$miRNA_ID
  expdat_miRNA_tmp = expdat_miRNA_tmp[,substr(colnames(expdat_miRNA_tmp),1,10)=="read_count"]
  colnames(expdat_miRNA_tmp) = substr(colnames(expdat_miRNA_tmp),12,nchar(colnames(expdat_miRNA_tmp)))
  expdat_miRNA_tmp = as.matrix(expdat_miRNA_tmp)
  mode(expdat_miRNA_tmp) = "integer"
  
  expdat_miRNA_tmp = SummarizedExperiment(
    assays = list(counts = expdat_miRNA_tmp),
    colData = data.frame(
      barcode = colnames(expdat_miRNA_tmp),
      batch_brcd_drvd = substr(colnames(expdat_miRNA_tmp),22,28),
      prjct = rep(PrjctSffx,dim(expdat_miRNA_tmp)[2])
    )
  )
  
  # samplewise filtering
  SmplsKeep_miRNA = is.element(colnames(assay(expdat_miRNA_tmp)), brcds_miRNA$brcds)
  expdat_miRNA_tmp = expdat_miRNA_tmp[,SmplsKeep_miRNA]
  
  gc()
  
  return(expdat_miRNA_tmp)
}

DNAmePreprocessing <- function(Prjct, brcds_DNAme, InputDir){
  #' Read and preprocess DNAme data from one TCGA project
  #'
  #' For one project in TCGA, read the DNAme data file.
  #' Then, beta value are transformed to M values. 
  #' Create a SummarizedExperiment object with the M transformed values.
  #' Select the samples that are in the brcds_DNAme object.
  #' 
  #' @param Prjct project name
  #' @param brcds_DNAme data.frame that contains barcodes of DNAme data
  #' @param InputDir Directory input name where unfiltered data are saved
  #' @@returns tidy DNAme data with the selected barcode sample names (given in 
  #' the brcds_DNAme object). This return is a SummarizedExperiment object.
  
  PrjctSffx = substr(Prjct,6,nchar(Prjct))
  
  print(PrjctSffx)
  
  gc()
  
  # import rds file
  expdat_DNAme_tmp = readRDS(file.path(InputDir,paste0(Prjct,"_DNAme.rds")))
  
  # tidy up DNAme
  expdat_DNAme_tmp = SummarizedExperiment(
    assays = list(MValues = sesame::BetaValueToMValue(assay(expdat_DNAme_tmp,1))),
    rowRanges = expdat_DNAme_tmp@rowRanges,
    colData = data.frame(
      barcode = colnames(expdat_DNAme_tmp),
      batch_brcd_drvd = substr(colnames(expdat_DNAme_tmp),22,28),
      prjct = rep(PrjctSffx,dim(expdat_DNAme_tmp)[2])
    )
  )
  
  # samplewise filtering
  SmplsKeep_DNAme = is.element(colnames(assay(expdat_DNAme_tmp)), brcds_DNAme$brcds)
  expdat_DNAme_tmp = expdat_DNAme_tmp[,SmplsKeep_DNAme]
  
  gc()
  
  return(expdat_DNAme_tmp)
}

SNVPreprocessing <- function(Prjct, brcds_SNV, InputDir){
  #' Read and preprocess SNV data from one TCGA project
  #'
  #' For one project in TCGA, read the SNV data file.
  #' Then, create a data.frame with raw data.
  #' Select the samples that are in the brcds_SNV object.
  #'
  #' @param Prjct project name
  #' @param brcds_SNV data.frame that contains barcodes of SNV data
  #' @param InputDir Directory input name where unfiltered data are saved
  #' @@returns tidy SNV data with the selected barcode sample names (given in 
  #' the brcds_SNV object). This return is data.frame object.
  
  PrjctSffx = substr(Prjct,6,nchar(Prjct))
  
  print(PrjctSffx)
  
  # import rds file
  expdat_SNV_tmp = readRDS(file.path(InputDir,paste0(Prjct,"_SNV.rds")))
  
  # tidy up SNV
  expdat_SNV_tmp = expdat_SNV_tmp[
    is.element(expdat_SNV_tmp$Variant_Classification,
               c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", 
                 "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
                 "Splice_Site", "Translation_Start_Site")
    ) & !is.na(expdat_SNV_tmp$Hugo_Symbol),]
  
  expdat_SNV_tmp = data.frame(
    Tumor_Sample_Barcode = expdat_SNV_tmp$Tumor_Sample_Barcode, 
    Hugo_Symbol = expdat_SNV_tmp$Hugo_Symbol,
    prjct = rep(PrjctSffx,dim(expdat_SNV_tmp)[1])
  )
  
  ### subset to only include data from these barcodes
  SmplsKeep_SNV = is.element(expdat_SNV_tmp$Tumor_Sample_Barcode, brcds_SNV$brcds)
  expdat_SNV_tmp = expdat_SNV_tmp[SmplsKeep_SNV,]
  
  gc()
  
  return(expdat_SNV_tmp)
}

## ----- READ (CALL PREVIOUS FUNCTIONS) AND MERGE FILES INTO ONE OBJECT

mRNAdataPreparation <- function(Prjcts, brcds_mRNA, InputDir){
  #'
  #' @param Prjct project name
  #'
  #'
  
  for(i in 1:length(Prjcts)){
    ## SELECT PROJECT
    Prjct <- Prjcts[i]
    ## READ DATA
    exdat_mRNA <- mRNAPreprocessing(Prjct, brcds_mRNA, InputDir)
    if(i == 1){
      expdat_mRNA_merged = exdat_mRNA
    }else{
      exdat_mRNA = exdat_mRNA[match(rownames(expdat_mRNA_merged),rownames(exdat_mRNA)),]
      expdat_mRNA_merged = SummarizedExperiment::cbind(expdat_mRNA_merged, exdat_mRNA)
    }
    rm(exdat_mRNA)
    gc()
  }
  return(expdat_mRNA_merged)
}

miRNAdataPreparation <- function(Prjcts, brcds_miRNA, InputDir){
  #'
  #' @param Prjct project name
  #'
  #'
  
  for(i in 1:length(Prjcts)){
    ## SELECT PROJECT
    Prjct <- Prjcts[i]
    ## READ DATA
    exdat_miRNA <- miRNAPreprocessing(Prjct, brcds_miRNA, InputDir)
    if(i == 1){
      expdat_miRNA_merged = exdat_miRNA
    }else{
      exdat_miRNA = exdat_miRNA[match(rownames(expdat_miRNA_merged),rownames(exdat_miRNA)),]
      expdat_miRNA_merged = SummarizedExperiment::cbind(expdat_miRNA_merged, exdat_miRNA)
    }
    rm(exdat_miRNA)
    gc()
  }
  return(expdat_miRNA_merged)
}

DNAmedataPreparation <- function(Prjcts, brcds_DNAme, InputDir){
  #'
  #' @param Prjct project name
  #'
  #'
  
  for(i in 1:length(Prjcts)){
    ## SELECT PROJECT
    Prjct <- Prjcts[i]
    ## READ DATA
    exdat_DNAme <- DNAmePreprocessing(Prjct, brcds_DNAme, InputDir)
    if(i == 1){
      exdat_DNAme_merged = exdat_DNAme
    }else{
      exdat_DNAme = exdat_DNAme[match(rownames(exdat_DNAme_merged),rownames(exdat_DNAme)),]
      exdat_DNAme_merged = SummarizedExperiment::cbind(exdat_DNAme_merged, exdat_DNAme)
    }
    rm(exdat_DNAme)
    gc()
  }
  return(exdat_DNAme_merged)
}

SNVdataPreparation <- function(Prjcts, brcds_SNV, InputDir){
  #'
  #' @param Prjct project name
  #'
  #'
  
  for(i in 1:length(Prjcts)){
    ## SELECT PROJECT
    Prjct <- Prjcts[i]
    ## READ DATA
    exdat_SNV <- SNVPreprocessing(Prjct, brcds_SNV, InputDir)
    if(i == 1){
      exdat_SNV_merged = exdat_SNV
    }else{
      exdat_SNV_merged = rbind(exdat_SNV_merged, exdat_SNV)
    }
    rm(exdat_SNV)
    gc()
  }
  
  ## convert to wide format matrix
  exdat_SNV_merged = exdat_SNV_merged %>%
    dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
    dplyr::summarize(n = n()) %>%
    as.data.frame()
  
  ## Reshape
  exdat_SNV_merged = reshape(exdat_SNV_merged,
                             idvar = "Hugo_Symbol",
                             timevar = "Tumor_Sample_Barcode",
                             direction = "wide")
  
  rownames(exdat_SNV_merged) = exdat_SNV_merged$Hugo_Symbol
  exdat_SNV_merged = exdat_SNV_merged[,substr(colnames(exdat_SNV_merged),1,2)=="n."]
  colnames(exdat_SNV_merged) = substr(colnames(exdat_SNV_merged),3,nchar(colnames(exdat_SNV_merged)))
  
  exdat_SNV_merged = as.matrix(exdat_SNV_merged)
  exdat_SNV_merged[is.na(exdat_SNV_merged)]=0
  exdat_SNV_merged = (exdat_SNV_merged > 0)+0
  mode(exdat_SNV_merged) = "integer"
  
  return(exdat_SNV_merged)
}

## MT - TO REMOVE
mergeExperimentalData <- function(expdat_list){
  #' Merge SE object of several TCGA projects into one SE object.
  #' 
  #' SE: SummarizedExperiment
  #' Create a row_range list for each data.
  #' Merge then, and remove duplicates.
  #' Create a single SE with all data and the merged row_ranges. 
  #' 
  #' @param expdat_list list of experimental data of different projects (SE object)
  #' @returns SummarizedExperiment object with the merged experimental data
  
  ## MT - Check if select shared features
  expdat_fltr <- mergeSEs(expdat_list, do.scale = FALSE, commonOnly = TRUE, addDatasetPrefix = FALSE)
  
  ## Create a rowrange list
  rowRanges_list <- lapply(expdat_list, function(expdat){
    rowRanges_df <- as.data.frame(expdat@rowRanges)
    rowRanges_df$rownames <- row.names(rowRanges_df)
    return(rowRanges_df)
  })
  rowRanges_merged_tmp <- do.call(rbind, c(rowRanges_list, make.row.names = FALSE))
  rowRanges_merged_tmp <- unique(rowRanges_merged_tmp)
  rownames(rowRanges_merged_tmp) <- rowRanges_merged_tmp$rownames
  rowRanges_merged <- GRanges(rowRanges_merged_tmp)
  
  ##
  expdat_merged <- SummarizedExperiment(assays = SimpleList(counts = assay(expdat_fltr)),
                                        rowRanges = rowRanges_merged[row.names(expdat_fltr)],
                                        colData = expdat_fltr@colData,
                                        metadata = expdat_fltr@metadata)
  
  return(expdat_merged)
}

## MT - TO REMOVE
mergeSEs_block <- function(expdat_list, SE){
  if(length(expdat_list) == 0){
    return(SE)
  }
  expdat = expdat_list[[1]]
  expdat_list = expdat_list[-1]
  gc()
  if(length(SE) < 1){
    SE = expdat
  }else{
    expdat = expdat[match(rownames(SE),rownames(expdat)),]
    SE = SummarizedExperiment::cbind(SE, expdat)
  }
  rm(expdat)
  gc()
  return(mergeSEs_block(expdat_list, SE, nb+1))
}

## MT - TO REMOVE
createRowRanges <- function(expdat_list){
  rowRanges_list <- lapply(expdat_list, function(expdat){
    rowRanges_df <- as.data.frame(expdat@rowRanges)
    rowRanges_df$rownames <- row.names(rowRanges_df)
    return(rowRanges_df)
  })
  
  ## Create a rowrange list
  rowRanges_merged_tmp <- do.call(rbind, c(rowRanges_list, make.row.names = FALSE))
  rowRanges_merged_tmp <- unique(rowRanges_merged_tmp)
  rownames(rowRanges_merged_tmp) <- rowRanges_merged_tmp$rownames
  rowRanges_merged <- GRanges(rowRanges_merged_tmp)
  
  rm(rowRanges_merged_tmp)
  gc()
  return(rowRanges_merged)
}

## MT - TO REMOVE
mergeExperimentalData_1 <- function(expdat_list){
  #' Merge SE object of several TCGA projects into one SE object.
  #'
  #' SE: SummarizedExperiment
  #' Create a row_range list for each data.
  #' Merge then, and remove duplicates.
  #' Create a single SE with all data and the merged row_ranges.
  #'
  #' @param expdat_list list of experimental data of different projects (SE object)
  #' @returns SummarizedExperiment object with the merged experimental data
  
  
  
  ## MT - Check if select shared features
  ## expdat_fltr <- mergeSEs(expdat_list, do.scale = FALSE, commonOnly = TRUE, addDatasetPrefix = FALSE)
  SE = list()
  block_size = 3
  SE = mergeSEs_block(expdat_list, SE, block_size)
  
  gc()
  
  return(SE)
}

## ----- FILTER DATA

mRNAFiltering <- function(expdat_mRNA){
  #' Filter mRNA experimental data
  #'
  #' Keep features (rows) with:
  #'  - missing value (NA) rate < 0.2
  #'  - non chromosome Y localisation
  #'  - rate of zero value < 0.9
  #'  - variability rate > 0
  #'  
  #' @param expdat_mRNA SE object of mRNA data
  #' @returns SE mRNA data object filtered
  
  FtrsKeep_mRNA = (rowMeans(is.na(assay(expdat_mRNA)))<0.2) &
    (!is.element(seqnames(expdat_mRNA),c("chrY"))) &
    (rowMeans(assay(expdat_mRNA)==0, na.rm = TRUE)<0.9) &
    (rowVars(assay(expdat_mRNA), na.rm = TRUE, useNames = FALSE)>0) &
    (!regexpr('[_]',rownames(assay(expdat_mRNA)))>0)
  
  FtrsKeep_mRNA[is.na(FtrsKeep_mRNA)]=FALSE
  
  expdat_mRNA = expdat_mRNA[FtrsKeep_mRNA,]
  
  return(expdat_mRNA)
}

miRNAFiltering <- function(expdat_miRNA){
  #' Filter miRNA experimental data
  #'
  #' Keep features (rows) with:
  #'  - missing value (NA) rate < 0.2
  #'  - rate of zero value < 0.9
  #'  - variability rate > 0
  #'  
  #' @param expdat_mRNA SE object of miRNA data
  #' @returns SE miRNA data object filtered
  
  FtrsKeep_miRNA = (rowMeans(is.na(assay(expdat_miRNA)))<0.2) &
    (rowMeans(assay(expdat_miRNA)==0, na.rm = TRUE)<0.9) &
    (rowVars(assay(expdat_miRNA), na.rm = TRUE)>0)
  
  FtrsKeep_miRNA[is.na(FtrsKeep_miRNA)]=FALSE
  
  expdat_miRNA = expdat_miRNA[FtrsKeep_miRNA,]
  
  return(expdat_miRNA)
}

DNAmeFiltering <- function(expdat_DNAme, TopD){
  #' Filter DNAme experimental data
  #'
  #' Keep features (rows) with:
  #'  - non chromosome Y and X localisation
  #'  - missing value (NA) rate < 0.2
  #'  - variability rate > 0
  #' Then, select the top of features based on the variability
  #'  
  #' @param expdat_DNAme SE object of DNAme data
  #' @param TopD number of features to keep
  #' @returns data.frame of the DNAme filtered data
  
  ## prefiltering of merged dataset
  FtrsKeep_DNAme = (!as.data.frame(rowData(expdat_DNAme))$MASK_general) &
    (!is.element(seqnames(expdat_DNAme),c("chrY","chrX"))) &
    (rowMeans(is.na(assay(expdat_DNAme)))<0.2) &
    (rowVars(assay(expdat_DNAme), na.rm = TRUE)>0)
  
  FtrsKeep_DNAme[is.na(FtrsKeep_DNAme)]=FALSE
  
  expdat_DNAme_tmp = expdat_DNAme[FtrsKeep_DNAme,]
  
  ## keep only most variable features
  FtrsKeep_DNAme_fltr = base::rank(-rowVars(assay(expdat_DNAme_tmp), na.rm = TRUE), ties.method = "first") <= TopD
  expdat_DNAme_fltr = expdat_DNAme_tmp[FtrsKeep_DNAme_fltr,] 
  
  expdat_DNAme_assay = assay(expdat_DNAme_fltr)
  
  return(expdat_DNAme_assay)
}

SNVFiltering <- function(expdat_SNV, TopD, SNVthreshold=0.01){
  #' Filter SNV experimental data
  #'
  #' Data are reshape to have same shape of other data (features in rows and
  #' samples in columns). Only necessary information are kept. 
  #' Keep features with variability rate > 0 and above the threshold. 
  #' Then, select the top of features based on the variability
  #'  
  #' @param expdat_SNV SE object of SNV data
  #' @param SNVthreshold threshold to keep features (default 0.01)
  #' @param TopD number of features to keep
  #' @returns data.frame of the SNV filtered data
  
  # ## convert to wide format matrix
  # expdat_SNV = expdat_SNV %>%
  #   dplyr::group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  #   dplyr::summarize(n = n()) %>%
  #   as.data.frame()
  # 
  # ## Reshape
  # expdat_SNV = reshape(expdat_SNV,
  #                      idvar = "Hugo_Symbol",
  #                      timevar = "Tumor_Sample_Barcode",
  #                      direction = "wide")
  # 
  # rownames(expdat_SNV) = expdat_SNV$Hugo_Symbol
  # expdat_SNV = expdat_SNV[,substr(colnames(expdat_SNV),1,2)=="n."]
  # colnames(expdat_SNV) = substr(colnames(expdat_SNV),3,18)
  # 
  # expdat_SNV_mtrx = as.matrix(expdat_SNV)
  # expdat_SNV_mtrx[is.na(expdat_SNV_mtrx)]=0
  # expdat_SNV_mtrx = (expdat_SNV_mtrx > 0)+0
  # mode(expdat_SNV_mtrx) = "integer"
  
  ## prefiltering of merged dataset
  FtrsKeep_SNV = (rowVars(expdat_SNV, na.rm = TRUE)>0) & 
    (rowMeans(expdat_SNV)>SNVthreshold)
  
  FtrsKeep_SNV[is.na(FtrsKeep_SNV)]=FALSE
  expdat_SNV_tmp = expdat_SNV[FtrsKeep_SNV,]
  
  ## keep X most variable and tidy up
  FtrsKeep_SNV_fltr = base::rank(-rowVars(expdat_SNV_tmp, na.rm = TRUE), ties.method = "first") <= TopD
  expdat_SNV_fltr = expdat_SNV_tmp[FtrsKeep_SNV_fltr,]
  
  return(expdat_SNV_fltr)
}

## ----- NORMALIZE DATA

GeoMeanFun = function(x){
  ## MT ADD DESCRIPTION
  #' Short description
  #' 
  #' Detailed description
  #' 
  #' @param x vector of numeric values
  #' @return mean of non zero values from x vector
  
  GeoMeans <- exp(sum(log(x[x > 0]))/length(x))
  return(GeoMeans)
}

countsNormalization <- function(expdat, GeoMeans){
  #' Normalize counts data
  #' 
  #' Normalize counts data using DESeq2 normalization. 
  #' Calculate geometric means for normalization if data = learning set
  #' Doesn't use geometric means for normalization if data = target set
  #' Use provided geometric means for normalization if transfer learning
  #' 
  #' @param expdat SE object of experimental data (could be miRNA or mRNA)
  #' @param GeoMeans "Trg", "Lrn" or vector of numerics
  #' @returns list of data.frame of the counts normalized and GeoMeans calculated
  
  ## create deseq object
  expdat_dds = DESeq2::DESeqDataSet(expdat, design = ~ 1)
  
  if(is.numeric(GeoMeans)){
    ## normalization
    expdat_dds_norm = DESeq2::estimateSizeFactors(expdat_dds, geoMeans = GeoMeans)
    GeoMeans = NULL
  }else if(GeoMeans == "Lrn"){
    ## calculate geometric means to use for normalization of both learning and target sets
    GeoMeans = apply(counts(expdat_dds),1,GeoMeanFun)
    ## normalization
    expdat_dds_norm = DESeq2::estimateSizeFactors(expdat_dds, geoMeans = as.vector(GeoMeans))
  } else if(GeoMeans == "Trg"){
    ## estimate size factors
    expdat_dds_norm = DESeq2::estimateSizeFactors(expdat_dds)
    GeoMeans = NULL
  }else{
    print("GeoMeans parameter should be 'Trg' or 'Lrn' or a numeric value")
    stop()
  }
  
  ## Extract normalized counts
  expdat_counts_norm <- list("counts" = counts(expdat_dds_norm, normalized = TRUE))
  
  ## Save GeoMeans
  expdat_counts_norm$GeoMeans <- GeoMeans
  
  return(expdat_counts_norm)
}

countsTransformation <- function(expdat_count, TopD){
  #' Log2 Transform and select top data based on variance
  #' 
  #' @param expdat_count data.frame of the counts
  #' @param TopD number of features to keep
  #' @returns data.frame of the log2 transformed and filtered data
  
  ## log transform and filter to keep only most variable
  expdat_counts_log = log2(expdat_count+1)
  FtrsKeep = base::rank(-rowVars(expdat_counts_log, na.rm = TRUE, useNames = FALSE), ties.method = "first") <= TopD
  expdat_counts_fltr = expdat_counts_log[FtrsKeep,]
  return(expdat_counts_fltr)
}

## --------------------------------------------------------------

## ---------------------------------- SAVE DATA AND METADATA ----

orderAndSaveData <- function(expdat, smpls, OutDir, fileName){
  #' Save data into a csv file
  #'
  #' @param expdat data.frame of experimental data
  #' @param smpls sample names vector with specific order
  #' @param OutDir output directory name
  #' @param fileName output file name (with .csv)
  #' @returns data.frame of the sorted data
  
  ## tidy up names 
  colnames(expdat) = substr(colnames(expdat),1,16) 
  
  ## ORDER COLUMNS
  expdat_sorted <- expdat[, smpls]
  
  ## save object
  write.table(expdat_sorted, file.path(OutDir, fileName), 
              quote = FALSE, sep=',', na='', row.names = FALSE, col.names = FALSE)
  
  return(expdat_sorted)
}

saveMetadata <- function(OutDir, Seed, smpls, expdat_list, brcds_list, PCA_list=NULL, Elbow_list=NULL, Prjcts, topD=NULL, GeoMeans_list=NULL){
  #' Save the metadata into rds and json files
  #' 
  #' @param OutDir output directory name
  #' @param Seed seed 
  #' @param smpls sample names vector with specific order
  #' @param expdat_list list of normalized and transformed data
  #' @param brcds_list list of samples data.frame of each data
  #' @param PCA_list list of calculated PCA values for each data
  #' @param Elbow_list list of calculated elbow values for each data
  #' @param Prjcts list of projects 
  #' @param topD top of features selected
  #' @param GeoMeans_list list of calculated GeoMeans for counts data
  
  ftrs_list <-  lapply(expdat_list, row.names)
  names(ftrs_list) <- paste("ftrs", names(ftrs_list), sep = "_")
  
  expdat_meta <- c(ftrs_list, brcds_list, PCA_list, Elbow_list, GeoMeans_list)
  expdat_meta$ElbowK_Total = sum(unlist(Elbow_list))
  expdat_meta$Seed = Seed
  expdat_meta$smpls = smpls
  expdat_meta$Prjcts = Prjcts
  expdat_meta$TopD = TopD
  
  ## SAVE METADATA TO RDS FILE
  saveRDS(expdat_meta,file.path(OutDir,'expdat_meta.rds'))
  
  ## SAVE METADATA TO JSON FILE
  expdat_meta.json = rjson::toJSON(expdat_meta)
  write(expdat_meta.json,file.path(OutDir,'expdat_meta.json'))
  
}

## --------------------------------------------------------------

## ------------------------------------ SUBSET DATA CREATION ----

selectBarcodeSampleNames <- function(smpls, expdat){
  #' Select samples from experimental data
  #'
  #' @param smpls sample names to select
  #' @param expdat experimental data
  #' @returns experimental data with selected samples
  
  brcds_sbstr <- substr(expdat$brcds,1,16)
  brcds_fltr <- expdat[is.element(brcds_sbstr, smpls),]
  return(brcds_fltr)
}

createSubsetOneProject <- function(SSnb, brcds_SS, if_SNV = FALSE, OutDir, GeoMeans, InputDir){
  #' Prepare data for one subset with only one project
  #' 
  #' 1. Extract barcodes for the corresponding subset
  #' 2. Prepare mRNA data
  #' 3. Prepare miRNA data
  #' 4. Prepapre DNAme data
  #' 5. (optional) Prepare SNV data
  #' 6. Save results and metadata
  #'
  #' @param SSnb number of subset to create
  #' @param brcds_SS barcodes of select samples for each subset
  #' @param if_SNV boolean to process SNV pipeline 
  #' @param OutDir folder to save preprocessed data
  
  SS = paste0('SS_',SSnb)
  
  print("[*]-------------------------------------")
  print(SS)
  
  ## directory to export data to
  SSOutDir = file.path(OutDir,SS)
  if(!dir.exists(SSOutDir)){
    dir.create(SSOutDir, showWarnings = FALSE, recursive = TRUE)
  }
  
  smpls_SS = brcds_SS$smpls_SS[[SS]]
  brcds_mRNA_SS = brcds_SS$brcds_mRNA_SS[[SS]]
  brcds_miRNA_SS = brcds_SS$brcds_miRNA_SS[[SS]]
  brcds_DNAme_SS = brcds_SS$brcds_DNAme_SS[[SS]]
  
  ## ---------- mRNA data preprocessing
  
  print("mRNA data preprocessing")
  
  ## Read experimental data files
  # expdat_mRNA = mRNAPreprocessing(Prjct = Prjct, brcds_mRNA = brcds_mRNA_SS, InputDir = InputDir)
  ## Filtering
  # expdat_mRNA_fltr <- mRNAFiltering(expdat_mRNA = expdat_mRNA)
  
  ## Read and merge files
  expdat_mRNA_merged <- mRNAdataPreparation(Prjcts = Prjct, brcds_mRNA = brcds_mRNA_SS, InputDir = InputDir)
  
  ## Normalization
  expdat_mRNA_norm <- countsNormalization(expdat = expdat_mRNA_merged, GeoMeans = GeoMeans)
  
  ## Transformation
  expdat_mRNA_log <- countsTransformation(expdat_count = expdat_mRNA_norm$counts, TopD = TopD)
  
  ## Save data
  expdat_mRNA_sorted <- orderAndSaveData(expdat = expdat_mRNA_log, smpls = smpls_SS, OutDir = SSOutDir, fileName = "mRNA.csv")
  
  ## CLEAN SPACE
  rm(expdat_mRNA_merged, expdat_mRNA_norm, expdat_mRNA_log)
  gc()
  
  ## ---------- miRNA data preprocessing
  
  print("miRNA data preprocessing")
  
  ## Read experimental data files
  # expdat_miRNA = miRNAPreprocessing(Prjct = Prjct, brcds_miRNA = brcds_miRNA_SS, InputDir = InputDir)
  ## Filtering
  # expdat_miRNA_fltr <- miRNAFiltering(expdat_miRNA = expdat_miRNA)
  
  ## Read and merge files
  expdat_miRNA_merged <- miRNAdataPreparation(Prjcts = Prjct, brcds_miRNA = brcds_miRNA_SS, InputDir = InputDir)
  
  ## Normalization
  expdat_miRNA_norm <- countsNormalization(expdat = expdat_miRNA_merged, GeoMeans = GeoMeans)
  
  ## Transformation
  expdat_miRNA_log <- countsTransformation(expdat_count = expdat_miRNA_norm$counts, TopD = TopD)
  
  ## Save data
  expdat_miRNA_sorted <- orderAndSaveData(expdat = expdat_miRNA_log, smpls = smpls_SS, OutDir = SSOutDir, fileName = "miRNA.csv")
  
  ## CLEAN SPACE
  rm(expdat_miRNA_merged, expdat_miRNA_norm, expdat_miRNA_log)
  gc()
  
  ## ---------- DNAme data preprocessing
  
  print("DNAme data preprocessing")
  
  ## Read experimental data files
  # expdat_DNAme <- DNAmePreprocessing(Prjct = Prjct, brcds_DNAme = brcds_DNAme_SS, InputDir = InputDir)
  
  ## Read and merge files 
  expdat_DNAme_merged <- DNAmedataPreparation(Prjcts = Prjct, brcds_DNAme = brcds_DNAme_SS, InputDir = InputDir)
  
  ## Filtering
  expdat_DNAme_fltr <- DNAmeFiltering(expdat_DNAme = expdat_DNAme_merged, TopD = TopD)
  
  ## Save data
  expdat_DNAme_sorted <- orderAndSaveData(expdat = expdat_DNAme_fltr, smpls = smpls_SS, OutDir = SSOutDir, fileName = "DNAme.csv")
  
  # CLEAN SPACE
  rm(expdat_DNAme_merged, expdat_DNAme_fltr)
  gc()
  
  ## ---------- SNV data preprocessing
  
  if(if_SNV){
    
    brcds_SNV_SS = brcds_SS$brcds_SNV_SS[[SS]]
    
    ## Read experimental data files
    # expdat_SNV <- SNVPreprocessing(Prjct = Prjct, brcds_SNV = brcds_SNV_SS, InputDir = InputDir)
    
    ## Reand and merge files
    expdat_SNV_merged <- SNVdataPreparation(Prjcts = Prjct, brcds_SNV = brcds_SNV_SS, InputDir = InputDir)
    
    ## Filtering
    expdat_SNV_fltr <- SNVFiltering(expdat_SNV = expdat_SNV_merged, TopD = TopD)
    
    ## Save data
    expdat_SNV_sorted <- orderAndSaveData(expdat = expdat_SNV_fltr, smpls = smpls_SS, OutDir = SSOutDir, fileName = "SNV.csv")
    
    ## CLEAN SPACE
    rm(expdat_SNV_merged, expdat_SNV_fltr)
    gc()
  }
  
  ## ---------- SAVE METADATA
  print("Save metadata")
  
  if(if_SNV){
    
    expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted, "SNV" = expdat_SNV_sorted)
    brcds_list <- list("brcds_mRNA" = brcds_mRNA_SS, "brcds_miRNA" = brcds_miRNA_SS, "brcds_DNAme" = brcds_DNAme_SS, "brcds_SNV" = brcds_SNV_SS)
  }
  else{
    expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted)
    brcds_list <- list("brcds_mRNA" = brcds_mRNA_SS, "brcds_miRNA" = brcds_miRNA_SS, "brcds_DNAme" = brcds_DNAme_SS)
  }
  
  saveMetadata(OutDir = SSOutDir, Seed = Seed, smpls = smpls_SS, expdat_list = expdat_list, brcds_list = brcds_list, Prjcts = Prjcts, GeoMeans_list = NULL)
}

createSubsetMultiProjects <- function(SSnb, brcds_SS, if_SNV = FALSE, OutDir, GeoMeans, InputDir){
  #' Prepare data for one subset for several projects
  #' 
  #' 1. Extract barcodes for the corresponding subset
  #' 2. Prepare mRNA data
  #' 3. Prepare miRNA data
  #' 4. Prepapre DNAme data
  #' 5. (optional) Prepare SNV data
  #' 6. Save results and metadata
  #'
  #' @param SSnb number of subset to create
  #' @param brcds_SS barcodes of select samples for each subset
  #' @param if_SNV boolean to process SNV pipeline 
  #' @param OutDir folder to save preprocessed data
  
  SS = paste0('SS_',SSnb)
  
  print("[*]-------------------------------------")
  print(SS)
  
  ## directory to export data to
  SSOutDir = file.path(OutDir,SS)
  if(!dir.exists(SSOutDir)){
    dir.create(SSOutDir, showWarnings = FALSE, recursive = TRUE)
  }
  
  smpls_SS = brcds_SS$smpls_SS[[SS]]
  brcds_mRNA_SS = brcds_SS$brcds_mRNA_SS[[SS]]
  brcds_miRNA_SS = brcds_SS$brcds_miRNA_SS[[SS]]
  brcds_DNAme_SS = brcds_SS$brcds_DNAme_SS[[SS]]
  
  ## ---------- mRNA data preprocessing
  
  print("mRNA data preprocessing")
  
  ## Read experimental data files
  # expdat_mRNA_list = sapply(Prjcts, mRNAPreprocessing, brcds_mRNA = brcds_mRNA_SS, InputDir = InputDir)
  ## Merge experimental data into only one SE object
  # expdat_mRNA_merged <- mergeExperimentalData(expdat_list = expdat_mRNA_list)
  
  ## Read and merge files
  expdat_mRNA_merged <- mRNAdataPreparation(Prjcts = Prjcts, brcds_mRNA = brcds_mRNA_SS, InputDir = InputDir)
  
  ## Filtering
  expdat_mRNA_fltr <- mRNAFiltering(expdat_mRNA = expdat_mRNA_merged)
  
  ## Normalization
  expdat_mRNA_norm <- countsNormalization(expdat = expdat_mRNA_fltr, GeoMeans = GeoMeans)
  
  ## Transformation
  expdat_mRNA_log <- countsTransformation(expdat_count = expdat_mRNA_norm$counts, TopD = TopD)
  
  ## Save data
  expdat_mRNA_sorted <- orderAndSaveData(expdat = expdat_mRNA_log, smpls = smpls_SS, OutDir = SSOutDir, fileName = "mRNA.csv")
  
  ## CLEAN SPACE
  rm(expdat_mRNA_merged, expdat_mRNA_fltr, expdat_mRNA_norm, expdat_mRNA_log)
  gc()
  
  ## ---------- miRNA data preprocessing
  
  print("miRNA data preprocessing")
  
  ## Read experimental data files
  # expdat_miRNA_list = sapply(Prjcts, miRNAPreprocessing, brcds_miRNA = brcds_miRNA_SS, InputDir = InputDir)
  ## Merge experimental data into only one SE object
  # expdat_miRNA_merged <- mergeSEs(expdat_miRNA_list, do.scale = FALSE, commonOnly = TRUE, addDatasetPrefix = FALSE)
  
  ## Read and merge files
  expdat_miRNA_merged <- miRNAdataPreparation(Prjcts = Prjcts, brcds_miRNA = brcds_miRNA_SS, InputDir = InputDir)
  
  ## Filtering
  expdat_miRNA_fltr <- miRNAFiltering(expdat_miRNA = expdat_miRNA_merged)
  
  ## Normalization
  expdat_miRNA_norm <- countsNormalization(expdat = expdat_miRNA_fltr, GeoMeans = GeoMeans)
  
  ## Transformation
  expdat_miRNA_log <- countsTransformation(expdat_count = expdat_miRNA_norm$counts, TopD = TopD)
  
  ## Save data
  expdat_miRNA_sorted <- orderAndSaveData(expdat = expdat_miRNA_log, smpls = smpls_SS, OutDir = SSOutDir, fileName = "miRNA.csv")
  
  ## CLEAN SPACE
  rm(expdat_miRNA_merged, expdat_miRNA_fltr, expdat_miRNA_norm, expdat_miRNA_log)
  gc()
  
  ## ---------- DNAme data preprocessing
  
  print("DNAme data preprocessing")
  
  ## Read experimental data files
  # expdat_DNAme_list <- sapply(Prjcts, DNAmePreprocessing, brcds_DNAme = brcds_DNAme_SS, InputDir = InputDir)
  ## Merge experimental data into only one SE object
  # expdat_DNAme_merged <- mergeExperimentalData(expdat_list = expdat_DNAme_list)
  
  ## Read and merge files
  expdat_DNAme_merged <- DNAmedataPreparation(Prjcts = Prjcts, brcds_DNAme = brcds_DNAme_SS, InputDir = InputDir)
  
  ## Filtering
  expdat_DNAme_fltr <- DNAmeFiltering(expdat_DNAme = expdat_DNAme_merged, TopD = TopD)
  
  ## Save data
  expdat_DNAme_sorted <- orderAndSaveData(expdat = expdat_DNAme_fltr, smpls = smpls_SS, OutDir = SSOutDir, fileName = "DNAme.csv")
  
  ## CLEAN SPACE
  rm(expdat_DNAme_merged, expdat_DNAme_fltr)
  gc()
  
  ## ---------- SNV data preprocessing
  
  if(if_SNV){
    
    brcds_SNV_SS = brcds_SS$brcds_SNV_SS[[SS]]
    
    ## Read experimental data files
    # expdat_SNV_list <- lapply(Prjcts, SNVPreprocessing, brcds_SNV = brcds_SNV_SS, InputDir = InputDir)
    ## Merge experimental data into only one SE object
    # expdat_SNV_merged <- do.call(rbind, expdat_SNV_list)
    
    ## Read and merge files
    expdat_SNV_merged <- SNVdataPreparation(Prjcts = Prjcts, brcds_SNV = brcds_SNV_SS, InputDir = InputDir)
    
    ## Filtering
    expdat_SNV_fltr <- SNVFiltering(expdat_SNV = expdat_SNV_merged, TopD = TopD)
    
    ## Save data
    expdat_SNV_sorted <- orderAndSaveData(expdat = expdat_SNV_fltr, smpls = smpls_SS, OutDir = SSOutDir, fileName = "SNV.csv")
    
    ## CLEAN SPACE
    rm(expdat_SNV_list, expdat_SNV_merged, expdat_SNV_fltr)
    gc()
  }
  
  ## ---------- SAVE METADATA
  print("Save metadata")
  
  if(if_SNV){
    expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted, "SNV" = expdat_SNV_sorted)
    brcds_list <- list("brcds_mRNA" = brcds_mRNA_SS, "brcds_miRNA" = brcds_miRNA_SS, "brcds_DNAme" = brcds_DNAme_SS, "brcds_SNV" = brcds_SNV_SS)
  }
  else{
    expdat_list <- list("mRNA" = expdat_mRNA_sorted, "miRNA" = expdat_miRNA_sorted, "DNAme" = expdat_DNAme_sorted)
    brcds_list <- list("brcds_mRNA" = brcds_mRNA_SS, "brcds_miRNA" = brcds_miRNA_SS, "brcds_DNAme" = brcds_DNAme_SS)
  }
  
  saveMetadata(OutDir = SSOutDir, Seed = Seed, smpls = smpls_SS, expdat_list = expdat_list, brcds_list = brcds_list, Prjcts = Prjcts)
  
  print("metadata SS")
  print(SSOutDir)
}

## --------------------------------------------------------------

