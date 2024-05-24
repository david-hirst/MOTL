## FUNCTIONS

## ---------------------------------- TCGA RELATED FUNCTIONS ----

TCGATargetDataPrefiltering <- function(view, brcds_SS, SS, YTrgFull, Fctrzn){
  #' 
  #' Filter out TCGA target subset data according variance
  #' 
  #' Extract data for the given sample names (brcds_SS)
  #' Remove features with variance equal to zero
  #' Match features between training set and learning set
  #' 
  #'
  #' @param view current view name data
  #' @param brcds_SS list of subset sample names for each subset number/view
  #' @param SS current subset number
  #' @param YTrgFull named list of training set data
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @returns the subset data for current view and SS number
  
  print(view)
  
  ## select samples and subset the YTrg
  brcds = brcds_SS[[paste0('brcds_',view,'_SS')]][[SS]]
  SmplsKeep = is.element(colnames(YTrgFull[[view]]), brcds$brcds)
  YTrgSS = YTrgFull[[view]][,SmplsKeep]
  print(paste0("YTrgSS dimensions: ", dim(YTrgSS)))
  
  ## prefiltering, only condition is that variance >0
  if(is.element(view,c('mRNA', 'DNAme', 'miRNA'))){
    FtrsKeep = rowVars(assay(YTrgSS), na.rm = TRUE)>0
  } else {
    FtrsKeep = rowVars(YTrgSS, na.rm = TRUE)>0
  }
  FtrsKeep[is.na(FtrsKeep)]=FALSE
  YTrgSS = YTrgSS[FtrsKeep,]
  print(paste0("YTrgSS dimensions after prefiltering: ", dim(YTrgSS)))
  
  ## harmonize features between Trg SS and Lrn data
  FtrsLrn = Fctrzn@features_metadata$feature[Fctrzn@features_metadata$view==view]
  FtrsCommon = FtrsLrn[is.element(FtrsLrn,rownames(YTrgSS))]
  FtrsKeep = is.element(rownames(YTrgSS),FtrsCommon)
  YTrgSS = YTrgSS[FtrsKeep,]
  YTrgSS = YTrgSS[match(FtrsCommon,rownames(YTrgSS)),]
  
  return(YTrgSS)
}

TCGATargetDataPreparation <- function(views, YTrgFull, brcds_SS, SS, Fctrzn, smpls, normalization = "Lrn", expdat_meta_Lrn, transformation = TRUE){
  #'
  #' Prepare TCGA data for transfer learning
  #' 
  #' Filter out features according variance
  #' Reshape data into matrices
  #' Order samples to have the same columns order between different data
  #' Normalize and/or transform counts data (DESeq2 normalization)
  #' If normalization = FALSE: no normalization
  #' If normalization = "Lrn": normalize with learning set geomeans
  #' If normalization = "Trg": normalize without geomeans
  #' 
  #' @param views list of data names
  #' @param YTrgFull named list of training data
  #' @param brcds_SS list of subset sample names for each subset number/view
  #' @param SS current subset number
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param smpls sample names
  #' @param normalization FALSE or "Lrn" or "Trg"
  #' @param expdat_meta_Lrn list of learning set factorization metadata
  #' @param transformation FALSE or TRUE
  #' @returns list of prepared subset data for the current subset number
  
  print("Feature prefiltering")
  
  ## Feature variance prefiltering and feature harmonization
  YTrgSSFull = sapply(views, function(view, brcds_SS, SS, YTrgFull, Fctrzn){
    YTrgSS = TCGATargetDataPrefiltering(view, brcds_SS, SS, YTrgFull, Fctrzn)
    return(YTrgSS)
  }, brcds_SS, SS, YTrgFull, Fctrzn)
  
  ## Reshape data
  YTrgSSFull = sapply(views, function(view, YTrgSSFull, smpls){
    YTrgSS = YTrgSSFull[[view]]
    if(is.element(view,c('mRNA', 'miRNA', 'DNAme'))){
      YTrgSS = assay(YTrgSS)
    }
    
    ## order columns
    colnames(YTrgSS) = substr(colnames(YTrgSS),1,16)
    YTrgSS = YTrgSS[,match(smpls,colnames(YTrgSS))]
    
    return(YTrgSS)
  }, YTrgSSFull, smpls)
  
  ## Normalization and transformation
  print("Normalization and transformation")
  YTrgSSFull = sapply(views, preprocessCountsData, YTrgSSFull, normalization, expdat_meta_Lrn, transformation)
  
  ## Return
  return(YTrgSSFull)
}

## --------------------------------------------------------------

## ----------------------------------------- INIT PARAMETERS ----

## Non-TCGA Target DATA

mRNA_addVersion <- function(expdat, Lrndat){
  #'
  #' get ensembl mRNA versions from learning dataset
  #' and attach to the correspnding mRNA ensembl id in the target dataset
  #'
  #' @param expdat the mRNA matrix from the target dataset with genes in rows
  #' @param Lrndat the mRNA W matrix from the learning dataset factorization with genes in rows
  #' @returns the target mRNA matrix with versions attached
  
  tmp = as.data.frame(do.call(rbind,strsplit(rownames(Lrndat),"[.]")))
  # match to stripped ids from target set
  tmp = data.frame(V1 = rownames(expdat)) %>%
    dplyr::left_join(tmp, by = c('V1')) %>%
    as.data.frame()
  # rename target dataset features
  rownames(expdat) = paste0(tmp$V1,'.',tmp$V2)
 
  # return tidied up matrix
  return(expdat)
  
}

TargetDataPrefiltering <- function(view, YTrg_list, Fctrzn, smpls){
  #'
  #' Prepare the target data for a given view
  #' 
  #' Remove the features with variance equal to zero
  #' Harmonize features between the target data and the learning data 
  #' (keep only shared features)
  #' Order columns based on give sample order (smpls)
  #'
  #' @param view current view data name
  #' @param YTrg_list named list of data
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param smpls ordered vector of sample names
  #' @returns prepared data for the current view
  
  YTrg = YTrg_list[[view]]
  
  ## prefiltering, only condition is that variance >0
  FtrsKeep = rowVars(YTrg, na.rm = TRUE)>0
  FtrsKeep[is.na(FtrsKeep)]=FALSE
  YTrg = YTrg[FtrsKeep,]
  print(paste0("YTrg dimensions after prefiltering: ", dim(YTrg)))
  
  ## harmonize features between Trg and Lrn data
  FtrsLrn = Fctrzn@features_metadata$feature[Fctrzn@features_metadata$view==view]
  FtrsCommon = FtrsLrn[is.element(FtrsLrn,rownames(YTrg))]
  FtrsKeep = is.element(rownames(YTrg),FtrsCommon)
  YTrg = YTrg[FtrsKeep,]
  YTrg = YTrg[match(FtrsCommon,rownames(YTrg)),]
  
  ## order columns
  YTrg = YTrg[,match(smpls,colnames(YTrg))]
  
  return(YTrg)
  
}

GeoMeans_Lrn_init <- function(view, expdat_meta_Lrn, YTrgFtrs){
  #'
  #' Retrieve the calculated geomeans of the learning set
  #' 
  #' @param view current view data name
  #' @param expdat_meta_Lrn list of learning set factorization metadata
  #' @param YTrgFtrs feature names of the current view
  #' @returns calculated geomeans of the learning set
  
  if (is.element(view,c('mRNA', 'miRNA'))){
    GeoMeans_Lrn = expdat_meta_Lrn[[paste0('GeoMeans_',view)]]
    FtrsKeep = is.element(names(GeoMeans_Lrn),YTrgFtrs)
    GeoMeans_Lrn = GeoMeans_Lrn[FtrsKeep]
    GeoMeans_Lrn = GeoMeans_Lrn[match(YTrgFtrs,names(GeoMeans_Lrn))]
  } else {
    GeoMeans_Lrn = numeric()
  }
  return(GeoMeans_Lrn)
  
}

preprocessCountsData <- function(view, YTrg_list, normalization = FALSE, expdat_meta_Lrn, transformation = FALSE){
  #'
  #' Normalize and/or transform counts data
  #' 
  #' If normalization = FALSE: no normalization
  #' If normalization = "Lrn": normalize with learning set geomeans
  #' If normalization = "Trg": normalize without geomeans
  #'
  #' @param view current data view name
  #' @param YTrg_list named list of training data
  #' @param normalization FALSE or "Lrn" or "Trg" 
  #' @param expdat_meta_Lrn list of learning set factorization metadata
  #' @param transformation boolean
  #' @returns transformed/normalized counts data for the current view
  
  ## Select current data
  YTrg = YTrg_list[[view]]
  
  if(view %in% c("mRNA", "miRNA")){
    
    print(view)
    
    ## Normalization
    ## if Lrn, retreive the learning set geomeans and normalize with it
    ## if Trg, normalize without geomeans
    if(is.character(normalization)){
      if(normalization == "Lrn"){
        print("Normalize with the Learning set GeoMeans")
        GeoMeans <- GeoMeans_Lrn_init(view, expdat_meta_Lrn, rownames(YTrg))
      }
      if(normalization == "Trg"){
        print("Normalize without GeoMeans")
        GeoMeans = normalization
      }
      YTrg = SummarizedExperiment(assays = list(YTrg))
      YTrg = countsNormalization(expdat = YTrg, GeoMeans = GeoMeans)
      YTrg = YTrg$counts
    }else{print("No normalization")}
    
    ## Transformation
    if(transformation){
      print("Log transform data")
      YTrg = countsTransformation(expdat_count = YTrg, TopD = nrow(YTrg))
    }else{print("No transformation")}
  }
  
  ## Return
  return(YTrg)
}

TargetDataPreparation <- function(views, YTrg_list, Fctrzn, smpls, expdat_meta_Lrn, normalization = FALSE, transformation = FALSE){
  #'
  #' Prepare data for the transfer learning
  #' 
  #' Filter out features according variance
  #' Normalize/transform counts data
  #'
  #' @param views list of data names
  #' @param YTrg_list named list of training data
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param smpls sample names
  #' @param normalization FALSE or "Lrn" or "Trg"
  #' @param expdat_meta_Lrn list of learning set factorization metadata
  #' @param transformation FALSE or TRUE
  #' @returns list of prepared subset data for the current subset number
  
  ## Feature variance prefiltering and feature harmonization
  print("Feature prefiltering")
  YTrg = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  
  ## Normalization and transformation
  print("Normalization and transformation")
  YTrg = sapply(views, preprocessCountsData, YTrg, normalization, expdat_meta_Lrn, transformation)
  
  ## Return 
  return(YTrg)
}

initTransferLearningParamaters <- function(YTrg, views, expdat_meta_Lrn, Fctrzn, likelihoods){
  #'
  #' Transfer Learning parameters initialization
  #' 
  #' Extract the factorized learning set weight intercepts (from MOFA)
  #' Extract the factorized learning set weights (from MOFA)
  #' Extract the factorized learning set squared weights (from MOFA)
  #' Extract the learning set Tau and log(Tau) (from MOFA)
  #' For each parameter, commun features (YTrgFtrs) with the YTrg are selected
  #' YTrg, Tau and TauLn matrices are transposed
  #' 
  #' @param YTrg named list of data (data should have the same order of the columns)
  #' @param views vector of data names
  #' @param expdat_meta_Lrn list of learning set factorization metadata
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param likelihoods list of data types
  #' @returns list of initialized parameters for transfer learning
  #' (YTrg, Fctrzn_Lrn_W0, Fctrzn_Lrn_WFctrzn_Lrn_WSq, Tau, TauLn)
  #' 
  
  ## Feature names in each data
  YTrgFtrs <- lapply(YTrg, rownames)
  
  ## FACTORIZED LEARNING WEIGHTS MATRIX ZERO
  print("Factorized learning set weight intercepts")
  Fctrzn_Lrn_W0 <- sapply(views, function(view, Fctrzn, YTrgFtrs){
    Fctrzn_Lrn_W0 = Fctrzn@expectations[["W0"]][[view]]
    FtrsKeep = is.element(names(Fctrzn_Lrn_W0),YTrgFtrs[[view]])
    Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0[FtrsKeep]
    Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0[match(YTrgFtrs[[view]],names(Fctrzn_Lrn_W0))]
    return(Fctrzn_Lrn_W0)
  }, Fctrzn, YTrgFtrs)
  
  ## FACTORIZED LEARNING WEIGHTS MATRIX
  print("Factorized learning set weights")
  Fctrzn_Lrn_W <- sapply(views, function(view, Fctrzn, YTrgFtrs){
    Fctrzn_Lrn_W = Fctrzn@expectations[["W"]][[view]]
    FtrsKeep = is.element(rownames(Fctrzn_Lrn_W),YTrgFtrs[[view]])
    Fctrzn_Lrn_W = Fctrzn_Lrn_W[FtrsKeep,]
    Fctrzn_Lrn_W = Fctrzn_Lrn_W[match(YTrgFtrs[[view]],rownames(Fctrzn_Lrn_W)),]
    return(Fctrzn_Lrn_W)
  }, Fctrzn, YTrgFtrs)
  
  ## FACTORIZED LEARNING WEIGHTS MATRIX SQUARED
  print("Factorized learning set squared weights")
  Fctrzn_Lrn_WSq <- sapply(views, function(view, Fctrzn, YTrgFtrs){
    Fctrzn_Lrn_WSq = Fctrzn@expectations[["WSq"]][[view]]
    FtrsKeep = is.element(rownames(Fctrzn_Lrn_WSq),YTrgFtrs[[view]])
    Fctrzn_Lrn_WSq = Fctrzn_Lrn_WSq[FtrsKeep,]
    Fctrzn_Lrn_WSq = Fctrzn_Lrn_WSq[match(YTrgFtrs[[view]],rownames(Fctrzn_Lrn_WSq)),]
    return(Fctrzn_Lrn_WSq)
  }, Fctrzn, YTrgFtrs)
  
  ## TAU PARAMETER
  print("Tau")
  Tau <- sapply(views, function(view, Fctrzn, YTrg, YTrgFtrs){
    Tau = Fctrzn@expectations[["Tau"]][[view]]$group0
    FtrsKeep = is.element(rownames(Tau),YTrgFtrs[[view]])
    Tau = Tau[FtrsKeep,]
    Tau = Tau[match(YTrgFtrs[[view]],rownames(Tau)),]
    Tau = matrix(rowMeans(Tau, na.rm = TRUE),
                 nrow = dim(YTrg[[view]])[1], 
                 ncol = dim(YTrg[[view]])[2], 
                 byrow = FALSE)
    rownames(Tau) = rownames(YTrg[[view]])
    return(Tau)
  }, Fctrzn, YTrg, YTrgFtrs)
  
  ## LOG TAU PARAMETER
  print("LOG Tau")
  TauLn <- sapply(views, function(view, likelihoods, Fctrzn, YTrg, YTrgFtrs){
    if (likelihoods[view]=="gaussian"){
      TauLn = Fctrzn@expectations[["TauLn"]][[view]]
      FtrsKeep = is.element(names(TauLn),YTrgFtrs[[view]])
      TauLn = TauLn[FtrsKeep]
      TauLn = TauLn[match(YTrgFtrs[[view]],names(TauLn))]
      TauLn = matrix(TauLn,
                     nrow = dim(YTrg[[view]])[1], 
                     ncol = dim(YTrg[[view]])[2], 
                     byrow = FALSE)
      rownames(TauLn) = rownames(YTrg[[view]])
    } else{
      TauLn = numeric()
    }
    return(TauLn)
  }, likelihoods, Fctrzn, YTrg, YTrgFtrs)
  
  ## transpose matrices where necessary - to make them samples x features
  YTrg <- lapply(YTrg, t)
  Tau = lapply(Tau, t)
  TauLn = lapply(TauLn, t)
  
  ## return
  TL_param = list(
    "YTrg" = YTrg,
    "Fctrzn_Lrn_W0" = Fctrzn_Lrn_W0,
    "Fctrzn_Lrn_W" = Fctrzn_Lrn_W,
    "Fctrzn_Lrn_WSq" = Fctrzn_Lrn_WSq,
    "Tau" = Tau,
    "TauLn" = TauLn
  )
  return(TL_param)
}

## LEARNING DATA

Tau_init <- function(viewsLrn, Fctrzn, InputModel){
  #' 
  #' Initialization of the Tau values for each view
  #' 
  #' @param viewsLrn list of view in learning set
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param InputModel factorization model object of learning set (from MOFA)
  #' @returns named list of Tau matrices
  
  ## Extract Tau from the factorization of the learning set
  Tau <- h5read(InputModel, "expectations/Tau")
  Tau <- Tau[match(viewsLrn,names(Tau))]
  
  ## For each view, transfer rownames into the corresponding Tau matrix
  for (i in 1:length(viewsLrn)){
    view = viewsLrn[i]
    rownames(Tau[[view]]$group0)=rownames(Fctrzn@expectations[["W"]][[view]])
  }
  
  # Return a named list of Tau matrix
  return(Tau)
}

## MT - maybe change name into TauLn_init ?
## MT - view = viewsLrn ?
TauLn_calculation <- function(view, likelihoodsLrn, Fctrzn, LrnFctrnDir, LrnSimple = TRUE){
  #'
  #' Initialization of the log(Tau) values for each learning view
  #'
  #' @param view current learning view name
  #' @param likelihoodsLrn type of data in the learning set
  #' @param LrnSimple if TRUE then E[W^2] and E[LnTau] are calculated from E[W] and E[Tau] otherwise imported
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param LrnFctrnDir directory where log(Tau) values are saved
  #' @returns log(Tau) matrix for the current view
  
  if (likelihoodsLrn[view]=="gaussian"){
    if (LrnSimple){
      TauLn = log(Fctrzn@expectations[["Tau"]][[view]]$group0[,1])
    } else {
      TauLn = as.vector(read.csv(file.path(LrnFctrnDir,paste0("TauLn_",view,".csv")), header=FALSE)$V1)
      names(TauLn)=rownames(Fctrzn@expectations[["W"]][[view]])
    }
  }else {
    TauLn = numeric()
  }
  return(TauLn)
}

## MT - maybe change name into WSq_init ?
## MT - explain LrnSimple input
WSq_calculation <- function(view, Fctrzn, LrnFctrnDir, LrnSimple = TRUE){
  #'
  #' load or calculate E[W^2] values
  #' factors were ordered in the same way as for other latent variables
  #' if any factors are dropped due to being inactive, they are at the end of the dataset
  #' so can filter based on dimension of W
  #' 
  #' @param view current view name
  #' @param LrnSimple 
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param LrnFctrnDir directory where WSq values are saved
  #' @returns weights squared matrix for the current view
  
  if (LrnSimple){
    WSq = (Fctrzn@expectations[["W"]][[view]])^2
  } else {
    WSq = read.csv(file.path(LrnFctrnDir,paste0("WSq_",view,".csv")),header=FALSE)
    WSq = as.matrix(WSq)[,1:dim(Fctrzn@expectations[["W"]][[view]])[2]]
    rownames(WSq)=rownames(Fctrzn@expectations[["W"]][[view]])
  }
  return(WSq)
}

## MT - maybe change name into W0_init ?
W0_calculation <- function(view, CenterTrg, Fctrzn, LrnFctrnDir){
  #'
  #' Initialization of the weight intercept
  #' 
  #' @param view current view data name
  #' @param CenterTrg use (FALSE) or not (TRUE) use estimated intercepts
  #' @param Fctrzn learning factorization model object (from MOFA)
  #' @param LrnFctrnDir directory where estimated intercept values are saved
  #' @returns weight intercepts values of tjhe current data
  
  if (CenterTrg){
    W0 = Fctrzn@expectations[["W"]][[view]][,1]*0
  } else {
    EstInts = readRDS(file.path(LrnFctrnDir,"EstimatedIntercepts.rds"))
    EstInts = EstInts$Intercepts
    W0 = EstInts[[view]]
  }
  return(W0)
}

intercepts_calculation <- function(seed, expdat_meta, Fctrzn, FctrznDir, LrnDir){
  #'
  #'
  #' @param seed
  #' @param expdat_meta learning set metadata
  #' @param Fctrzn
  #' @param FctrznDir
  #' @param LrnDir
  
  print("Estimation of the intercept")
  
  fit_start_time = Sys.time()
  
  ## Init seed
  mode(Seed) = 'integer'
  set.seed(Seed)
  
  ## Extract data from factorization model object
  views = Fctrzn@data_options$views
  names(views) <- views
  likelihoods = Fctrzn@model_options$likelihoods
  M = Fctrzn@dimensions$M
  D = Fctrzn@dimensions$D
  
  # loop through the views and estimate the intercept
  # for gaussian data its just the mean
  # for other data will try mle with a naive estimator as backup
  
  intercepts_list <- lapply(views, function(view, likelihoods, D, YTmp, expdat_meta, Fctrzn){
    
    print(view)
    
    # likelihood = likelihoods[which(names(likelihoods)==view)]
    # DTmp = D[which(names(D)==view)]
    likelihood = likelihoods[view]
    DTmp = D[view]
    
    # YTmp = read.table(file = file.path(ExpDataDir, paste0(view,'.csv')), sep = ",")
    YTmp <- as.data.frame(data.table::fread(file = file.path(LrnDir, paste0(view,'.csv')), sep = ","))
    YTmp = t(as.matrix(YTmp))
    rownames(YTmp) = expdat_meta$smpls
    colnames(YTmp) = expdat_meta[[which(names(expdat_meta) == paste0("ftrs_",view))]] 
    
    ZWTmp = Fctrzn@expectations$Z$group0 %*% 
      t(Fctrzn@expectations$W[[which(names(Fctrzn@expectations$W)==view)]])
    
    # mean(colnames(YTmp)==colnames(ZWTmp))
    # mean(rownames(YTmp)==rownames(ZWTmp))
    
    invisible(gc())
    
    if(likelihood=="gaussian"){
      InterceptsNaive = colMeans(YTmp, na.rm = TRUE)
      Intercepts = InterceptsNaive
      InterceptsMethod = rep("Naive",length(Intercepts))
      names(InterceptsMethod) = names(InterceptsNaive)
    } else if(likelihood == "poisson"){
      
      ## naive intercept based on approximation to feature means of ZW
      InterceptsNaive = as.vector(log(-1 + exp(colMeans(YTmp))))
      
      ## mle estimate of intercept for ZW
      ## if optimiser fails for a feature will return the naive estimate
      
      Intercepts_df <- do.call(rbind, lapply(c(1:DTmp), function(d, YTmp, ZWTmp){
        ## compute for each feature vector
        YTmp_d = YTmp[,d]
        YTmp_d_keep = !is.na(YTmp_d)
        YTmp_d = YTmp_d[YTmp_d_keep]
        
        ZWTmp_d = ZWTmp[YTmp_d_keep,d]
        
        ## NLL function to optimize
        # nLL = function(interceptMLE) -sum(stats::dpois(YLrn[,d], log(1 + exp(ZWLrn[,d]+interceptMLE)), log = TRUE))
        nLL = function(interceptMLE) -sum(log(
          stats::dpois(YTmp_d, log(1 + exp(ZWTmp_d+interceptMLE)))[stats::dpois(YTmp_d[,d], log(1 + exp(ZWTmp_d+interceptMLE)))!=0]
        ))
        
        ## try to solve it and use the result otherwise use the naive estimate
        # interceptMLEfit = try(as.vector(stats4::mle(nLL, start=list(interceptMLE=0))@coef[1]))
        interceptMLEfit = try(as.vector(stats4::mle(nLL, start=list(interceptMLE=InterceptsNaive[d]))@coef[1]))
        
        if (class(interceptMLEfit)=="try-error"){
          InterceptsTmp = InterceptsNaive[d]
          InterceptsMethodTmp = "Naive"
        } else {
          InterceptsTmp = interceptMLEfit
          InterceptsMethodTmp = "MLE"
        }
        
        intercept <- data.frame("intercept" = InterceptsTmp, "Method" = InterceptsMethodTmp, row.names = names(InterceptsNaive)[d])
        
        return(intercept)
      }, YTmp, ZWTmp))
      
      Intercepts = setNames(Intercepts_df$intercept, row.names(Intercepts_df))
      InterceptsMethod = setNames(Intercepts_df$Method, row.names(Intercepts_df))
      
    } else if(likelihood=="bernoulli"){
      
      ## naive intercept based on approximation to feature means of ZW
      InterceptsNaive = log(colMeans(YTmp, na.rm = TRUE)/(1-colMeans(YTmp, na.rm = TRUE)))
      
      ## mle estimate of intercept for ZW
      ## if optimiser fails for a feature will return the naive estimate
      
      # DTmp = 10
      
      Intercepts_df <- do.call(rbind, lapply(c(1:DTmp), function(d, YTmp, ZWTmp){
        ## compute for each feature vector
        
        YTmp_d = YTmp[,d]
        YTmp_d_keep = !is.na(YTmp_d)
        YTmp_d = YTmp_d[YTmp_d_keep]
        
        ZWTmp_d = ZWTmp[YTmp_d_keep,d]
        
        ## NLL function to optimize
        nLL = function(InterceptMLE) -sum(log(
          dbinom(YTmp_d, size=1, plogis(ZWTmp_d+InterceptMLE))[dbinom(YTmp_d, size=1, plogis(ZWTmp_d+InterceptMLE))!=0]
        ))
        
        ## try to solve it and use the result otherwise use the naive estimate
        
        interceptMLEfit = try(stats4::mle(nLL, start=list(InterceptMLE=InterceptsNaive[d]))@coef[1])
        
        if (class(interceptMLEfit)=="try-error"){
          InterceptsTmp = InterceptsNaive[d]
          InterceptsMethodTmp = "Naive"
        } else {
          InterceptsTmp = interceptMLEfit
          InterceptsMethodTmp = "MLE"
        }
        intercept <- data.frame("intercept" = InterceptsTmp, "Method" = InterceptsMethodTmp, row.names = names(InterceptsNaive)[d])
        
        return(intercept)
        
      }, YTmp, ZWTmp))
      
      Intercepts = setNames(Intercepts_df$intercept, row.names(Intercepts_df))
      InterceptsMethod = setNames(Intercepts_df$Method, row.names(Intercepts_df))
      
    } else{
      InterceptsNaive = numeric()
      Intercepts = numeric()
      InterceptsMethod = character()
    }
    
    intercepts_list <- list(
      "InterceptsNaive" = InterceptsNaive,
      "Intercepts" = Intercepts,
      "InterceptsMethod" = InterceptsMethod)
    
    return(intercepts_list)
    
  }, likelihoods, D, YTmp, expdat_meta, Fctrzn)
  
  ## save the intercepts in the relevant factorization folder
  
  fit_end_time = Sys.time()
  
  EstimatedIntercepts <- list(
    "Seed" = Seed,
    "InterceptsNaive" = lapply(views, function(v, intercepts_list){return(intercepts_list[[v]][["InterceptsNaive"]])}, intercepts_list), 
    "Intercepts" = lapply(views, function(v, intercepts_list){return(intercepts_list[[v]][["Intercepts"]])}, intercepts_list),
    "InterceptsMethod" = lapply(views, function(v, intercepts_list){return(intercepts_list[[v]][["InterceptsMethod"]])}, intercepts_list),
    "fit_start_time" = fit_start_time,
    "fit_end_time" = fit_end_time
  )
  
  saveRDS(EstimatedIntercepts,file.path(FctrznDir,"EstimatedIntercepts.rds"))
  
  print("finished")
  
  
  
}

## --------------------------------------------------------------


## ----------------------------- TRANSFER LEARNING FUNCTIONS ----

Zeta_calculation <- function(view, likelihoods, E_ZWSq, E_ZE_W){
  #'
  #' Calculate the Zeta matrix for the current data view
  #' 
  #' Zeta values used for non-gaussian data
  #' for poisson Zeta_nd = E[\sum_{k} z_{n,k} w_{d,k}] so Zeta = ZMu %*% t(W)
  #' for bernoulli Zeta_nd = sqrt(E[(\sum_{k} z_{n,k} w_{d,k})^2])
  #' 
  #' @param view current data view name
  #' @param likelihoods list of data types
  #' @param E_ZWSq 
  #' @param E_ZE_W
  #' @returns Zeta matrix for the current data view
  
  if (likelihoods[[view]]=="bernoulli"){
    Zeta = sqrt(E_ZWSq[[view]])
  } else{
    Zeta = E_ZE_W[[view]]
  }
  
  return(Zeta)
}

Tau_calculation <- function(view, likelihoods, Zeta, Tau){
  #'
  #' Update Tau values which only change for bernoulli data
  #' 
  #' @param view current view data name
  #' @param likelihoods list of data types
  #' @param Zeta Zeta matrix for the current view
  #' @param Tau list of tau matrices
  #' @returns (updated) Tau matrix for the current view data
  
  if (likelihoods[[view]]=="bernoulli"){
    Tau = (1/2)*(1/Zeta[[view]])*tanh(Zeta[[view]]/2)
  } else {
    Tau = Tau[[view]]
  }
  return(Tau)
}

YGauss_calculation <- function(view, likelihoods, YTrg, Zeta, Tau, CenterTrg){
  #' 
  #' Initialize or update pseudo Y values (YGauss)
  #'
  #' For gaussian data this is just the (centered) Y values which are fixed
  #' For non gaussian these are transformed y values that change after each update of z
  #' the y pseudo values are centered at each step if the centering option is selected
  #' For gaussian data this is done for It>=0, for others it is It>0
  #' 
  #' @param view current view name
  #' @param likelihoods list of data types
  #' @param YTrg current data matrix
  #' @param Zeta list of Zeta matrices
  #' @param Tau list of Tau matrices
  #' @param CenterTrg use (FALSE) or not (TRUE) use estimated intercepts
  #' @returns pseudo Y values for the current view
  
  if (likelihoods[[view]]=="poisson"){
    YGauss = Zeta[[view]] - plogis(Zeta[[view]])*(1-YTrg[[view]]/(log(1+exp(Zeta[[view]])) + PoisRateCstnt))/Tau[[view]]
  } else if (likelihoods[[view]]=="bernoulli"){
    YGauss = (2*YTrg[[view]] - 1) / (2*Tau[[view]])
  } else {
    YGauss = YTrg[[view]]
  }
  
  if (CenterTrg){
    YGauss = sweep(YGauss,2,as.vector(colMeans(YGauss, na.rm = TRUE)),"-")
  }
  
  return(YGauss)
}

ZVar_calculation <- function(view, Tau, Fctrzn_Lrn_WSq){
  #' 
  #' Calculation of the Z variances for the current data
  #' 
  #' Z variances using initialised / updated tau values and W^2 values
  #' based on the appendix of the mofa paper and Github code
  #' 
  #' @param view current data view name
  #' @param Tau list of Tau matrices
  #' @param Fctrzn_Lrn_WSq Factorized learning set squared weights
  #' @returns calculated Z variances matrix for the current data
  
  ZVar_m = Tau[[view]] %*% Fctrzn_Lrn_WSq[[view]]
  return(ZVar_m)
}

ZMu_calculation <- function(view, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss){
  #'
  #' Z mu calculation for the current data
  #' 
  #' @param view current data view name
  #' @param k feature index in the current data
  #' @param Fctrzn_Lrn_W list of Factorized learning set weight matrices
  #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
  #' @param Tau list of Tau matrices
  #' @param ZMu_0 list of ZMu intercepts matrices
  #' @param ZMu list of ZMu matrices
  #' @param YGauss list of pseudo Y value matrices
  #' @returns ZMu values for the current view
  
  ZMu_tmp1 = matrix(Fctrzn_Lrn_W[[view]][,k], nrow=dim(Tau[[view]])[1],ncol=dim(Tau[[view]])[2],byrow = TRUE)
  ZMu_tmp1 = Tau[[view]] * ZMu_tmp1
  ZMu_tmp2 = cbind(ZMu_0,ZMu[,-k]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,-k]))
  ZMu_tmp2 = YGauss[[view]] - ZMu_tmp2
  ZMu_tmp3 = ZMu_tmp1 * ZMu_tmp2
  ZMu_tmp3 = rowSums(ZMu_tmp3, na.rm = TRUE)
  
  return(ZMu_tmp3)
}

ELBO_calculation <- function(view, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrg, YGauss){
  #'
  #' Calculate the ELBO value for the current view/iterations
  #'
  #' likelihoods
  #' for poisson and bernoulli it is the bound which is used
  #' for gaussian it is expanded gaussian log likelihood
  #' it seems in MOFA they do not allow for the centering that is done for the 
  #' pseudo data for non gaussian data
  #' they used the raw uncentered data for the elbo - this seems strange
  #' although i guess it acts as a lower bound for the lower bound ...
  #'
  #' @param view current view dara name
  #' @param likelihoods list of data types
  #' @param Tau list of Tau matrices
  #' @param TauLn list of log(Tau) matrices
  #' @param E_ZWSq
  #' @param E_ZE_W 
  #' @param Zeta list of Zeta matrices
  #' @param YTrg list of data
  #' @param YGauss list of pseudo Y value matrices
  #' @returns 
  
  if(likelihoods[[view]]=="poisson"){
    # b_nd is an upper bound for -log(p(y|x))
    # b_nd = k_nd/2 * (x_nd - zeta_nd)^2 + (x_nd - zeta_nd)f'(zeta_nd) + f(zeta_nd)
    # f'(a) = (1/(1+e^(-a)))(1 - y/log(1+e^a))
    # f(a) = log(1+e^a) - ylog(log(1+e^a))
    # The elbo component is -b_nd as this is a lower bound for log(p(y|x))
    ## A CONSTANT IS ADDED TO RATE CALCULATIONS HERE AS PER MOFA CODE TO AVOID ERRORS
    
    ELBO_L_tmpA = 0.5 * Tau[[view]] * (E_ZWSq[[view]] - 2 * E_ZE_W[[view]] * Zeta[[view]] + Zeta[[view]]^2)
    ELBO_L_tmpB = (E_ZE_W[[view]] - Zeta[[view]]) * plogis(Zeta[[view]]) * (1 - YTrg[[view]]/(log(1+exp(Zeta[[view]]))+ PoisRateCstnt))
    ELBO_L_tmpC = (log(1+exp(Zeta[[view]]))+ PoisRateCstnt) - YTrg[[view]] * log((log(1+exp(Zeta[[view]]))+ PoisRateCstnt))
    ELBO_L_tmp = -sum(ELBO_L_tmpA + ELBO_L_tmpB + ELBO_L_tmpC, na.rm = TRUE)
    
  } else if (likelihoods[[view]]=="bernoulli"){
    # Basedthe MOFA paper, MOFA code and eq(7) from jakoola paper
    # if g(a) = 1/(1 + e^(-a)) is the logistic (sigmoid) function
    # h_nd = (2 * y_nd - 1) * x_nd
    # lambda(a) = tanh(a/2)/(4*a)
    # b_nd = log(g(zeta_nd)) + (h_nd - zeta_nd)/2 - lambda(zeta_nd)(h_nd^2 - zeta_nd^2)
    # as tau_nd = 2 * lambda(zeta_nd) this becomes
    # b_nd = log(g(zeta_nd)) + (h_nd - zeta_nd)/2 - tau_nd/2 * (x_nd^2 - zeta_nd^2)
    # here b_nd is the lower bound for log(p(y|x))
    
    ELBO_L_tmpA = log(plogis(Zeta[[view]]))
    ELBO_L_tmpB = 0.5 * ((2 * YTrg[[view]] - 1) * E_ZE_W[[view]] - Zeta[[view]])
    ELBO_L_tmpC = 0.5 * Tau[[view]] * (E_ZWSq[[view]] - Zeta[[view]]^2)
    ELBO_L_tmp = sum(ELBO_L_tmpA + ELBO_L_tmpB - ELBO_L_tmpC, na.rm = TRUE) 
    
  } else {
    # gaussian log likelihood
    # log(f(y_nd|x_nd,tau_nd)) = 1/2 * (log(tau_nd) - log(2*pi) - tau_nd * (y_nd - x_nd)^2)
    ELBO_L_tmpA = TauLn[[view]] - log(2 * pi)
    ELBO_L_tmpB = Tau[[view]] * (YGauss[[view]]^2 - 2 * YGauss[[view]] * E_ZE_W[[view]] + E_ZWSq[[view]])
    ELBO_L_tmp = sum(0.5 * (ELBO_L_tmpA - ELBO_L_tmpB), na.rm = TRUE)
  }
  return(ELBO_L_tmp)
}

E_ZE_W_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W){
  #' 
  #' Calculate 
  #'
  #' @param view current view name
  #' @param ZMu_0 list of ZMu intercept matrices
  #' @param ZMu list of ZMu matrices
  #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
  #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
  #' @returns 
  
  E_ZE_W = cbind(ZMu_0,ZMu) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]]))
  return(E_ZE_W)
}

E_Z_SqE_W_Sq_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W){
  #' 
  #' Calculate 
  #'
  #' @param view current view name
  #' @param ZMu_0 list of ZMu intercept matrices
  #' @param ZMu list of ZMu matrices
  #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
  #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
  #' @returns 
  
  E_Z_SqE_W_Sq = (cbind(ZMu_0,ZMu)^2) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]])^2)
  return(E_Z_SqE_W_Sq)
}

E_ZSqE_WSq_update <- function(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq){
  #' 
  #' Calculate 
  #'
  #' @param view current view name
  #' @param ZMu_0 list of ZMu intercept matrices
  #' @param ZMuSq list of ZMu squared matrices
  #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
  #' @param Fctrzn_Lrn_WSq  list of factorized learning set weight squared matrices
  #' @returns 
  
  E_ZSqE_WSq = cbind(ZMu_0^2,ZMuSq) %*% t(cbind(Fctrzn_Lrn_W0[[view]]^2,Fctrzn_Lrn_WSq[[view]]))
  return(E_ZSqE_WSq)
}

E_ZWSq_update <- function(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq){
  #' 
  #' Calculate 
  #'
  #' @param view current view name
  #' @param E_ZE_W
  #' @param ZMuSq list of ZMu squared matrices
  #' @param E_Z_SqE_W_Sq
  #' @param E_ZSqE_WSq
  #' @returns 
  
  E_ZWSq = (E_ZE_W[[view]]^2) - E_Z_SqE_W_Sq[[view]] + E_ZSqE_WSq[[view]]
  return(E_ZWSq)
}

VarExplFun <- function(views, YGauss, ZMu_0, Fctrzn_Lrn_W0, ZMu, Fctrzn_Lrn_W){
  #' 
  #' Calculate the variance explained by each factor for eah view
  #' 
  #' @param views list of view data names
  #' @param YGauss list of pseudo Y value matrices
  #' @param ZMu_0 list of ZMu intercept matrices
  #' @param ZMu list of ZMu matrices
  #' @param Fctrzn_Lrn_W0 list of factorized learning set weight intercept matrices
  #' @param Fctrzn_Lrn_W list of factorized learning set weight matrices
  #' @returns variance explained matrix
  
  SS_tmp <- sapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0){
    SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0,ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm=TRUE)
    return(SS_tmp)
  }, YGauss, ZMu_0, Fctrzn_Lrn_W0, simplify = FALSE)
  
  VarExpl = sapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
    factorNames = colnames(ZMu)
    var_expl_tmp = sapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
      RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,factorName]))))^2, na.rm=TRUE)
      var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
      return(var_expl_tmp)
    }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
    return(var_expl_tmp)
  }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
  
  return(VarExpl)
}

transferLearning_function <- function(TL_param, MaxIterations, MinIterations, minFactors, 
                                      StartDropFactor, FreqDropFactor, StartELBO, FreqELBO, DropFactorTH, ConvergenceIts, ConvergenceTH,
                                      CenterTrg, PoisRateCstnt = 0.0001, ss_start_time = NULL, outputDir = "./"){
  #'
  #' Transfer Learning with Variational Inference
  #'
  #' @param TL_param list of initialized parameters for transfer learning
  #' @param MaxIterations maximum number of iteration
  #' @param MinIterations minimum number of iteration
  #' @param minFactors minimum number of factors
  #' @param StartDropFactor after which iteration to start dropping factors
  #' @param FreqDropFactor how often to drop factors
  #' @param StartELBO which iteration to start checking ELBO on
  #' @param FreqELBO how often to assess the ELBO
  #' @param DropFactorTH factor with lowest max variance, that is less than this, is dropped
  #' @param ConvergenceIts
  #' @param ConvergenceTH 
  #' @param CenterTrg Center Trg with own means or use estimated intercepts
  #' @param PoisRateCstnt amount to add to the poison rate function to avoid errors
  #' @param ss_start_time analyse start time (steps before transfer learning)
  #' @param outputDir output directory name
  #' @returns list of transfer learning data
  
  ss_fit_start_time = Sys.time()
  
  ## RETREIVE PARAMETERS
  YTrgSS <- TL_param$YTrg
  Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
  Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
  Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
  Tau <- TL_param$Tau
  TauLn <- TL_param$TauLn
  smpls = rownames(TL_param$YTrg[[1]])
  
  ## INIT PARAMETERS
  ELBO = numeric()
  convergence_token = 0
  
  ## FOR EACH ITERATION
  for (It in 0:MaxIterations){
    
    PrintMessage = paste0("It=",It)
    
    ## Drop factors explaining variance below threshold using MOFA formula:  
    # SS = np.square(Y[m][gg, :]).sum()
    # Res = np.sum((Y[m][gg, :] - Ypred) ** 2.0)
    # r2[g][m, k] = 1.0 - Res / SS as per paper
    BegK = dim(Fctrzn_Lrn_W[[1]])[2]
    if ((BegK > minFactors) & (It > 1) & (It > StartDropFactor) & (((It-StartDropFactor-1) %% FreqDropFactor) == 0)){
      
      print("Drop factors")
      
      VarExpl = VarExplFun(views = views, YGauss = YGauss, ZMu_0 = ZMu_0, Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0, ZMu = ZMu, Fctrzn_Lrn_W = Fctrzn_Lrn_W)
      
      # SS_tmp <- sapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0){
      #   SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0,ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm=TRUE)
      #   return(SS_tmp)
      # }, YGauss, ZMu_0, Fctrzn_Lrn_W0, simplify = FALSE)
      
      # VarExpl = sapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
      #   factorNames = colnames(ZMu)
      #   var_expl_tmp = sapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
      #     RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,factorName]))))^2, na.rm=TRUE)
      #     var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
      #     return(var_expl_tmp)
      #     }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
      #   return(var_expl_tmp)
      #   }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
      
      var_expl_max <- apply(VarExpl,1,max)
      
      ## drop factor with lowest max variance explained if below the threshold
      if (min(var_expl_max)<DropFactorTH){
        fctrs_to_keep = !base::rank(var_expl_max, ties.method = "first")==1
        ZMu = ZMu[,fctrs_to_keep]
        ZMuSq = ZMuSq[,fctrs_to_keep]
        for (i in 1:length(views)){
          Fctrzn_Lrn_W[[i]] = Fctrzn_Lrn_W[[i]][,fctrs_to_keep]
          Fctrzn_Lrn_WSq[[i]] = Fctrzn_Lrn_WSq[[i]][,fctrs_to_keep]
        }
        
        ## recalculate expectations based on new number of factors
        ## explanation of these formulae are further down the script
        E_ZE_W <- sapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        E_Z_SqE_W_Sq <- sapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
        E_ZSqE_WSq <- sapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
        E_ZWSq <- sapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
      }
    }
    
    if (It>0){
      print("Zeta, Tau and YGauss")
      
      ## Zeta values used for non-gaussian data
      Zeta <- sapply(views, Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)
      
      ## Update Tau values which only change for bernoulli data
      Tau <- sapply(views, Tau_calculation, likelihoods, Zeta, Tau)
      
      ## Initialise / update Pseudo Y values - called YGauss here
      YGauss <- sapply(views, YGauss_calculation, likelihoods, YTrgSS, Zeta, Tau, CenterTrg)
    }
    
    ## Z variances using initialised / updated tau values and W^2 values
    print("Zeta")
    ZVar <- Reduce("+", lapply(views, ZVar_calculation, Tau, Fctrzn_Lrn_WSq))
    ZVar = 1/(ZVar + 1)
    
    ## Initialize or update ZMu values
    if (It == 0){
      print("Init Z")
      # initialise with means of learning set Z
      ZMu = matrix(data=as.vector(colMeans(Fctrzn@expectations$Z$group0)),
                   nrow = dim(ZVar)[1], ncol = dim(ZVar)[2],
                   byrow = TRUE)
      rownames(ZMu) = smpls
      colnames(ZMu) = colnames(Fctrzn_Lrn_W[[1]])
      
      # vector of 1s to act as multiplier of W intercept term
      ZMu_0 = rep(1,dim(ZVar)[1])
    } else {
      print("update Z")
      # variational updates
      for (k in 1:dim(ZMu)[2]){
        ZMu_tmp <- Reduce("+", lapply(views, ZMu_calculation, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss))
        ZMu[,k] = ZMu_tmp*ZVar[,k] #update factor value for subsequent calculation
      }
    }
    
    PrintMessage = paste0(PrintMessage,' K=',dim(ZMu)[2])
    
    ## Z^2 values
    print("Z squared")
    ZMuSq = ZVar + ZMu^2
    
    ## Some pre calculations - results used in various parts: E_ZE_W and ZZWW
    ## E_ZE_W_nd = E[\sum_{k} z_{n,k} w_{d,k}]
    ## E_ZWSq_nd = E[(\sum_{k} z_{n,k} w_{d,k})^2] = E[\sum_k \sum_j z_{n,k} w_{d,k} z_{n,j} w_{d,j}]
    ## If A = square(ZMu%*%t(W)) - square(ZMu)%*%square(t(W))
    ## And B = ZMuSq%*%t(WSq)
    ## Then A_nd = \sum_{k} \sum_{j != k} E[z_{n,k}]E[w_{n,k}]E[z_{n,j}]E[w_{n,j}]
    ## And B_nd = \sum_{k} E[(z_{n,k})^2]E[(w_{d,k})^2]
    ## E_ZWSq_nd = (A + B)_nd = (square(ZMu%*%t(W)) - square(ZMu)%*%square(t(W)) + ZMuSq%*%t(WSq))_nd
    
    print("Expected")
    E_ZE_W <- sapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    E_Z_SqE_W_Sq <- sapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
    E_ZSqE_WSq <- sapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
    E_ZWSq <- sapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)    
    
    ## Calculate and check the ELBO
    
    if ((It > 0) & (It >= StartELBO) & (((It-StartELBO) %% FreqELBO) == 0)){
      
      print("ELBO")
      
      ELBO_L <- do.call(sum, lapply(X = views, FUN = ELBO_calculation, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrgSS, YGauss))
      
      ## The ELBO components for the variational and prior distributions for Z
      ELBO_P = sum(0.5 * (-log(2*pi) - ZMuSq))
      ELBO_Q = sum(0.5 * (-1 - log(2*pi) - log(ZVar)))
      ELBO = c(ELBO, ELBO_L + ELBO_P - ELBO_Q)
      
      PrintMessage = paste0(PrintMessage,' ELBO=', round(ELBO[length(ELBO)],2))
      
      # Calculate and check the change in ELBO
      # I originally didn't allow negative changes in ELBO before convergence but MOFA does so I now allow it
      if (length(ELBO)>=2){
        ELBO_delta = 100*abs((ELBO[length(ELBO)] - ELBO[length(ELBO)-1])/ELBO[1])
        if ((ELBO_delta < ConvergenceTH)){
          convergence_token = convergence_token + 1
        } else {
          convergence_token = 0
        }
        PrintMessage = paste0(PrintMessage,' ELBO_delta=', round(ELBO_delta,8), '%')
      }
      
    }
    
    print(PrintMessage)
    
    ## can the algorithm be stopped?
    if ((It >= 2) & (It >= MinIterations) & (convergence_token >= ConvergenceIts)){
      print("converged")
      break
    }
    
  }
  
  ## Variance explained calculation with final factors
  VarExpl = VarExplFun(views = views, YGauss = YGauss, ZMu_0 = ZMu_0, Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0, ZMu = ZMu, Fctrzn_Lrn_W = Fctrzn_Lrn_W)
  
  # export the data for further analysis
  ss_end_time = Sys.time()
  TL_data = list(
    'YTrgSS' = YTrgSS,
    'YGauss' = YGauss,
    'ZMu_0' = ZMu_0,
    'ZMu' = ZMu,
    'Fctrzn_Lrn_W0' = Fctrzn_Lrn_W0,
    'Fctrzn_Lrn_W' = Fctrzn_Lrn_W,
    'ELBO' = ELBO,
    'VarExpl' = VarExpl,
    'ss_start_time' = ss_start_time,
    'ss_fit_start_time' = ss_fit_start_time,
    'ss_end_time' = ss_end_time
  )
  saveRDS(TL_data, file.path(outputDir,"TL_data.rds"))
  
  ## delete before next subset in case i've missed something in the loop
  rm(list = c("YTrgSS", "YGauss", "ZMu_0", "ZMu", "Fctrzn_Lrn_W0", "Fctrzn_Lrn_W", "ELBO", "ELBO_P", "ELBO_Q", "ELBO_L", 
              'E_ZE_W', 'E_Z_SqE_W_Sq', 'E_ZSqE_WSq', 'E_ZWSq', "ZMuSq", "Fctrzn_Lrn_WSq", "Tau", "TauLn", "ZVar", "Zeta"))
  
  invisible(gc())
  
  return(TL_data)
  
}


## --------------------------------------------------------------
