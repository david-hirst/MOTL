## FUNCTIONS

## ---------------------------------- TCGA RELATED FUNCTIONS ----

TCGATargetDataPrefiltering <- function(view, brcds_SS, SS, YTrgFull, Fctrzn){
  #' 
  #' Prepare target subset data
  #' 
  #' Extract data for the given sample names (brcds_SS)
  #' Remove features with variance equal to zero
  #' 
  #'
  #' @param view current view name data
  #' @param brcds_SS list of subset sample names for each subset number/view
  #' @param SS current subset number
  #' @param YTrgFull
  #' @param Fctrzn
  #' @returns YTrgSS
  
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
  #'
  #'
  #'
  #'
  
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

TargetDataPrefiltering <- function(view, YTrg_list, Fctrzn, smpls){
  #'
  #' Prepare the target data for a given view
  #' 
  #' Remove the features with variance equal to zero
  #' Harmonize features between the target data and the learning data 
  #' (keep only commun features)
  #' Order columns based on give sample order (smpls)
  #'
  #' @param view current view data name
  #' @param YTrg_list named list of data
  #' @param Fctrzn
  #' @param smpls ordered vector of sample names
  #' @returns 
  
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
  #'
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
  #' @param view current data view name
  #' @param YTrg_list list of training data
  #' @param normalization FALSE, or "Lrn", or "Trg"
  #' @param expdat_meta_Lrn
  #' @param transformation boolean
  #' @returns
  
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
  #'
  #'
  #'
  #'
  
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
  #' Extract calculated geomeans from learning set for counts data
  #' Extract the factorized learning set weight intercepts (from MOFA)
  #' Extract the factorized learning set weights (from MOFA)
  #' Extract the factorized learning set squared weights (from MOFA)
  #' Extract the learning set Tau and log(Tau) (from MOFA)
  #' For each parameter, commun features (YTrgFtrs) with the YTrg are selected
  #' YTrg, Tau and TauLn matrices are transposed
  #' 
  #' @param YTrg named list of data (data should have the same order of the columns)
  #' @param views vector of data names
  #' @param expdat_meta_Lrn list of learning set metadata 
  #' @param Fctrzn 
  #' @param likelihoods
  #' @returns TL_param list of parameters 
  #' (YTrg, GeoMeans_Lrn, Fctrzn_Lrn_W0, Fctrzn_Lrn_WFctrzn_Lrn_WSq, Tau, TauLn)
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
  ## MT to check
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

## --------------------------------------------------------------


## ----------------------------- TRANSFER LEARNING FUNCTIONS ----

## MT - question for david, where is the poisson calculation ?
## MT - Zeta matrix ?
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
  #' @returns Zeta Zeta matrix for the current data view
  
  if (likelihoods[[view]]=="bernoulli"){
    Zeta = sqrt(E_ZWSq[[view]])
  } else{
    Zeta = E_ZE_W[[view]]
  }
  
  return(Zeta)
}

Tau_init <- function(viewsLrn, Fctrzn, InputModel){
  #' 
  #' Initialization of the Tau values for each view
  #' 
  #' @param viewsLrn
  #' @param Fctrzn
  #' @param InputModel
  #' @returns 
  
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

## MT - Tau matrix ?
Tau_calculation <- function(view, likelihoods, Zeta, Tau){
  #'
  #' Update Tau values which only change for bernoulli data
  #' 
  #' @param view current view data name
  #' @param likelihoods list of data types
  #' @param Zeta Zeta matrix for the current view
  #' @param Tau 
  #' @returns Tau the (updated) Tau matrix for the current view data
  
  if (likelihoods[[view]]=="bernoulli"){
    Tau = (1/2)*(1/Zeta[[view]])*tanh(Zeta[[view]]/2)
  } else {
    Tau = Tau[[view]]
  }
  return(Tau)
}

TauLn_calculation <- function(view, likelihoodsLrn, Fctrzn, LrnFctrnDir, LrnSimple = TRUE){
  #'
  #' @param view
  #' @param likelihoodsLrn
  #' @param LrnSimple
  #' @param Fctrzn
  #' @param LrnFctrnDir
  #' @returns
  
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
  #' @param Zeta
  #' @param Tau
  #' @param CenterTrg
  #' @returns YGauss
  
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
  #' @param Tau
  #' @param Fctrzn_Lrn_WSq 
  #' @returns ZVar_m the calculated Z variances matrix for the current data
  
  ZVar_m = Tau[[view]] %*% Fctrzn_Lrn_WSq[[view]]
  return(ZVar_m)
}

ZMu_calculation <- function(view, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss){
  #'
  #' Z mu calculation for the current data
  #' 
  #' @param view current data view name
  #' @param k
  #' @param Fctrzn_Lrn_W
  #' @param Fctrzn_Lrn_W0
  #' @param Tau
  #' @param ZMu_0
  #' @param ZMu
  #' @param YGauss
  #' @returns 
  
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
  #' @param view
  #' @param likelihoods
  #' @param Tau
  #' @param TauLn
  #' @param E_ZWSq
  #' @param E_ZE_W
  #' @param Zeta
  #' @param YTrg
  #' @param YGauss
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
  #' @param view
  #' @param ZMu_0
  #' @param ZMu
  #' @param Fctrzn_Lrn_W0
  #' @param Fctrzn_Lrn_W
  #' @returns 
  
  E_ZE_W = cbind(ZMu_0,ZMu) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]]))
  return(E_ZE_W)
}

E_Z_SqE_W_Sq_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W){
  #' 
  #' Calculate 
  #'
  #' @param view
  #' @param ZMu_0
  #' @param ZMu
  #' @param Fctrzn_Lrn_W0
  #' @param Fctrzn_Lrn_W
  #' @returns 
  
  E_Z_SqE_W_Sq = (cbind(ZMu_0,ZMu)^2) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]])^2)
  return(E_Z_SqE_W_Sq)
}

E_ZSqE_WSq_update <- function(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq){
  #' 
  #' Calculate 
  #'
  #' @param view
  #' @param ZMu_0
  #' @param ZMuSq
  #' @param Fctrzn_Lrn_W0
  #' @param Fctrzn_Lrn_WSq
  #' @returns 
  
  E_ZSqE_WSq = cbind(ZMu_0^2,ZMuSq) %*% t(cbind(Fctrzn_Lrn_W0[[view]]^2,Fctrzn_Lrn_WSq[[view]]))
  return(E_ZSqE_WSq)
}

E_ZWSq_update <- function(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq){
  #' 
  #' Calculate 
  #'
  #' @param view
  #' @param E_ZE_W
  #' @param ZMuSq
  #' @param E_Z_SqE_W_Sq
  #' @param E_ZSqE_WSq
  #' @returns 
  
  E_ZWSq = (E_ZE_W[[view]]^2) - E_Z_SqE_W_Sq[[view]] + E_ZSqE_WSq[[view]]
  return(E_ZWSq)
}

WSq_calculation <- function(view, Fctrzn, LrnFctrnDir, LrnSimple = TRUE){
  #'
  #' load or calculate E[W^2] values
  #' factors were ordered in the same way as for other latent variables
  #' if any factors are dropped due to being inactive, they are at the end of the dataset
  #' so can filter based on dimension of W
  #' 
  #' @param view
  #' @param LrnSimple
  #' @param Fctrzn
  #' @param LrnFctrnDir
  #' @returns 
  
  if (LrnSimple){
    WSq = (Fctrzn@expectations[["W"]][[view]])^2
  } else {
    WSq = read.csv(file.path(LrnFctrnDir,paste0("WSq_",view,".csv")),header=FALSE)
    WSq = as.matrix(WSq)[,1:dim(Fctrzn@expectations[["W"]][[view]])[2]]
    rownames(WSq)=rownames(Fctrzn@expectations[["W"]][[view]])
  }
  return(WSq)
}

W0_calculation <- function(view, CenterTrg, Fctrzn, LrnFctrnDir){
  #'
  #' 
  #' @param view
  #' @param CenterTrg
  #' @param Fctrzn
  #' @param LrnFctrnDir
  #' @returns 
  
  if (CenterTrg){
    W0 = Fctrzn@expectations[["W"]][[view]][,1]*0
  } else {
    EstInts = readRDS(file.path(LrnFctrnDir,"EstimatedIntercepts.rds"))
    EstInts = EstInts$Intercepts
    W0 = EstInts[[view]]
  }
  return(W0)
}

transferLearning_function <- function(TL_param, MaxIterations, MinIterations, minFactors, 
                                      StartDropFactor, FreqDropFactor, StartELBO, FreqELBO, DropFactorTH, ConvergenceIts, ConvergenceTH,
                                      CenterTrg, PoisRateCstnt = 0.0001, outputDir = "./"){
  #'
  #' Transfer Learning with Variational Inference
  #'
  #' @param TL_param
  #' @param MaxIterations
  #' @param MinIterations
  #' @param minFactors
  #' @param StartDropFactor
  #' @param FreqDropFactor
  #' @param StartELBO
  #' @param FreqELBO
  #' @param DropFactorTH
  #' @param ConvergenceIts
  #' @param ConvergenceTH
  #' @param PoisRateCstnt
  #' @param outputDir
  #' @returns
  
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
    print(It)
    
    ## Drop factors explaining variance below threshold using MOFA formula:  
    # SS = np.square(Y[m][gg, :]).sum()
    # Res = np.sum((Y[m][gg, :] - Ypred) ** 2.0)
    # r2[g][m, k] = 1.0 - Res / SS as per paper
    BegK = dim(Fctrzn_Lrn_W[[1]])[2]
    if ((BegK > minFactors) & (It > 1) & (It > StartDropFactor) & (((It-StartDropFactor-1) %% FreqDropFactor) == 0)){
      
      print("Drop factors")
      
      SS_tmp <- sapply(views, function(view, YGauss, ZMu_0, Fctrzn_Lrn_W0){
        SS_tmp <- sum((YGauss[[view]] - (matrix(ZMu_0,ncol = 1) %*% t(Fctrzn_Lrn_W0[[view]])))^2, na.rm=TRUE)
        return(SS_tmp)
      }, YGauss, ZMu_0, Fctrzn_Lrn_W0, simplify = FALSE)
      
      #SS_tmp = lapply(split(SS_tmp, names(SS_tmp)), unname)
      
      VarExpl = sapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
        factorNames = colnames(ZMu)
        var_expl_tmp = sapply(factorNames, function(factorName, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
          RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,factorName]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,factorName]))))^2, na.rm=TRUE)
          var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
          return(var_expl_tmp)
          }, view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
        return(var_expl_tmp)
        }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp)
      var_expl_max <- apply(VarExpl,1,max)
      print(head(VarExpl))
      
      # var_expl_max <- numeric()
      # for (k in 1:dim(ZMu)[2]){
      #   var_expl_tmp <- unlist(lapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
      #     RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,k]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,k]))))^2, na.rm=TRUE)
      #     var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
      #     return(var_expl_tmp)
      #   }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp))
      #   var_expl_max <- c(var_expl_max, max(var_expl_tmp))
      # }
      
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
  
  # add names where necessary
  # rownames(ZMu) = smpls
  # colnames(ZMu) = colnames(Fctrzn_Lrn_W[[1]])
  
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
    'ss_start_time' = ss_start_time,
    'ss_fit_start_time' = ss_fit_start_time,
    'ss_end_time' = ss_end_time
  )
  saveRDS(TL_data, file.path(outputDir,"TL_data.rds"))
  
  ## delete before next subset in case i've missed something in the loop
  rm(list = c("YTrgSS", "YGauss", "ZMu_0", "ZMu", "Fctrzn_Lrn_W0", "Fctrzn_Lrn_W", "ELBO", "ELBO_P", "ELBO_Q", "ELBO_L", 
              'E_ZE_W', 'E_Z_SqE_W_Sq', 'E_ZSqE_WSq', 'E_ZWSq', "ZMuSq", "Fctrzn_Lrn_WSq", "Tau", "TauLn", "ZVar", "Zeta"))
  
  invisible(gc())
  
  TL_output <- list(
    "PoisRateCstnt" = PoisRateCstnt,
    "VarExpl" = VarExpl
  )
}


## --------------------------------------------------------------

## MT to remove
initTransferParameters_oldVersion <- function(views, brcds_SS, SS, YTrgFull, Fctrzn, likelihoods, expdat_meta_Lrn){
  #'
  #'
  #' @param view view name of the current data
  #' @param brcds_SS
  #' @param YTrgFull 
  #' @param Fctrzn
  #' @param likelihoods
  #' 
  
  ## Prepare YTrgSS
  print("Prepare YTrgSS")
  YTrgSS <- sapply(views, YTrgSS_preparation, brcds_SS = brcds_SS, SS = SS, YTrgFull = YTrgFull, Fctrzn = Fctrzn)
  FtrsCommon <- lapply(YTrgSS, rownames)
  
  ## GeoMeans
  print("Geomeans of learning set")
  GeoMeans_Lrn <- sapply(views, function(view, expdat_meta_Lrn, FtrsCommon){
    if (is.element(view,c('mRNA', 'miRNA'))){
      GeoMeans_Lrn = expdat_meta_Lrn[[paste0('GeoMeans_',view)]]
      FtrsKeep = is.element(names(GeoMeans_Lrn),FtrsCommon[[view]])
      GeoMeans_Lrn = GeoMeans_Lrn[FtrsKeep]
      GeoMeans_Lrn = GeoMeans_Lrn[match(FtrsCommon[[view]],names(GeoMeans_Lrn))]
    } else {
      GeoMeans_Lrn = numeric()
    }
    return(GeoMeans_Lrn)
  }, expdat_meta_Lrn, FtrsCommon)
  
  ## FACTORIZED LEARNING WEIGHTS MATRIX ZERO
  print("Factorized learning set weight intercepts")
  Fctrzn_Lrn_W0 <- sapply(views, function(view, Fctrzn, FtrsCommon){
    Fctrzn_Lrn_W0 = Fctrzn@expectations[["W0"]][[view]]
    FtrsKeep = is.element(names(Fctrzn_Lrn_W0),FtrsCommon[[view]])
    Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0[FtrsKeep]
    Fctrzn_Lrn_W0 = Fctrzn_Lrn_W0[match(FtrsCommon[[view]],names(Fctrzn_Lrn_W0))]
    return(Fctrzn_Lrn_W0)
  }, Fctrzn, FtrsCommon)
  
  ## FACTORIZED LEARNING WEIGHTS MATRIX
  print("Factorized learning set weights")
  Fctrzn_Lrn_W <- sapply(views, function(view, Fctrzn, FtrsCommon){
    Fctrzn_Lrn_W = Fctrzn@expectations[["W"]][[view]]
    FtrsKeep = is.element(rownames(Fctrzn_Lrn_W),FtrsCommon[[view]])
    Fctrzn_Lrn_W = Fctrzn_Lrn_W[FtrsKeep,]
    Fctrzn_Lrn_W = Fctrzn_Lrn_W[match(FtrsCommon[[view]],rownames(Fctrzn_Lrn_W)),]
    return(Fctrzn_Lrn_W)
  }, Fctrzn, FtrsCommon)
  
  ## FACTORIZED LEARNING WEIGHTS MATRIX SQUARED ?
  print("Factorized learning set squared weights")
  Fctrzn_Lrn_WSq <- sapply(views, function(view, Fctrzn, FtrsCommon){
    Fctrzn_Lrn_WSq = Fctrzn@expectations[["WSq"]][[view]]
    FtrsKeep = is.element(rownames(Fctrzn_Lrn_WSq),FtrsCommon[[view]])
    Fctrzn_Lrn_WSq = Fctrzn_Lrn_WSq[FtrsKeep,]
    Fctrzn_Lrn_WSq = Fctrzn_Lrn_WSq[match(FtrsCommon[[view]],rownames(Fctrzn_Lrn_WSq)),]
    return(Fctrzn_Lrn_WSq)
  }, Fctrzn, FtrsCommon)
  
  ## TAU PARAMETER
  print("Tau")
  Tau <- sapply(views, function(view, Fctrzn, YTrgSS, FtrsCommon){
    Tau = Fctrzn@expectations[["Tau"]][[view]]$group0
    FtrsKeep = is.element(rownames(Tau),FtrsCommon[[view]])
    Tau = Tau[FtrsKeep,]
    Tau = Tau[match(FtrsCommon[[view]],rownames(Tau)),]
    Tau = matrix(rowMeans(Tau, na.rm = TRUE),
                 nrow = dim(YTrgSS[[view]])[1], 
                 ncol = dim(YTrgSS[[view]])[2], 
                 byrow = FALSE)
    rownames(Tau) = rownames(YTrgSS[[view]])
    return(Tau)
  }, Fctrzn, YTrgSS, FtrsCommon)
  
  
  ## LOG TAU PARAMETER
  print("LOG Tau")
  TauLn <- sapply(views, function(view, likelihoods, Fctrzn, YTrgSS, FtrsCommon){
    if (likelihoods[view]=="gaussian"){
      TauLn = Fctrzn@expectations[["TauLn"]][[view]]
      FtrsKeep = is.element(names(TauLn),FtrsCommon[[view]])
      TauLn = TauLn[FtrsKeep]
      TauLn = TauLn[match(FtrsCommon[[view]],names(TauLn))]
      TauLn = matrix(TauLn,
                     nrow = dim(YTrgSS[[view]])[1], 
                     ncol = dim(YTrgSS[[view]])[2], 
                     byrow = FALSE)
      rownames(TauLn) = rownames(YTrgSS[[view]])
    } else{
      TauLn = numeric()
    }
    return(TauLn)
  }, likelihoods, Fctrzn, YTrgSS, FtrsCommon)
  
  ## PREPROCESS DATA
  YTrgSS <- sapply(views, function(view, YTrgSS, GeoMeans_Lrn){
    YTrgSS <- preprocessTransfer(view = view, YTrgSS = YTrgSS[[view]], GeoMeans_Lrn = GeoMeans_Lrn[[view]])
    return(YTrgSS)
  }, YTrgSS, GeoMeans_Lrn)
  
  ## transpose matrices where necessary - to make them samples x features
  YTrgSS <- sapply(YTrgSS, t)
  Tau = sapply(Tau, t)
  TauLn = sapply(TauLn, t)
  
  ## return
  TL_param = list(
    "YTrgSS" = YTrgSS,
    "GeoMeans_Lrn" = GeoMeans_Lrn,
    "Fctrzn_Lrn_W0" = Fctrzn_Lrn_W0,
    "Fctrzn_Lrn_W" = Fctrzn_Lrn_W,
    "Fctrzn_Lrn_WSq" = Fctrzn_Lrn_WSq,
    "Tau" = Tau,
    "TauLn" = TauLn
  )
  return(TL_param)
}

## old name: preprocessTransfer
## MT to remove also
preprocessTCGAData4TransferLearning <- function(view, YTrgSS_list, GeoMeans_Lrn, smpls){
  #'
  #' Prepare TCGA data for the transfert learning (for one given view)
  #' 
  #' Counts data are normalized and transformed (mRNA and miRNA)
  #' DNAme SE object is transformed into dataframe
  #' SNV data are already well shaped
  #' Columns data (samples) are then reshaped and sorted according the smpls vector
  #'
  #' @param YTrgSS list of data of one target subset
  #' @param view view name of the current data
  #' @param GeoMeans_Lrn list of geomeans calculated for learning data
  #' @param smpls vector of sample names
  #' @returns YTrgSS data well shaped for a specific view
  
  ## preprocess and transform
  if (is.element(view,c('mRNA', 'miRNA'))){
    YTrgSS = countsNormalization(expdat = YTrgSS_list[[view]], GeoMeans = GeoMeans_Lrn)
    YTrgSS = countsTransformation(expdat_count = YTrgSS$counts, TopD = nrow(YTrgSS$counts))
  } else if (view=='DNAme'){
    YTrgSS = assay(YTrgSS_list[[view]])
  }
  
  ## order columns
  colnames(YTrgSS) = substr(colnames(YTrgSS),1,16)
  YTrgSS = YTrgSS[,match(smpls,colnames(YTrgSS))]
  
  ## Return
  return(YTrgSS)
}

### MT to remove
TCGAinitTransferLearningParamaters <- function(YTrgFull, views, brcds_SS, SS, expdat_meta_Lrn, Fctrzn, likelihoods){
  #'
  #' Transfer Learning parameters initialization
  #' 
  #' 
  #' @param YTrgFull named list of data
  #' @param views vector of data names
  #' 
  #' 
  
  ## Prepare YTrgSS
  print("Prepare YTrgSS")
  YTrgSS <- sapply(views, TCGATrgSubsetPreparation, brcds_SS = brcds_SS, SS = SS, YTrgFull = YTrgFull, Fctrzn = Fctrzn)
  YTrgFtrs <- lapply(YTrgSS, rownames)
  
  ## PREPROCESS DATA
  YTrgSS <- sapply(views, function(view, YTrgSS, GeoMeans_Lrn){
    YTrgSS <- preprocessTCGAData4TransferLearning(view = view, YTrgSS = YTrgSS[[view]], GeoMeans_Lrn = GeoMeans_Lrn[[view]])
    return(YTrgSS)
  }, YTrgSS, GeoMeans_Lrn)
  
  TL_param = initTransferLearningParamaters(YTrgSS, views, expdat_meta_Lrn, Fctrzn, likelihoods)
  
  return(TL_param)
}
