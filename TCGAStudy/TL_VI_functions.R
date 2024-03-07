## FUNCTIONS

preprocessTransfer <- function(YTrgSS, view, GeoMeans_Lrn){
  #'
  #'
  #' @param YTrgSS
  #' @param views
  #' @param GeoMeans_Lrn
  #'
  
  ## preprocess and transform
  if (is.element(view,c('mRNA', 'miRNA'))){
    YTrgSS = countsNormalization(expdat = YTrgSS, GeoMeans = GeoMeans_Lrn)
    YTrgSS = countsTransformation(expdat_count = YTrgSS$counts, TopD = nrow(YTrgSS$counts))
  } else if (view=='DNAme'){
    YTrgSS = assay(YTrgSS)
  }
  
  ## order columns
  colnames(YTrgSS) = substr(colnames(YTrgSS),1,16)
  YTrgSS = YTrgSS[,match(smpls,colnames(YTrgSS))]
  
  ## Return
  return(YTrgSS)
}

YTrgSS_preparation <- function(view, brcds_SS, SS, YTrgFull, Fctrzn){
  #'
  #' 
  #'
  
  print(view)
  
  ## select samples and subset the Ydat
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

initTransferParameters <- function(views, brcds_SS, SS, YTrgFull, Fctrzn, likelihoods, expdat_meta_Lrn){
  #'
  #'
  #' @param view
  #' @param brcds_SS
  #' @param YTrgFull 
  #' @param Fctrzn
  #' @param likelihoods
  #' 
  
  ## INIT PARAMETER
  # YTrgSS = vector("list")
  # GeoMeans_Lrn = vector("list")
  # Fctrzn_Lrn_W0 = vector("list")
  # Fctrzn_Lrn_W = vector("list")
  # Fctrzn_Lrn_WSq = vector("list")
  # Tau = vector("list")
  # TauLn = vector("list")
  
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
  
  
  ## LEARNING TAU PARAMETER
  print("Learning set Tau")
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

Zeta_calculation <- function(view, likelihoods, E_ZWSq, E_ZE_W){
  #' Zeta values used for non-gaussian data
  #' for poisson Zeta_nd = E[\sum_{k} z_{n,k} w_{d,k}] so Zeta = ZMu %*% t(W)
  #' for bernoulli Zeta_nd = sqrt(E[(\sum_{k} z_{n,k} w_{d,k})^2])
  #' 
  
  if (likelihoods[[view]]=="bernoulli"){
    Zeta = sqrt(E_ZWSq[[view]])
  } else{
    Zeta = E_ZE_W[[view]]
  }
  
  return(Zeta)
}

Tau_calculation <- function(view, likelihoods, Zeta, Tau){
  #'Update Tau values which only change for bernoulli data
  if (likelihoods[[view]]=="bernoulli"){
    Tau = (1/2)*(1/Zeta[[view]])*tanh(Zeta[[view]]/2)
  } else {
    Tau = Tau[[view]]
  }
  return(Tau)
}

YGauss_calculation <- function(view, likelihoods, YTrgSS, Zeta, Tau, CenterTrg){
  #'
  #' Initialise / update Pseudo Y values - called YGauss here
  #' for gaussian data this is just the (centered) Y values which are fixed
  #' for non gaussian these are transformed y values that change after each update of z
  #' the y pseudo values are centered at each step if the centering option is selected
  #' for gaussian data this is done for It>=0, for others it is It>0
  
  if (likelihoods[[view]]=="poisson"){
    YGauss = Zeta[[view]] - plogis(Zeta[[view]])*(1-YTrgSS[[view]]/(log(1+exp(Zeta[[view]])) + PoisRateCstnt))/Tau[[view]]
  } else if (likelihoods[[view]]=="bernoulli"){
    YGauss = (2*YTrgSS[[view]] - 1) / (2*Tau[[view]])
  } else {
    YGauss = YTrgSS[[view]]
  }
  
  if (CenterTrg){
    YGauss = sweep(YGauss,2,as.vector(colMeans(YGauss, na.rm = TRUE)),"-")
  }
  
  return(YGauss)
}

ZVar_calculation <- function(view, Tau, Fctrzn_Lrn_WSq){
  #' Z variances using initialised / updated tau values and W^2 values
  #' based on the appendix of the mofa paper and Github code
  ZVar_m = Tau[[view]] %*% Fctrzn_Lrn_WSq[[view]]
  return(ZVar_m)
}

ZMu_calculation <- function(view, k, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss){
  #'
  #'
  #'
  ZMu_tmp1 = matrix(Fctrzn_Lrn_W[[view]][,k], nrow=dim(Tau[[view]])[1],ncol=dim(Tau[[view]])[2],byrow = TRUE)
  ZMu_tmp1 = Tau[[view]] * ZMu_tmp1
  ZMu_tmp2 = cbind(ZMu_0,ZMu[,-k]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,-k]))
  ZMu_tmp2 = YGauss[[view]] - ZMu_tmp2
  ZMu_tmp3 = ZMu_tmp1 * ZMu_tmp2
  ZMu_tmp3 = rowSums(ZMu_tmp3, na.rm = TRUE)
  return(ZMu_tmp3)
}

ELBO_calculation <- function(view, likelihoods, Tau, TauLn, E_ZWSq, E_ZE_W, Zeta, YTrgSS, YGauss){
  #'
  #'
  #' likelihoods
  #' for poisson and bernoulli it is the bound which is used
  #' for gaussian it is expanded gaussian log likelihood
  #' it seems in MOFA they do not allow for the centering that is done for the 
  #' pseudo data for non gaussian data
  #' they used the raw uncentered data for the elbo - this seems strange
  #' although i guess it acts as a lower bound for the lower bound ...
  #'
  #'
  
  if(likelihoods[[view]]=="poisson"){
    # b_nd is an upper bound for -log(p(y|x))
    # b_nd = k_nd/2 * (x_nd - zeta_nd)^2 + (x_nd - zeta_nd)f'(zeta_nd) + f(zeta_nd)
    # f'(a) = (1/(1+e^(-a)))(1 - y/log(1+e^a))
    # f(a) = log(1+e^a) - ylog(log(1+e^a))
    # The elbo component is -b_nd as this is a lower bound for log(p(y|x))
    ## A CONSTANT IS ADDED TO RATE CALCULATIONS HERE AS PER MOFA CODE TO AVOID ERRORS
    
    ELBO_L_tmpA = 0.5 * Tau[[view]] * (E_ZWSq[[view]] - 2 * E_ZE_W[[view]] * Zeta[[view]] + Zeta[[view]]^2)
    ELBO_L_tmpB = (E_ZE_W[[view]] - Zeta[[view]]) * plogis(Zeta[[view]]) * (1 - YTrgSS[[view]]/(log(1+exp(Zeta[[view]]))+ PoisRateCstnt))
    ELBO_L_tmpC = (log(1+exp(Zeta[[view]]))+ PoisRateCstnt) - YTrgSS[[view]] * log((log(1+exp(Zeta[[view]]))+ PoisRateCstnt))
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
    ELBO_L_tmpB = 0.5 * ((2 * YTrgSS[[view]] - 1) * E_ZE_W[[view]] - Zeta[[view]])
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
  E_ZE_W = cbind(ZMu_0,ZMu) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]]))
  return(E_ZE_W)
}

E_Z_SqE_W_Sq_update <- function(view, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W){
  E_Z_SqE_W_Sq = (cbind(ZMu_0,ZMu)^2) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]])^2)
  return(E_Z_SqE_W_Sq)
}

E_ZSqE_WSq_update <- function(view, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq){
  E_ZSqE_WSq = cbind(ZMu_0^2,ZMuSq) %*% t(cbind(Fctrzn_Lrn_W0[[view]]^2,Fctrzn_Lrn_WSq[[view]]))
  return(E_ZSqE_WSq)
}

E_ZWSq_update <- function(view, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq){
  E_ZWSq = (E_ZE_W[[view]]^2) - E_Z_SqE_W_Sq[[view]] + E_ZSqE_WSq[[view]]
  return(E_ZWSq)
}

transferLearning_function <- function(TL_param, MaxIterations, MinIterations, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, minFactors, 
                                      StartDropFactor, FreqDropFactor, StartELBO, FreqELBO, DropFactorTH, ConvergenceIts, ConvergenceTH, 
                                      PoisRateCstnt){
  #'
  #' Transfer Learning with Variational Inference
  #'
  #' @param MaxIterations
  #' @param Fctrzn_Lrn_W
  #' @param Fctrzn_Lrn_W0
  #' @param M
  #' @param YGauss
  #' @param ZMu
  #' 
  #' 
  #
  
  ss_fit_start_time = Sys.time()
  
  ## RETREIVE PARAMETERS
  YTrgSS <- TL_param$YTrgSS
  GeoMeans_Lrn <- TL_param$GeoMeans_Lrn
  Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
  Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
  Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
  Tau <- TL_param$Tau
  TauLn <- TL_param$TauLn
  
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
      }, YGauss, ZMu_0, Fctrzn_Lrn_W0)
      
      var_expl_max <- numeric()
      for (k in 1:dim(ZMu)[2]){
        var_expl_tmp <- unlist(lapply(views, function(view, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp){
          RSS_tmp = sum((YGauss[[view]] - (cbind(ZMu_0,ZMu[,k]) %*% t(cbind(Fctrzn_Lrn_W0[[view]],Fctrzn_Lrn_W[[view]][,k]))))^2, na.rm=TRUE)
          var_expl_tmp = 1-(RSS_tmp/SS_tmp[[view]])
          return(var_expl_tmp)
        }, YGauss, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W, SS_tmp))
        var_expl_max <- c(var_expl_max, max(var_expl_tmp))
      }
      
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
    
    print("Exepected")
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
  rownames(ZMu) = smpls
  colnames(ZMu) = colnames(Fctrzn_Lrn_W[[1]])
  
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
  saveRDS(TL_data, file.path(TL_SSOutDir,"TL_data.rds"))
  
  ## delete before next subset in case i've missed something in the loop
  rm(list = c("YTrgSS", "YGauss", "ZMu_0", "ZMu", "Fctrzn_Lrn_W0", "Fctrzn_Lrn_W", "ELBO", "ELBO_P", "ELBO_Q", "ELBO_L", 
              'E_ZE_W', 'E_Z_SqE_W_Sq', 'E_ZSqE_WSq', 'E_ZWSq', "ZMuSq", "GeoMeans_Lrn", "Fctrzn_Lrn_WSq", "Tau", "TauLn", "ZVar", "Zeta"))
  
  invisible(gc())
}

TauLn_calculation <- function(view, likelihoodsLrn, LrnSimple, Fctrzn, LrnFctrnDir){
  #'
  #'
  #'
  
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

WSq_calculation <- function(view, LrnSimple, Fctrzn, LrnFctrnDir){
  #' load or calculate E[W^2] values
  #' factors were ordered in the same way as for other latent variables
  #' if any factors are dropped due to being inactive, they are at the end of the dataset
  #' so can filter based on dimension of W
  
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
  #'
  
  if (CenterTrg){
    W0 = Fctrzn@expectations[["W"]][[view]][,1]*0
  } else {
    EstInts = readRDS(file.path(LrnFctrnDir,"EstimatedIntercepts.rds"))
    EstInts = EstInts$Intercepts
    W0 = EstInts[[view]]
  }
  return(W0)
}

Tau_init <- function(viewsLrn, Fctrzn, InputModel){
  Tau <- h5read(InputModel, "expectations/Tau")
  Tau <- Tau[match(viewsLrn,names(Tau))]
  
  for (i in 1:length(viewsLrn)){
    view = viewsLrn[i]
    rownames(Tau[[view]]$group0)=rownames(Fctrzn@expectations[["W"]][[view]])
  }
  return(Tau)
}


