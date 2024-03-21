
#############
##### Evaluation of factorizations
#############

## Packages

library(MOFA2)
library(rhdf5)
library(rjson)
library(dplyr)

#######
### options and functions
#######

Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)


Prjcts = c('LAML_PAAD', 'LAML_SKCM','PAAD_SKCM',
           'LAML_PAAD_SKCM')

TopD = 5000 ## how many features were retained

SS_size = c(5, 5, 5,
            5) ## how many samples per cancer type in each sub set
names(SS_size)=Prjcts

TrgFullK = c(100, 100, 100,
             100) # Starting K for the full trg set factorization
names(TrgFullK)=Prjcts

TrgFullTH = '01TH' #string with digits to right of decimal for drop threshold followed by TH

TrgSSK = c(10, 10, 10,
           15) # Starting K for the trg set direct factorization
names(TrgSSK) = Prjcts

FctrznMethods = c('Direct', 'TL_VI', 'Random_TL')
RandomKRange = 'same K as TL_VI' # hard coded can ignore
FeatWseScl_W = TRUE ## scale weight vectors featurewise? otherwise scaled factorwise

## function for scaling weight vectors and ZW vectors
ScaleW = function(x) {
  x/norm(matrix(x),type="F")
}

CorMethod = "pearson" ## either "pearson" or "spearman"

GetBestHit = function(x) {
  which(base::rank(-x, ties.method = "first")==1)
}

GetBestScore = function(x){
  max(x,na.rm = TRUE)
}

corCutOff = 0.05 ## pvalue theshold for including correlations
pvalueTH = 0.05 ## THRESHOLD FOR DIFFERENTIALLY ACTIVE factors

### 
## loop through projects to evaluate
###

# empty vectors to append scores to
prjct_name = character()
ss_iteration = integer()
fctrzn_method = character()
score_name = character()
score = numeric()
# empty list for storing factor matches
factor_matches = vector('list')

for (Prjct in Prjcts){
  
  print(Prjct)

  ### directories  
  TrgFullRootDir = file.path(paste0('Trg_',Prjct,'_Full_',TopD,'D'))
  TrgFullDir = file.path(TrgFullRootDir, paste0('Fctrzn_',TrgFullK[Prjct],'K_',TrgFullTH))
  TrgSSDir = paste0('Trg_',Prjct,'_SS',SS_size[Prjct],'_',TopD,'D')
  
  ### import Trg Full factorization 
  InputModel = file.path(TrgFullDir,"Model.hdf5")
  TrgFull_Fctrzn = load_model(file = InputModel)
  TrgFull_Views = TrgFull_Fctrzn@data_options$views
  
  ## Import TrgSS meta data  
  TrgSS_meta_data = readRDS(file.path(TrgSSDir,'expdat_meta.rds'))
  Trg_smpls = readRDS(file.path(TrgSSDir,'brcds_SS.rds'))
  if (regexpr('[_]',Prjct)>0){
    Trg_prjcts = Trg_smpls$prjcts_SS
  }
  Trg_smpls = Trg_smpls$smpls_SS
  
  ## if appropriate get significant factors from full factorization
  if (regexpr('[_]',Prjct)>0){
    # get Z matrix
    TrgFull_Z_all = TrgFull_Fctrzn@expectations$Z$group0
    # get project for each sample and put in sample order as in Z matrix
    TrgFull_meta_data = readRDS(file.path(TrgFullRootDir,'expdat_meta.rds'))    
    TrgFull_prjcts_all = TrgFull_meta_data$brcds_mRNA$prjct   
    names(TrgFull_prjcts_all) = substr(TrgFull_meta_data$brcds_mRNA$brcds,1,16)
    TrgFull_prjcts_all = as.factor(TrgFull_prjcts_all[rownames(TrgFull_Z_all)])
    
    # loop through factors to get pvalues using either wilcox rank sum or kruskal test
    TrgFull_Z_pv = numeric()
    for (fctr in 1:ncol(TrgFull_Z_all)){      
      if (length(levels(TrgFull_prjcts_all)) > 2){
        DifFactTest = kruskal.test(TrgFull_Z_all[,fctr] ~ TrgFull_prjcts_all)
      } else {
        DifFactTest = wilcox.test(TrgFull_Z_all[,fctr] ~ TrgFull_prjcts_all)
      }
      TrgFull_Z_pv = c(TrgFull_Z_pv, DifFactTest$p.value)     
    }
    TrgFull_Z_pv = p.adjust(TrgFull_Z_pv, method="BH")
    TrgFull_SigFct = which(TrgFull_Z_pv<pvalueTH)
    
    GT_positives = length(TrgFull_SigFct)  
  }
  
  ## loop through subsets and factorization methods
  
  for (ss in 1:TrgSS_meta_data$SS_count){
    
    SS = paste0('SS_',ss)
    print(SS)
    ## get scores from factorization of full trg set for ss samples
    TrgSS_smpls = Trg_smpls[[SS]]
    TrgFull_Z = TrgFull_Fctrzn@expectations$Z$group0
    TrgFull_smpls_kp = is.element(rownames(TrgFull_Z),TrgSS_smpls)
    TrgFull_Z = TrgFull_Z[TrgFull_smpls_kp,]
    TrgFull_Z = TrgFull_Z[match(TrgSS_smpls,rownames(TrgFull_Z)),]
    
    for (fm in FctrznMethods){
      
      print(fm)
    
      ## set up directory
      if (fm == 'Direct'){
        TrgSSFctrznDir = paste0('Fctrzn_',TrgSSK[Prjct],'K')
      } else if (fm == 'Random_TL'){
        TrgSSFctrznDir = 'TL_VI'
      } else {
        TrgSSFctrznDir = fm
      }
      TrgSSFctrznDir = file.path(TrgSSDir,SS,TrgSSFctrznDir)
      
      ## import factorization fit
      if (fm == 'Direct'){
        TrgSS_Fctrzn = load_model(file = file.path(TrgSSFctrznDir,"Model.hdf5"))
        TrgSS_Z = TrgSS_Fctrzn@expectations$Z$group0
        TrgSS_W_tmp = TrgSS_Fctrzn@expectations$W
  
      } else {
        TrgSS_Fctrzn = readRDS(file.path(TrgSSFctrznDir,"TL_data.rds"))
        TrgSS_Z = TrgSS_Fctrzn$ZMu
        TrgSS_W_tmp = TrgSS_Fctrzn$Fctrzn_Lrn_W
      }
      
      ## harmonize samples
      TrgSS_Z = TrgSS_Z[match(TrgSS_smpls,rownames(TrgSS_Z)),]
      
      ## randomize if required
      if (fm == 'Random_TL'){
        # RandomK = sample(RandomKRange,1)
        RandomK = dim(TrgSS_Z)[2]
        TrgSS_Z = matrix(
          rnorm(dim(TrgSS_Z)[1]*RandomK,0,1),
          nrow=dim(TrgSS_Z)[1], ncol=RandomK
          )
        rownames(TrgSS_Z)=TrgSS_smpls
      }
      
      ## harmonize views and features
      
      TrgSS_Views = names(TrgSS_W_tmp)
      Shared_Views = TrgSS_Views[is.element(TrgSS_Views,TrgFull_Views)]
      
      TrgFull_W = vector("list")
      TrgSS_W = vector("list")
      
      for (view in Shared_Views){
        
        TrgFull_W[[view]] = TrgFull_Fctrzn@expectations$W[[view]]
        TrgFull_Features = rownames(TrgFull_W[[view]])
        
        ## randomize if required
        if (fm == 'Random_TL'){
          
          # TrgSS_Features = TrgFull_Features
          # Shared_Features = TrgFull_Features
          
          TrgSS_Features = sample(TrgFull_Features,
                                  sum(is.element(rownames(TrgSS_W_tmp[[view]]),
                                                 TrgFull_Features)),
                                  replace = FALSE
                                  )
          Shared_Features = TrgSS_Features
          
          TrgSS_W[[view]] = matrix(
            rnorm(length(Shared_Features)*RandomK,0,1),
            nrow=length(Shared_Features), ncol=RandomK
          )
          
          rownames(TrgSS_W[[view]])=Shared_Features
          
        } else {
          
          TrgSS_W[[view]] = TrgSS_W_tmp[[view]]
          TrgSS_Features = rownames(TrgSS_W[[view]])
          Shared_Features = TrgSS_Features[is.element(TrgSS_Features,TrgFull_Features)]
          
          TrgSS_W[[view]] = TrgSS_W[[view]][is.element(TrgSS_Features,Shared_Features),]
          TrgSS_W[[view]] = TrgSS_W[[view]][match(Shared_Features,rownames(TrgSS_W[[view]])),]
          
        }
        
        TrgFull_W[[view]] = TrgFull_W[[view]][is.element(TrgFull_Features,Shared_Features),]
        TrgFull_W[[view]] = TrgFull_W[[view]][match(Shared_Features,rownames(TrgFull_W[[view]])),]
        
      }
      
      
      ## standardize and concatenate weight vectors
      for (view in Shared_Views){
      
        if (FeatWseScl_W){
          TrgSS_W[[view]] = t(apply(TrgSS_W[[view]],1, ScaleW))
          TrgFull_W[[view]] = t(apply(TrgFull_W[[view]],1, ScaleW))
        } else {
          TrgSS_W[[view]] = apply(TrgSS_W[[view]],2, ScaleW)
          TrgFull_W[[view]] = apply(TrgFull_W[[view]],2, ScaleW)
        }
        
        TrgSS_W[[view]][is.na(TrgSS_W[[view]])]=0
        TrgSS_W[[view]] = base::scale(TrgSS_W[[view]], center = TRUE, scale = FALSE)
        
        TrgFull_W[[view]][is.na(TrgFull_W[[view]])]=0
        TrgFull_W[[view]] = base::scale(TrgFull_W[[view]], center = TRUE, scale = FALSE)
        
      }
      
      TrgSS_W = do.call(rbind,TrgSS_W)
      TrgFull_W = do.call(rbind,TrgFull_W)
      
  
      ###
      ## best hits, FMeasure and F1 score if appropriate
      ###
      
      # correlation matrices and values
      Z_cors = cor(x = TrgSS_Z, y = TrgFull_Z, method = CorMethod)
      W_cors = cor(x = TrgSS_W, y = TrgFull_W, method = CorMethod)
      
      # correlation pvalues
      Z_cors_pv = Z_cors
      for (r in 1:dim(Z_cors)[1]){
        for (c in 1:dim(Z_cors)[2]){
          Z_cors_pv[r,c] = cor.test(x = TrgSS_Z[,r], y = TrgFull_Z[,c], 
                                    alternative = c("two.sided"), 
                                    method = CorMethod)$p.value
        }
      }
      Z_cors_pv[Z_cors_pv>=corCutOff]=NA
      Z_cors_pv[!is.na(Z_cors_pv)]=1
      
      W_cors_pv = W_cors
      for (r in 1:dim(W_cors)[1]){
        for (c in 1:dim(W_cors)[2]){
          W_cors_pv[r,c] = cor.test(x = TrgSS_W[,r], y = TrgFull_W[,c], 
                                    alternative = c("two.sided"), 
                                    method = CorMethod)$p.value
        }
      }
      W_cors_pv[W_cors_pv>=corCutOff]=NA
      W_cors_pv[!is.na(W_cors_pv)]=1
            
      ## only significant correlations are considered
      Z_cors = Z_cors*Z_cors_pv     
      W_cors = W_cors*W_cors_pv
            
      # Z best hits      
      TrgSS_Z_BHs = data.frame(
        SS = 1:nrow(Z_cors),
        Full = apply(Z_cors,1,GetBestHit)
      )
      TrgSS_Z_BHs = TrgSS_Z_BHs[rowMeans(is.na(Z_cors))<1,]
      
      TrgFull_Z_BHs = data.frame(
        SS = apply(Z_cors,2,GetBestHit),
        Full = 1:ncol(Z_cors)
      )
      TrgFull_Z_BHs = TrgFull_Z_BHs[colMeans(is.na(Z_cors))<1,]
      
      
      # W best hits       
      TrgSS_W_BHs = data.frame(
        SS = 1:nrow(W_cors),
        Full = apply(W_cors,1,GetBestHit)
      )
      TrgSS_W_BHs = TrgSS_W_BHs[rowMeans(is.na(W_cors))<1,]
      
      TrgFull_W_BHs = data.frame(
        SS = apply(W_cors,2,GetBestHit),
        Full = 1:ncol(W_cors)
      )
      TrgFull_W_BHs = TrgFull_W_BHs[colMeans(is.na(W_cors))<1,]

      #ZW best hits
      TrgSS_ZW_BHs = TrgSS_Z_BHs %>%
        dplyr::inner_join(TrgSS_W_BHs, by = c('SS','Full')) %>%
        as.data.frame()
      
      TrgFull_ZW_BHs = TrgFull_Z_BHs %>%
        dplyr::inner_join(TrgFull_W_BHs, by = c('SS','Full')) %>%
        as.data.frame()
           
      ### recovery, relevance and FMeasure based on correlations
      Z_cors[Z_cors<0]=0
      Z_cors[is.na(Z_cors)]=0
            
      Z_Recovery = mean(apply(Z_cors,2,GetBestScore))+0.0001
      Z_Relevance = mean(apply(Z_cors,1,GetBestScore))+0.0001
      Z_FMeasure = 2/((1/Z_Recovery)+(1/Z_Relevance))

      W_cors[W_cors<0]=0
      W_cors[is.na(W_cors)]=0
      
      W_Recovery = mean(apply(W_cors,2,GetBestScore))+0.0001
      W_Relevance = mean(apply(W_cors,1,GetBestScore))+0.0001
      W_FMeasure = 2/((1/W_Recovery)+(1/W_Relevance))
       
      ## F1 scores and true positives if samples from multiple cancers
      if (regexpr('[_]',Prjct)>0){
        
        TrgSS_prjcts = as.factor(Trg_prjcts[[SS]][rownames(TrgSS_Z)])
       
        ## significant factors from SS factorization
        TrgSS_Z_pv = numeric()
        for (fctr in 1:ncol(TrgSS_Z)){
          if (length(levels(TrgSS_prjcts)) > 2){
            DifFactTest = kruskal.test(TrgSS_Z[,fctr] ~ TrgSS_prjcts)
            } else {
              DifFactTest = wilcox.test(TrgSS_Z[,fctr] ~ TrgSS_prjcts)
            }
          TrgSS_Z_pv = c(TrgSS_Z_pv, DifFactTest$p.value)
        }
        TrgSS_SigFct = which(TrgSS_Z_pv<pvalueTH)
        TrgSS_SigFct = unique(TrgSS_W_BHs$Full[is.element(TrgSS_W_BHs$SS,TrgSS_SigFct)])
        
        ## compare to factorization of full target dataset
        
        Inf_positives = length(TrgSS_SigFct)
        True_positives = sum(is.element(TrgSS_SigFct,TrgFull_SigFct))
        if (Inf_positives==0){
          Precision=0
        } else {
          Precision = True_positives/Inf_positives
        }
        if (GT_positives==0){
          Recall=0
        } else {
          Recall = True_positives/GT_positives
        }
        if ((Precision+Recall)==0){
          F1=0
        } else {
          F1 = (2*Precision*Recall)/(Precision+Recall)
        }
        
        prjct_name = c(prjct_name,rep(Prjct,6))
        ss_iteration = c(ss_iteration,rep(ss,6))
        fctrzn_method = c(fctrzn_method,rep(fm,6))
        score_name = c(score_name,
        'GT_positives','Inf_positives','True_positives','Precision','Recall','F1'
        )
        score = c(score,
        GT_positives,Inf_positives,True_positives,Precision,Recall,F1
        )

        factor_matches[['True_positives']][[Prjct]][[SS]][[fm]] = TrgSS_SigFct[is.element(TrgSS_SigFct,TrgFull_SigFct)]
        
        rm(list=c('TrgSS_prjcts','TrgSS_SigFct'))
        
      }
      
      ## record other scores and factor matches
      prjct_name = c(prjct_name, rep(Prjct,14))
      ss_iteration = c(ss_iteration, rep(ss,14))
      fctrzn_method = c(fctrzn_method, rep(fm,14)) 

      score_name = c(score_name,
      'SS_K', 'Full_K', 
      'Z_BHs_SSFull', 'W_BHs_SSFull', 'ZW_BHs_SSFull',
      'Z_BHs_FullSS', 'W_BHs_FullSS', 'ZW_BHs_FullSS',
      'Z_FMeasure', 'W_FMeasure',
      'Z_Recovery', 'Z_Relevance', 'W_Recovery', 'W_Relevance'
      ) 
           
      score = c(score,
      dim(Z_cors)[1], dim(Z_cors)[2],
      length(unique(TrgSS_Z_BHs$Full)), length(unique(TrgSS_W_BHs$Full)), length(unique(TrgSS_ZW_BHs$Full)),
      length(unique(TrgFull_Z_BHs$SS)), length(unique(TrgFull_W_BHs$SS)), length(unique(TrgFull_ZW_BHs$SS)),
      Z_FMeasure, W_FMeasure,
      Z_Recovery, Z_Relevance, W_Recovery, W_Relevance 
      )

      factor_matches[['Z_BHs_SSFull']][[Prjct]][[SS]][[fm]] = TrgSS_Z_BHs
      factor_matches[['W_BHs_SSFull']][[Prjct]][[SS]][[fm]] = TrgSS_W_BHs
      factor_matches[['ZW_BHs_SSFull']][[Prjct]][[SS]][[fm]] = TrgSS_ZW_BHs
      
      ## tidy up before next iteration      
      rm(list=c('TrgSSFctrznDir','TrgSS_Z','TrgSS_W_tmp','TrgSS_Views','Shared_Views', 
                'W_cors','Z_cors',
                'W_cors_pv','Z_cors_pv',
                'TrgSS_Z_BHs','TrgFull_Z_BHs', 'TrgSS_ZW_BHs',
                'TrgSS_W_BHs','TrgFull_W_BHs', 'TrgFull_ZW_BHs',
                'Z_Recovery', 'Z_Relevance', 'Z_FMeasure',
                'W_Recovery', 'W_Relevance', 'W_FMeasure'
                ))
      
    }
    
    rm(list=c('TrgSS_smpls', 'TrgFull_Z', 'TrgFull_smpls_kp'))
    invisible(gc())
    
  }
  
  rm(list=c('TrgFullDir','TrgSSDir','InputModel', 'TrgFull_Fctrzn', 'TrgFull_Views', 
            'TrgSS_meta_data' ,'Trg_smpls',
            'TrgFull_meta_data', 'TrgFull_Z_all', 'TrgFull_SigFct', 'GT_positives'))
  invisible(gc())

}

## combine results into a data frame

FctrznEvaluationResults = data.frame(
  prjct_name = prjct_name,
  ss_iteration = ss_iteration,
  fctrzn_method = fctrzn_method,
  score_name = score_name,
  score = score
)

FctrznEvaluation = list(
  'FctrznEvaluationResults' = FctrznEvaluationResults,
  'factor_matches' = factor_matches,
  'Seed' = Seed,
  'Prjcts' = Prjcts,
  'TopD' = TopD,
  'SS_size' = SS_size,
  'TrgFullK' = TrgFullK,
  'TrgFullTH' = TrgFullTH,
  'TrgSSK' = TrgSSK,
  'FctrznMethods' = FctrznMethods,
  'RandomKRange' = RandomKRange,
  'FeatWseScl_W' = FeatWseScl_W,
  'ScaleW' = ScaleW,
  'CorMethod' = CorMethod,
  'corCutOff' = corCutOff,
  'pvalueTH' = pvalueTH,
  'GetBestHit' = GetBestHit,
  'GetBestScore' = GetBestScore
)

saveRDS(FctrznEvaluation, file.path('Results',
                                    paste0('FctrznEvaluation_',TopD,'D_',TrgFullTH,'.rds'))

)

print("finished")
