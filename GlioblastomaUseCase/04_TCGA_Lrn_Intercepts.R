
################
#### Import factorization and input data and learn intercepts
###############

## MOFA doesn't learn intercepts but it could be good to have them
## this script can be used to learn intercepts via MLE

library(ggplot2)
library(MOFA2)
library(rhdf5)
library(rjson)

fit_start_time = Sys.time()

Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

### Set directories

ExpDataDir = 'Lrn_5000D'
FctrznDir = file.path(ExpDataDir,'Fctrzn_100K_001TH')

## import meta data and model

expdat_meta = readRDS(file.path(ExpDataDir,'expdat_meta.rds'))

InputModel = file.path(FctrznDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

views = Fctrzn@data_options$views
likelihoods = Fctrzn@model_options$likelihoods
M = Fctrzn@dimensions$M
D = Fctrzn@dimensions$D

# loop through the views and estimate the intercept
# for gaussian data its just the mean
# for other data will try mle with a naive estimator as backup

InterceptsNaive = vector("list")
Intercepts = vector("list")
InterceptsMethod = vector("list")

for (m in 1:M){
  
  InterceptsTmp = numeric()
  InterceptsMethodTmp = character()
  
  view = views[m]
  print(view)
  
  likelihood = likelihoods[which(names(likelihoods)==view)]
  DTmp = D[which(names(D)==view)]
 
  
  YTmp = read.table(file = file.path(ExpDataDir, paste0(view,'.csv')), 
                    sep = ",")
  YTmp = t(as.matrix(YTmp))
  rownames(YTmp) = expdat_meta$smpls
  colnames(YTmp) = expdat_meta[[which(names(expdat_meta) == paste0("ftrs_",view))]] 
  
  ZWTmp = Fctrzn@expectations$Z$group0 %*% 
    t(Fctrzn@expectations$W[[which(names(Fctrzn@expectations$W)==view)]])
  
  invisible(gc())
  
  if (likelihood=="gaussian"){
    
    InterceptsNaive[[m]] = colMeans(YTmp, na.rm = TRUE)
    Intercepts[[m]] = InterceptsNaive[[m]]
    InterceptsMethod[[m]] = rep("Naive",length(Intercepts[[m]]))
    names(InterceptsMethod[[m]]) = names(InterceptsNaive[[m]])
    
  } else if (likelihood=="bernoulli"){
    
    ## naive intercept based on approximation to feature means of ZW
    InterceptsNaive[[m]] = log(
      colMeans(YTmp, na.rm = TRUE)/(1-colMeans(YTmp, na.rm = TRUE))
    )
    
    ## mle estimate of intercept for ZW
    ## if optimiser fails for a feature will return the naive estimate
    
    for (d in 1:DTmp){
      
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
      
      interceptMLEfit = try(stats4::mle(nLL, start=list(InterceptMLE=InterceptsNaive[[m]][d]))@coef[1])
      
      if (class(interceptMLEfit)=="try-error"){
        InterceptsTmp = c(InterceptsTmp,InterceptsNaive[[m]][d])
        InterceptsMethodTmp = c(InterceptsMethodTmp,"Naive")
      } else {
        InterceptsTmp = c(InterceptsTmp,interceptMLEfit)
        InterceptsMethodTmp = c(InterceptsMethodTmp,"MLE")
      }
      
      names(InterceptsTmp)[d] = names(InterceptsNaive[[m]])[d]
      names(InterceptsMethodTmp)[d] = names(InterceptsNaive[[m]])[d]
      
    }
    
    Intercepts[[m]] = InterceptsTmp
    InterceptsMethod[[m]] = InterceptsMethodTmp
    
    
  } else {
    
    InterceptsNaive[[m]] = numeric()
    Intercepts[[m]] = numeric()
    InterceptsMethod[[m]] = character()
    
  }
  
}

names(InterceptsNaive) = views
names(Intercepts) = views
names(InterceptsMethod) = views

## save the intercepts in the relevant factorization folder

fit_end_time = Sys.time()

EstimatedIntercepts = list(
  Seed = Seed,
  InterceptsNaive = InterceptsNaive,
  Intercepts = Intercepts,
  InterceptsMethod = InterceptsMethod,
  fit_start_time = fit_start_time,
  fit_end_time = fit_end_time
)

saveRDS(EstimatedIntercepts,file.path(FctrznDir,"EstimatedIntercepts.rds"))

print("finished")

###
### some testing
###

# EstimatedIntercepts = readRDS(file.path(FctrznDir,"EstimatedIntercepts.rds"))
# Intercepts = EstimatedIntercepts$Intercepts
# InterceptsMethod = EstimatedIntercepts$InterceptsMethod
# InterceptsNaive = EstimatedIntercepts$InterceptsNaive
# 
# 
# hist(Intercepts$mRNA)
# sum(is.na(Intercepts$mRNA))
# 
# hist(Intercepts$miRNA)
# sum(is.na(Intercepts$miRNA))
# 
# hist(Intercepts$DNAme)
# sum(is.na(Intercepts$DNAme))
# 
# table(InterceptsMethod$SNV)
# sum(is.na(Intercepts$SNV))
# 
# rnd_idx = sample(1:length(InterceptsNaive$SNV),1000,replace=FALSE)
# plot(x=InterceptsNaive$SNV[rnd_idx],y=Intercepts$SNV[rnd_idx])
# abline(a=0,b=1)
