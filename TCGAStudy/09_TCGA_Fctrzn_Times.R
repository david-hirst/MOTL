
#############
##### Factorization times
#############

## Packages

library(MOFA2)
library(rhdf5)
library(rjson)
library(dplyr)
library(knitr)
library(kableExtra)

#######
### options and functions
#######
TopD = 5000 ## how many features to retain
Lrn_set = paste0('Lrn','_',TopD,'D')
LrnK = 100
Prjcts = c('LAML_PAAD', 'LAML_SKCM','PAAD_SKCM',
           'LAML_PAAD_SKCM')
SS_size = c(5, 5, 5,
            5) ## how many samples in each sub set
names(SS_size)=Prjcts

FctrznMethods = c('Direct', 'TL_VI')

### 
## loop through projects and factorization to record times
###

prjct_name = character()
fctrzn_method = character()
fctrzn_seconds = numeric()
SS_count = integer()

for (Prjct in Prjcts){
  
  TrgSSDir = paste0('Trg_',Prjct,'_SS',SS_size[Prjct],'_',TopD,'D')
  TrgSS_meta_data = readRDS(file.path(TrgSSDir,'expdat_meta_SS.rds'))
  
  for (fm in FctrznMethods){
    
    if (fm == 'Direct'){
      ## add preprocessing times by subset -  stored separately for direct
      for (ss in 1:TrgSS_meta_data$SS_count){
        SS = paste0('SS_',ss)
        TrgSSPrprcssDir = file.path(TrgSSDir,SS)
        TrgSSPrprcss_meta_data = readRDS(file.path(TrgSSPrprcssDir,'expdat_meta.rds'))
        prjct_name = c(prjct_name,Prjct)
        fctrzn_method = c(fctrzn_method,fm)
        fctrzn_seconds = c(fctrzn_seconds,
                        as.numeric(difftime(TrgSSPrprcss_meta_data$script_end_time,
                                            TrgSSPrprcss_meta_data$script_start_time, 
                                            units="secs"))
                        )
        SS_count = c(SS_count,TrgSS_meta_data$SS_count)
        rm(list=c('TrgSSPrprcss_meta_data'))
      }
      ## add factorization time for all subsets
      TrgSSFctrzn_meta_data = rjson::fromJSON(file = file.path(TrgSSDir,'FctrMeta.json'))
      prjct_name = c(prjct_name,Prjct)
      fctrzn_method = c(fctrzn_method,fm)
      fctrzn_seconds = c(fctrzn_seconds,
                      TrgSSFctrzn_meta_data$script_end_time - TrgSSFctrzn_meta_data$script_start_time
      )
      SS_count = c(SS_count,TrgSS_meta_data$SS_count)
      rm(list=c('TrgSSFctrzn_meta_data'))
    } 
    if (fm == 'TL_VI'){
      TrgSSFctrzn_meta_data = readRDS(file.path(TrgSSDir,'TL_VI_meta_data.rds'))
      prjct_name = c(prjct_name,Prjct)
      fctrzn_method = c(fctrzn_method,fm)
      fctrzn_seconds = c(fctrzn_seconds,
                         as.numeric(difftime(TrgSSFctrzn_meta_data$script_end_time,
                                             TrgSSFctrzn_meta_data$loop_start_time, 
                                             units="secs"))
      )
      SS_count = c(SS_count,TrgSS_meta_data$SS_count)
      rm(list=c('TrgSSFctrzn_meta_data'))
    }
  }
  
  rm(list=c('TrgSS_meta_data'))
}


##
## add in learning dataset times
##

## barcode filtering
Lrn_brcds_times = readRDS(file.path('Lrn_barcodes','brcds_times.rds'))
prjct_name = c(prjct_name,Lrn_set)
fctrzn_method = c(fctrzn_method,'Direct')
fctrzn_seconds = c(fctrzn_seconds,
                   as.numeric(difftime(Lrn_brcds_times$script_end_time, 
                                       Lrn_brcds_times$script_start_time,
                                       units="secs"))
)
SS_count = c(SS_count,1)

## preprocessing
Lrn_prprcss_times = readRDS(file.path(Lrn_set,'expdat_meta.rds'))
prjct_name = c(prjct_name,Lrn_set)
fctrzn_method = c(fctrzn_method,'Direct')
fctrzn_seconds = c(fctrzn_seconds,
                   as.numeric(difftime(Lrn_prprcss_times$script_end_time, 
                                       Lrn_prprcss_times$script_start_time,
                                       units="secs"))
)
SS_count = c(SS_count,1)

# factorization
Lrn_fctrzn_times = rjson::fromJSON(file = file.path(Lrn_set,paste0('Fctrzn_',LrnK,'K'),
                                                    'FctrMeta.json'))
prjct_name = c(prjct_name,Lrn_set)
fctrzn_method = c(fctrzn_method,'Direct')
fctrzn_seconds = c(fctrzn_seconds,
                   Lrn_fctrzn_times$script_end_time - Lrn_fctrzn_times$script_start_time)
SS_count = c(SS_count,1)

#### 
#### combine and summarize
#####

Fctrzn_Times_DF = data.frame(
  prjct_name = prjct_name,
  fctrzn_method = fctrzn_method,
  fctrzn_seconds = fctrzn_seconds,
  SS_count = SS_count
)

Fctrzn_Times_DF = Fctrzn_Times_DF %>%
  dplyr::group_by(prjct_name,fctrzn_method) %>%
  dplyr::summarize(fctrzn_seconds = sum(fctrzn_seconds),
                   SS_count = mean(SS_count)) %>%
  dplyr::mutate(fctrzn_seconds_ave = round((fctrzn_seconds/SS_count),0)) %>%
  as.data.frame()

Fctrzn_Times_DF = Fctrzn_Times_DF[,c(1,2,5)]

Fctrzn_Times_DF = reshape(
  Fctrzn_Times_DF,
  idvar = "prjct_name",
  timevar = "fctrzn_method",
  v.names = "fctrzn_seconds_ave",
  direction="wide"
)

colnames(Fctrzn_Times_DF) = substr(colnames(Fctrzn_Times_DF),
                                   regexpr('[.]',colnames(Fctrzn_Times_DF))+1,
                                   nchar(colnames(Fctrzn_Times_DF)))

write.table(Fctrzn_Times_DF,
            file.path('Results','Fctrzn_Times.csv'),
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = TRUE)



