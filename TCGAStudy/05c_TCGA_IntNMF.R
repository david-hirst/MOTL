#############
### Run IntNMF on target datasets
#############


library(rjson)
library(InterSIM)

Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

Prjcts = c('LAML_PAAD', 'LAML_SKCM','PAAD_SKCM','LAML_PAAD_SKCM')
TopD = 5000 ## how many features to retain
SS_size = 5 ## how many samples in each sub set

#### loop through projects

for (Prjct in Prjcts){

    print(Prjct)

    BasePath = paste0('Trg_',Prjct,'_SS',SS_size,'_',TopD,'D')
    

    expdat_meta_base = readRDS(file.path(BasePath,"expdat_meta.rds"))
    if (expdat_meta_base$if_SNV){ 
        viewsTrg = c('mRNA','miRNA','DNAme','SNV')
    } else {
         viewsTrg = c('mRNA','miRNA','DNAme')
    }

    print(viewsTrg)

    SS_count = expdat_meta_base$SS_count

    for (ss in 1:SS_count){
        
        SS = paste0('SS_',ss)
        print(SS)

        InputDir = file.path(BasePath, SS)
        OutputDir = file.path(InputDir,'IntNMF')
        if (!dir.exists(OutputDir)){
            dir.create(OutputDir)
            }

        ## import meta data    
        expdat_meta = readRDS(file.path(InputDir, "expdat_meta.rds"))

        ## import the target dataset

        YData = vector("list")
        for (vt in viewsTrg){
            YData[[vt]] = read.table(file.path(InputDir,paste0(vt,'.csv')), header = FALSE, sep=',')
            YData[[vt]] = as.matrix(YData[[vt]])
            colnames(YData[[vt]]) = as.vector(expdat_meta[["smpls"]])
            rownames(YData[[vt]]) = as.vector(expdat_meta[[paste0('ftrs_',vt)]])
            YData[[vt]] = YData[[vt]][as.vector(rowSums(is.na(YData[[vt]]))==0),] # can't handle missing values
            YData[[vt]] = t(YData[[vt]]) # transform to put features as columns
            YData[[vt]] = sweep(YData[[vt]],2,as.vector(apply(YData[[vt]],2,"min")),"-") #nonnegative values
        }

        ## get a scaling factor for each matrix
        YData_scale = numeric()
        for (vt in seq_along(viewsTrg)){
        YData_scale = c(YData_scale,
                        mean((YData[[vt]] - mean(YData[[vt]]))^2)
                        )
        }
        YData_scale = max(YData_scale)/YData_scale

        ## what to put for k
        KtoUse = length(expdat_meta[["smpls"]])

        ## factorize the target dataset with IntNMF
        fit_start_time = Sys.time()
        
        IntNMF_fctrzn = IntNMF::nmf.mnnals(
        dat = YData,
        k = KtoUse, 
        wt = YData_scale
        )
        
        fit_end_time = Sys.time()

        # retain view names for evaluation
        names(IntNMF_fctrzn$H) = names(YData)

        # export data

        IntNMF_data = list(
        "KtoUse" = KtoUse,
        "IntNMF_fctrzn" = IntNMF_fctrzn,
        "fit_start_time" = fit_start_time,
        "fit_end_time" = fit_end_time
        )

        saveRDS(IntNMF_data, file.path(OutputDir, "IntNMF_data.rds"))


    }

}

print("finished")
