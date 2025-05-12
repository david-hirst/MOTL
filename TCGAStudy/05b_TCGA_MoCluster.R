#############
### Run MoCluster on target datasets
#############

library(mogsa)

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
        OutputDir = file.path(InputDir,'MoCluster')
        if (!dir.exists(OutputDir)){
            dir.create(OutputDir)
            }

        ## import meta data    
        expdat_meta = readRDS(file.path(InputDir, "expdat_meta.rds"))

        YData = vector("list")
        for (vt in viewsTrg){
            YData[[vt]] = read.table(file.path(InputDir,paste0(vt,'.csv')), header = FALSE, sep=',')
            YData[[vt]] = as.matrix(YData[[vt]])
            colnames(YData[[vt]]) = as.vector(expdat_meta[["smpls"]])
            rownames(YData[[vt]]) = as.vector(expdat_meta[[paste0('ftrs_',vt)]])
            YData[[vt]] = YData[[vt]][as.vector(rowSums(is.na(YData[[vt]]))==0),] # can't handle missing values
        }

        ## what to put for k
        KtoUse = length(expdat_meta[["smpls"]])

        ## factorize the target dataset with MoCluster
        fit_start_time = Sys.time()
        
        MoCluster_fctrzn = mogsa::mbpca(
        x = YData,
        ncomp = KtoUse, 
        method = "blockLoading",
        scale=TRUE
        )
        
        fit_end_time = Sys.time()

        MoCluster_data = list(
        "KtoUse" = KtoUse,
        "MoCluster_fctrzn" = MoCluster_fctrzn,
        "fit_start_time" = fit_start_time,
        "fit_end_time" = fit_end_time
        )

        saveRDS(MoCluster_data, file.path(OutputDir, "MoCluster_data.rds"))


    }

}

print("finished")
