##############
#### preprocess target data and then factorize with MOTL
##############

#############################################
### Part 0: Libraries and packages
#############################################

## had to install venn diagram package here as conda failed

library(SummarizedExperiment)
library(DESeq2)
library(sesame)
library(dplyr)
library(MOFA2)
library(rhdf5)
library(fgsea)
library(msigdbr)
library(VennDiagram)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)

source('TCGA_preprocessedData_functions.R')
source('TL_VI_functions.R')

## SET SEED
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)

#############################################
### Part 1: Import and prepare data from the article
#############################################

##
### parameters and directories
##

TopD = 5000
LrnK = 100
LrnVarTH = '001TH'

LrnDir = paste0('Lrn_',TopD,'D')
LrnFctrnDir = file.path(LrnDir,paste0('Fctrzn_',LrnK,'K_',LrnVarTH))

InDir = file.path('RawData')
InRawCountsDir = file.path(InDir,'CountData')
HGCCgroups = c("Normal","CL","MS","PN")
OutDir = file.path(paste0('Trg_',paste(HGCCgroups,collapse='_'),'_',TopD,"D"),"MOTL")
if(!dir.exists(OutDir)){
  dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
}

# MOFA direct directory
MOFAdirectDir = file.path(paste0('Trg_',paste(HGCCgroups,collapse='_'),'_',TopD,"D"),
                          'MOFA', 'Fctrzn_13K_01TH')

####
## functions
####

mRNA_raw_import = function(x){
  # import raw counts line by line for each sample
  mRNA_trg_tmp = strsplit(
    readLines(file.path(InRawCountsDir, paste0(x,'.genes.results')), n = -1), 
    split = "\t"
  )
  
  correct_length = 7
  mRNA_trg_tmp_chk = lapply(mRNA_trg_tmp, function(x) length(x) == correct_length)
  mRNA_trg_tmp = mRNA_trg_tmp[unlist(mRNA_trg_tmp_chk)]
  mRNA_trg_tmp = do.call(rbind, mRNA_trg_tmp)
  colnames(mRNA_trg_tmp) = mRNA_trg_tmp[1,]
  mRNA_trg_tmp = mRNA_trg_tmp[-1,]
  
  mRNA_trg_tmp = data.frame(
    gene_id = mRNA_trg_tmp[,'gene_id'],
    sample = x,
    expected_count = round(as.numeric(mRNA_trg_tmp[,'expected_count']))
  )
  
  return(mRNA_trg_tmp)
}

mRNA_stripVersion = function(expdat){
  # strip version and gene symbol from ensembl ids
  vStart = regexpr('[.]',rownames(expdat))
  suffixStart = regexpr('_PAR_Y',rownames(expdat))
  suffixEnd = suffixStart -1 + attr(regexpr('_PAR_Y',rownames(expdat)),"match.length")
  rownames(expdat) = paste0(
    substr(rownames(expdat),1,vStart-1),
    substr(rownames(expdat), suffixStart, suffixEnd)
  )
  # return tidied up matrix
  return(expdat)
}


##
### sample metadata
##

smpls.df = read.table(
  file.path(InDir,"phenodata_RNAseq.tsv"),
  header = TRUE
)
smpls.df$Sample = sub("[.]","-",smpls.df$Sample)
# select samples for factorization
smpls = smpls.df$Sample[is.element(smpls.df$HGCC,HGCCgroups)]
smpls = sample(smpls, length(smpls), replace = FALSE)

# filter and sort sample table and add cancer vs control variable
smpls.df = smpls.df[match(smpls,smpls.df$Sample),]
smpls.df$CaseControl = smpls.df$HGCC
smpls.df$CaseControl[smpls.df$CaseControl!='Normal']='Case'

# specify grouping variable for testing downstream
sample_group = as.factor(smpls.df$CaseControl)
# sample_group = as.factor(smpls.df$HGCC)

#####
# mRNA
####

## import data and make raw counts matrix

smpls_raw = sub("-",".",smpls)
expdat_mRNA = lapply(X = smpls_raw, FUN = mRNA_raw_import)
expdat_mRNA = do.call(rbind, expdat_mRNA)
expdat_mRNA = reshape(expdat_mRNA,
                      idvar = "gene_id",
                      timevar = "sample",
                      direction = "wide")
rownames(expdat_mRNA) = expdat_mRNA$gene_id
expdat_mRNA = as.matrix(expdat_mRNA[,-1])
colnames(expdat_mRNA) = substr(colnames(expdat_mRNA),16,nchar(colnames(expdat_mRNA)))

### tidy up feature and sample names
colnames(expdat_mRNA) = sub("[.]","-",colnames(expdat_mRNA))
expdat_mRNA = mRNA_stripVersion(expdat_mRNA)

# filter and order samples
expdat_mRNA = expdat_mRNA[,smpls]

#####
# DNAme
####

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

# filter and order samples
expdat_DNAme = expdat_DNAme[,smpls]

#############################################
### Part 2: Use the TCGA TL preprocessing functions
#############################################

## TL PARAMETERS
# LrnSimple = TRUE ## if TRUE then E[W^2] and E[LnTau] are calculated form E[W] and E[Tau] otherwise imported
CenterTrg = FALSE # Center Trg with own means or use estimated intercepts?
if (CenterTrg){
  TLDirName = 'TL_VI_C'
} else {
  TLDirName = 'TL_VI'
}

# ------------------------------------------- LEARNING SET ----

script_start_time = Sys.time()

# IMPORT FACTORIZATION RESULTS

## METADATA
expdat_meta_Lrn = readRDS(file.path(LrnDir,"expdat_meta.rds"))

## FATORIZATION MODEL OBJECT
InputModel = file.path(LrnFctrnDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

viewsLrn = Fctrzn@data_options$views
likelihoodsLrn = Fctrzn@model_options$likelihoods
MLrn = Fctrzn@dimensions$M

## LOAD IN THE EXPECTATIONS FOR TAU AND REORDER FOR CONSISTENCY
Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)

## E[log(tau)] CALCULATION (ONLY RELEVANT FOR GAUSSIAN DATA)
Fctrzn@expectations[["TauLn"]] = sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnDir)

## load or calculate E[W^2] values
Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir)

## import the W intercepts if not centering, create vectors containing 0s if centering
Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg, Fctrzn, LrnFctrnDir)


# -------------------------------------------------------------

# --------------------------------------------- TARGET SET ----

ss_start_time = Sys.time()

## add versions to target dataset mRNA names that are consistant with Lrn set

# mRNA_addVersion = function(expdat, Lrndat){
#   # get versions from learning dataset
#   tmp = as.data.frame(do.call(rbind,strsplit(rownames(Lrndat),"[.]")))
#   # match to stripped ids from target set
#   tmp = data.frame(V1 = rownames(expdat)) %>%
#     dplyr::left_join(tmp, by = c('V1')) %>%
#     as.data.frame()
#   # rename target dataset features
#   rownames(expdat) = paste0(tmp$V1,'.',tmp$V2)
#   # return tidied up matrix
#   return(expdat)
# }

expdat_mRNA = mRNA_addVersion(expdat = expdat_mRNA, Lrndat = Fctrzn@expectations$W$mRNA)

## create a list of all target dataset omics data
YTrg_list = list(
  mRNA = expdat_mRNA,
  DNAme = expdat_DNAme
)

## harmonize the vectors of views and likelihoods
viewsTrg = names(YTrg_list)
views = viewsLrn[is.element(viewsLrn,viewsTrg)]
likelihoods = likelihoodsLrn[views]

## Prepare YTrg data
YTrg_list = TargetDataPreparation(views = views, YTrg_list = YTrg_list, 
                                  Fctrzn = Fctrzn, smpls = smpls, expdat_meta_Lrn = expdat_meta_Lrn,
                                  normalization = 'Lrn', transformation = TRUE)

## INIT PARAMETERS
TL_param = initTransferLearningParamaters(YTrg = YTrg_list, views = views, 
                                          expdat_meta_Lrn = expdat_meta_Lrn, 
                                          Fctrzn = Fctrzn, likelihoods = likelihoods
                                          )

#############################################
### Part 3: Use the Transfer learning function
#############################################

# defining output folder here - this should be an input to the funciton
TL_OutDir = file.path(OutDir, TLDirName)
if(!dir.exists(TL_OutDir)){
  dir.create(TL_OutDir, showWarnings = FALSE, recursive = TRUE)
}

## PARAMETERS
minFactors = 13 ## floor when dropping factors - number of samples in evaluations
StartDropFactor = 1 # after which iteration to start dropping factors
FreqDropFactor = 1 ## how often to drop factors
StartELBO = 1 #which iteration to start checking ELBO on, excl initiation iteration
FreqELBO = 5 #how often to assess the ELBO
DropFactorTH = 0.01 # factor with lowest max variance, that is less than this, is dropped
MaxIterations = 10000
MinIterations = 2 # as per MOFA with defaults in python, >= 2 (hardcoded floor), exl initial setup
ConvergenceIts = 2 # as per MOFA python defaults
ConvergenceTH = 0.0005 # default for MOFA - fast option
# PoisRateCstnt = 0.0001 ## amount to add to the poison rate function to avoid errors

## run the MOTL function

TL_data = transferLearning_function(TL_param = TL_param, MaxIterations = MaxIterations, MinIterations =  MinIterations, 
                                    minFactors = minFactors, StartDropFactor = StartDropFactor, FreqDropFactor = FreqDropFactor, 
                                    StartELBO = StartELBO, FreqELBO = FreqELBO, DropFactorTH = DropFactorTH, 
                                    ConvergenceIts = ConvergenceIts, ConvergenceTH = ConvergenceTH, 
                                    CenterTrg = CenterTrg, ss_start_time = ss_start_time, outputDir = TL_OutDir)

script_end_time = Sys.time()
#### save the metadata

TL_meta_data = list(
  'Seed' = Seed,
  'CenterTrg' = CenterTrg,
  'StartDropFactor' = StartDropFactor,
  'FreqDropFactor' = FreqDropFactor,
  'DropFactorTH' = DropFactorTH,
  'minFactors' = minFactors,
  'MinIterations' = MinIterations,
  'MaxIterations' = MaxIterations,
  'StartELBO' = StartELBO,
  'FreqELBO' = FreqELBO,
  'ConvergenceTH' = ConvergenceTH,
  'ConvergenceIts' = ConvergenceIts,
  'script_start_time' = script_start_time,
  'script_end_time' = script_end_time
)
saveRDS(TL_meta_data, file.path(OutDir,paste0(TLDirName,'_meta_data.rds')))

#############################################
### Part 3: Find associated pathways and processes
#############################################

## import saved rds file
TL_data = readRDS(file.path(TL_OutDir,"TL_data.rds"))

YGauss = TL_data$YGauss
ZMu_0 = TL_data$ZMu_0
ZMu = TL_data$ZMu
Fctrzn_Lrn_W0 = TL_data$Fctrzn_Lrn_W0
Fctrzn_Lrn_W = TL_data$Fctrzn_Lrn_W
views = names(YGauss)
factorNames = colnames(ZMu)
VarExpl = TL_data$VarExpl

##
## test for differentially active factors
##


# test if can use wilcox
minSamples = min(unlist(table(sample_group)))
# testing function
DifFactTester = function(x, minSamples) {
  if (minSamples >= 4) {
    wilcox.test(x ~ sample_group)$p.value
  } else {
    t.test(x ~ sample_group)$p.value
  }
}

Z_pv = apply(X = ZMu, MARGIN = 2, FUN = DifFactTester, minSamples = minSamples)
Z_pv = p.adjust(Z_pv, method="BH")
Z_pv_sig = which(Z_pv < 0.05)
names(Z_pv_sig) = names(Z_pv)[Z_pv_sig]

##
## retrieve variance explained for mRNA
##

Z_mRNA_rel = which(VarExpl[,'mRNA']>0.01)
names(Z_mRNA_rel) = rownames(VarExpl)[Z_mRNA_rel]

####
## GSEA mRNA weights for factors that are differentially active and relevant to mRNA
####

#  scale and tidy row names of mRNA W matrix
GroupsTested = paste(levels(sample_group),collapse='_vs_')
ScaleW = function(x) {
  x/norm(matrix(x),type="F")
}
W_mRNA = t(apply(Fctrzn_Lrn_W[['mRNA']],1, ScaleW))
W_mRNA = mRNA_stripVersion(W_mRNA)

# factors and genestes to test
kToTest = names(Z_pv_sig)[is.element(names(Z_pv_sig),names(Z_mRNA_rel))]
GenesetCatToTest = c("C2", "C2", "C5", "C5")
GenesetSubCatToTest = c("CP:KEGG", "CP:REACTOME", "GO:BP", "GO:CC")

ToTest = data.frame(
  k = rep(kToTest, each=length(GenesetCatToTest)),
  category = rep(GenesetCatToTest,times=length(kToTest)),
  subcategory = rep(GenesetSubCatToTest,times=length(kToTest))
)


SigGeneSets = apply(X = ToTest, MARGIN = 1, FUN = function(ToTestVec, WMat = W_mRNA, Groups = GroupsTested){
  
  k = as.character(ToTestVec[1])
  category = as.character(ToTestVec[2])
  subcategory = as.character(ToTestVec[3])
  rankings = WMat[,k]
  names(rankings) = rownames(WMat)
  
  gene_sets = msigdbr(
    species = "Homo sapiens", 
    category = category, 
    subcategory = subcategory
  )
  
  gene_sets = gene_sets %>% 
    dplyr::distinct(gs_name, ensembl_gene) %>% 
    as.data.frame()
  
  gene_sets = split(
    x = gene_sets$ensembl_gene, 
    f = gene_sets$gs_name
  )
  # set seed
  Seed = 1234567
  mode(Seed) = 'integer'
  set.seed(Seed)
  
  # run fgsea
  fgseaRes = fgsea(
    pathways = gene_sets,
    stats = rankings,
    minSize = 15,
    maxSize = 500
  )
  
  fgseaRes = fgseaRes %>% 
    dplyr::select(pathway, padj, ES, NES, size) %>% 
    as.data.frame()
  
  fgseaRes = fgseaRes[fgseaRes$padj<0.01 & !is.na(fgseaRes$padj),]
  fgseaRes$k = rep(k,nrow(fgseaRes))
  fgseaRes$subcategory = rep(subcategory,nrow(fgseaRes))
  fgseaRes$Groups = rep(Groups,nrow(fgseaRes))
  
  return(fgseaRes)
  
})

SigGeneSets = do.call(rbind, SigGeneSets)

## export as csv
write.table(SigGeneSets,
            file.path(TL_OutDir, paste0('SigGeneSets','_',GroupsTested,'.csv')),
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = TRUE)

#############################################
### Part 4: compare to direct MOFA
#############################################

# MOTL sig pathways
SigGeneSets_MOTL = read.table(file.path(TL_OutDir, paste0('SigGeneSets','_',GroupsTested,'.csv')),
                              header = TRUE, sep=',')
SigGeneSets_MOTL$Method = rep('MOTL',nrow(SigGeneSets_MOTL))

# MOFA sig pathways
SigGeneSets_MOFAdirect = read.table(file.path(MOFAdirectDir, paste0('SigGeneSets','_',GroupsTested,'.csv')),
                                    header = TRUE, sep=',')

SigGeneSets_MOFAdirect$Method = rep('Direct',nrow(SigGeneSets_MOFAdirect))

## venn diagram by category

GeneSetsPlotsDir = file.path(TL_OutDir,'GeneSetsPlots')
if(!dir.exists(GeneSetsPlotsDir)){
  dir.create(GeneSetsPlotsDir, showWarnings = FALSE, recursive = TRUE)
}

for (PlotSubCat in unique(SigGeneSets_MOTL$subcategory)){
  VennX = list(
    MOTL = unique(SigGeneSets_MOTL$pathway[SigGeneSets_MOTL$subcategory==PlotSubCat]),
    Direct = unique(SigGeneSets_MOFAdirect$pathway[SigGeneSets_MOFAdirect$subcategory==PlotSubCat])
  )
  
  VennDiagram::venn.diagram(
    x = VennX,
    filename = file.path(GeneSetsPlotsDir, paste0('CompareSigGeneSets_',PlotSubCat,'_',GroupsTested,'.png')),
    disable.logging = TRUE,
    category.names = names(VennX),
    imagetype="png"
  )
}

# venn diagram all

VennX = list(
  MOTL = unique(SigGeneSets_MOTL$pathway),
  Direct = unique(SigGeneSets_MOFAdirect$pathway)
)

VennDiagram::venn.diagram(
  x = VennX,
  filename = file.path(GeneSetsPlotsDir, paste0('CompareSigGeneSets_',GroupsTested,'.png')),
  disable.logging = TRUE,
  category.names = names(VennX),
  imagetype="png"
)

#############################################
### Part 5: check if overlap is significant
#############################################

categories = c("C2", "C2", "C5", "C5")
subcategories = c("CP:KEGG", "CP:REACTOME", "GO:BP", "GO:CC")
minGSsize = 15
maxGSsize = 500 

k_all = 0
m_all = 0
q_all = 0
n_all = 0

for (sc in 1:length(subcategories)){
  
  category = categories[sc]
  subcategory = subcategories[sc]
  
  gene_sets = msigdbr(
    species = "Homo sapiens", 
    category = category, 
    subcategory = subcategory
  ) %>% 
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(genes = n_distinct(ensembl_gene))%>% 
    as.data.frame()
  
  gene_sets = gene_sets$gs_name[(gene_sets$genes >= minGSsize) & (gene_sets$genes <= maxGSsize)]
  
  ## test for significant overlap
  k = unique(SigGeneSets_MOTL$pathway[SigGeneSets_MOTL$subcategory==subcategory])
  m = unique(SigGeneSets_MOFAdirect$pathway[SigGeneSets_MOFAdirect$subcategory==subcategory])
  q = k[is.element(k,m)]
  n = gene_sets[!is.element(gene_sets,m)]
  
  hypertest = phyper(q=length(q)-1, m=length(m), n=length(n), k=length(k), lower.tail = FALSE)
  
  print(paste0("For ",subcategory,' genesets, the overlap was ',length(q),
               ', with P(X >= x): ',round(hypertest,6)))
  
  
  k_all = k_all + length(k)
  m_all = m_all + length(m)
  q_all = q_all + length(q)
  n_all = n_all + length(n)
  
}

hypertest_all = phyper(q=q_all-1, m=m_all, n=n_all, k=k_all, lower.tail = FALSE)
print(paste0('For all genesets, the overlap was ',q_all,
             ', with P(X >= x): ',round(hypertest_all,6)))

########
## Create a combined CSV of genesets
#######

gseaAllLookup = data.frame(
  BaseDir = c(
    'Trg_Normal_CL_5000D',
    'Trg_Normal_MS_5000D', 
    'Trg_Normal_PN_5000D', 
    'Trg_Normal_CL_MS_PN_5000D'
    ),
  MOFADir = c(
    'Fctrzn_7K_01TH',
    'Fctrzn_7K_01TH',
    'Fctrzn_7K_01TH',
    'Fctrzn_13K_01TH'
  ),
  Comparison = c(
    'Normal_vs_CL',
    'Normal_vs_MS',
    'Normal_vs_PN',
    'Normal_vs_ALL'
  )
)

for (i in 1:dim(gseaAllLookup)[1]){
  # MOTL sig pathways
  SigGeneSets_MOTL = read.table(file.path(gseaAllLookup$BaseDir[i], 'MOTL', 'TL_VI', 
                                          'SigGeneSets_Case_vs_Normal.csv'),
                                header = TRUE, sep=',')
  SigGeneSets_MOTL$Method = rep('MOTL',nrow(SigGeneSets_MOTL))
  SigGeneSets_MOTL$Groups = rep(gseaAllLookup$Comparison[i],nrow(SigGeneSets_MOTL))
  # MOFA sig pathways
  SigGeneSets_MOFAdirect = read.table(file.path(gseaAllLookup$BaseDir[i], 'MOFA', gseaAllLookup$MOFADir[i], 
                                                'SigGeneSets_Case_vs_Normal.csv'),
                                      header = TRUE, sep=',')
  
  SigGeneSets_MOFAdirect$Method = rep('MOFAdirect',nrow(SigGeneSets_MOFAdirect))
  SigGeneSets_MOFAdirect$Groups = rep(gseaAllLookup$Comparison[i],nrow(SigGeneSets_MOFAdirect))
  # combine
  if(i==1){
    SigGeneSets_Combined = rbind(SigGeneSets_MOTL,SigGeneSets_MOFAdirect)
  } else{
    SigGeneSets_Combined = rbind(SigGeneSets_Combined,SigGeneSets_MOTL)
    SigGeneSets_Combined = rbind(SigGeneSets_Combined,SigGeneSets_MOFAdirect)
  }
}

## export as csv
write.table(SigGeneSets_Combined,
            file.path('Results', 'SigGeneSets_Combined.csv'),
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = TRUE)

########
## Heatmap of Z values
#######

#order the samples
smpls.df = smpls.df[order(smpls.df$CaseControl,smpls.df$HGCC,smpls.df$Sample),]
smpls.df$CaseControl[smpls.df$CaseControl=='Case']='pd-GBSC'

# motl factorization
TL_fctrzn = readRDS(file.path(TL_OutDir,"TL_data.rds"))
MOFAdirect_fctrzn = load_model(file = file.path(MOFAdirectDir,"Model.hdf5"))

# mofa factorization
ZMu_List = list(
  Direct = t(MOFAdirect_fctrzn@expectations$Z$group0[smpls.df$Sample,]),
  MOTL = t(TL_fctrzn$ZMu[smpls.df$Sample,])
)


heatMapAnno = data.frame(
  SubType = as.factor(smpls.df$HGCC),
  Type = as.factor(smpls.df$CaseControl)
)
heatMapAnno$SubType = relevel(heatMapAnno$SubType,ref='Normal')
heatMapAnno$Type = relevel(heatMapAnno$Type,ref='Normal')
rownames(heatMapAnno) = smpls.df$Sample

annot_colors <- list(SubType = brewer.pal(5, "Set2")[c(1,3,4,5)],
                   Type = brewer.pal(5, "Set2")[1:2])
names(annot_colors$SubType) <- levels(heatMapAnno$SubType)
names(annot_colors$Type) <- levels(heatMapAnno$Type)

heatmaps_list = lapply(X = ZMu_List, FUN = pheatmap, 
                  scale = "row",
                  border_color = NA,
                  cellwidth = 9,
                  cellheight = 3,
                  fontsize = 6,
                  angle_col = 45,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  # clustering_method = "single",
                  treeheight_row = 0,
                  treeheight_col = 10,
                  cutree_cols = 1,
                  legend = FALSE,
                  show_rownames = FALSE,
                  show_colnames = TRUE,
                  annotation_col = heatMapAnno,
                  annotation_colors = annot_colors,
                  annotation_names_col = TRUE,
                  annotation_legend = TRUE)
heatmaps_list = lapply(heatmaps_list, FUN = function(x){x[[4]]})

heatmaps_g = grid.arrange(arrangeGrob(grobs = heatmaps_list, ncol=1))
ggsave(filename = file.path('Results', 'heatmaps_0clusters.png'), 
       plot = heatmaps_g,
       width = 7.5, height = 10, units = "cm")


write.table(smpls.df,
            file.path('Results', 'SampleTable.csv'),
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = TRUE)
