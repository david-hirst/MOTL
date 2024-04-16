############
#### PLot factorization results
##################

## libraries
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(nlme)
library(emmeans) #emmeans installed with R not with conda

### which configs

TopD = 5000
TrgFullTH = '01TH'

## import the scores
FctrznEvaluation = readRDS(file.path('Results',
                                     paste0('FctrznEvaluation_',TopD,'D_',TrgFullTH,'.rds'))
                           )

FctrznEvaluationResults = FctrznEvaluation$FctrznEvaluationResults

## rename if necessary
FctrznEvaluationResults$fctrzn_method[FctrznEvaluationResults$fctrzn_method=="TL_VI"]="MOTL"
FctrznEvaluationResults$score_name[FctrznEvaluationResults$score_name=="Z_FMeasure"]="FM_Z"
FctrznEvaluationResults$score_name[FctrznEvaluationResults$score_name=="W_FMeasure"]="FM_W"

## tables of means
score_sum = FctrznEvaluationResults[is.element(FctrznEvaluationResults$score_name,
                                        c('SS_K', 'Full_K', 
                                          'Z_Recovery', 'Z_Relevance', 'W_Recovery', 'W_Relevance',
                                          'FM_Z', 'FM_W',
                                          'GT_positives', 'Inf_positives', 'True_positives',
                                          'Precision', 'Recall', 'F1')) &
                                          is.element(FctrznEvaluationResults$fctrzn_method,
                                        c('Direct','MOTL')) &
                                          is.element(FctrznEvaluationResults$prjct_name,
                                        c('LAML_PAAD','LAML_SKCM','PAAD_SKCM','LAML_PAAD_SKCM')),] %>%
  dplyr::group_by(prjct_name, score_name, fctrzn_method) %>%
  dplyr::summarize(ave = round(mean(score),3)) %>%
  as.data.frame()

score_sum = reshape(score_sum,
                    idvar = c("prjct_name","score_name"),
                    timevar = "fctrzn_method",
                    direction = "wide")

colnames(score_sum) = substr(colnames(score_sum),
                                   regexpr('[.]',colnames(score_sum))+1,
                                   nchar(colnames(score_sum)))

## Save means as a csv
write.table(score_sum,
            file.path('Results','Fctrzn_Mean_Scores.csv'),
            quote = FALSE, sep=',', na='', row.names = FALSE, col.names = TRUE)

###
## Plot the results
###

PlotScoresLS = vector("list")
PlotPrjctsLS = vector("list")

PlotScoresLS[[1]] = c('FM_W', 'FM_Z', 'F1')
PlotScoresLS[[2]] = c('FM_W', 'FM_Z')
PlotScoresLS[[3]] = c('F1')

PlotPrjctsLS[[1]] = c('LAML_PAAD','LAML_SKCM','PAAD_SKCM','LAML_PAAD_SKCM')
PlotPrjctsLS[[2]] = PlotPrjctsLS[[1]]
PlotPrjctsLS[[3]] = PlotPrjctsLS[[1]]

PlotNames = c('Combined', 'FM', 'F1')

PlotHeights = c(10, 10, 10)
PlotWidths = c(10, 10, 7.5)

PlotMethods = c('Direct', 'MOTL')

for (ps in 1:length(PlotScoresLS)){

  PlotScores = PlotScoresLS[[ps]]
  PlotPrjcts = PlotPrjctsLS[[ps]]
  PlotName = PlotNames[ps]
  PlotHeight = PlotHeights[ps]
  PlotWidth = PlotWidths[ps]
  
  DFToPlot = FctrznEvaluationResults[
    is.element(FctrznEvaluationResults$score_name,PlotScores) &
      is.element(FctrznEvaluationResults$prjct_name,PlotPrjcts) &  
      is.element(FctrznEvaluationResults$fctrzn_method,PlotMethods)
    ,]
  
  # groupings and factors for colours and ordering
  DFToPlot$fctrzn_grp = DFToPlot$fctrzn_method
  DFToPlot$fctrzn_grp[!is.element(DFToPlot$fctrzn_grp,c('Direct','Random_TL'))]='other'
  DFToPlot$fctrzn_grp = as.factor(DFToPlot$fctrzn_grp)
  DFToPlot$fctrzn_method = as.factor(DFToPlot$fctrzn_method)
  DFToPlot$fctrzn_method = relevel(DFToPlot$fctrzn_method, ref = 'Direct')
  
  ggplot(DFToPlot, aes(x=fctrzn_method, y=score, fill=fctrzn_grp))+ 
    ylab("") +
    xlab("") +
    geom_violin(lwd=0.25) +
    stat_summary(fun=mean, geom="point", shape=18, color="black", size=1.5) +
    facet_grid(factor(prjct_name, levels=PlotPrjcts)~factor(score_name, levels=PlotScores), 
               scales = 'fixed')+
    coord_flip()+
    theme_bw()+
    theme(
      axis.text.x = element_text(size=6),
      axis.text.y = element_text(size=6),
      strip.text.x = element_text(size=6),
      # strip.text.y = element_text(angle = 0, size=6),
      strip.text.y = element_text(size=5),
      legend.position="none",
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank()
          ) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(n.breaks = 3)
  ggsave(file.path('Results',
                   paste0(PlotName,'_',TopD,'D_',TrgFullTH,'.pdf')),
         width = PlotWidth, height = PlotHeight, units = "cm"
         ) 
  ggsave(file.path('Results',
                   paste0(PlotName,'_',TopD,'D_',TrgFullTH,'.png')),
         width = PlotWidth, height = PlotHeight, units = "cm"
         )
}

print("plotting finished")

###
## Test significance
###

# generalized least squares
# accounting for correlations between scores from the same dataset

print("GLS tests")

PrjctsToTest = list(
  'multiple' = c('LAML_PAAD','LAML_SKCM','PAAD_SKCM','LAML_PAAD_SKCM'),
  'multiple_nonLP' = c('LAML_SKCM','PAAD_SKCM','LAML_PAAD_SKCM')
)
ScoresToTest = list(
  'F1' = names(PrjctsToTest),
  'FM_W' = names(PrjctsToTest),
  'FM_Z' = names(PrjctsToTest),
  'W_Recovery' = names(PrjctsToTest), 
  'W_Relevance' = names(PrjctsToTest)
)

MethodsToTest = c('Direct', 'MOTL')

# loop through scores and project groupings to test
for (STT in names(ScoresToTest)){
  print(STT)
  for (PGTT in ScoresToTest[[STT]]){
    # get the grouping of projects to test
    PsTT = PrjctsToTest[[PGTT]]
    print(PsTT)

    # subset the raw data
    DFToTest = FctrznEvaluationResults[
      is.element(FctrznEvaluationResults$score_name, STT) & 
      is.element(FctrznEvaluationResults$prjct_name, PsTT) & 
      is.element(FctrznEvaluationResults$fctrzn_method, MethodsToTest),]

    # same number of obs per prjct
    # removed as prjct name is a variable in the gls
    # PrjctObs = DFToTest %>%
    #   dplyr::group_by(prjct_name) %>%
    #   dplyr::summarise(max_obs = max(ss_iteration)) %>%
    #   as.data.frame()
    # PrjctObs = min(PrjctObs$max_obs)
    # DFToTest = DFToTest[DFToTest$ss_iteration<=PrjctObs,]

    # make the method a factor
    DFToTest$fctrzn_method = as.factor(DFToTest$fctrzn_method)
    DFToTest$fctrzn_method = relevel(DFToTest$fctrzn_method, ref = "Direct")

    # make datset id and project name factors
    DFToTest$ds_id = as.factor(paste0(DFToTest$prjct_name,"_",DFToTest$ss_iteration))
    DFToTest$prjct_name = as.factor(DFToTest$prjct_name)

    ## gls regression
    FctEval_gls = gls(model = score ~ fctrzn_method + prjct_name,
                      data = DFToTest,
                      correlation = corCompSymm(form = ~1|ds_id),
                      weights = varIdent(form = ~1|fctrzn_method)
    )

    ## test contrasts with emmans
    FctEval_gls_emm = emmeans(FctEval_gls, "fctrzn_method")
    FctEval_gls_ctr = contrast(FctEval_gls_emm, method = "pairwise", adjust = "tukey")

    ## print results
    print(FctEval_gls_ctr)

  }
}

print("testing finished")
