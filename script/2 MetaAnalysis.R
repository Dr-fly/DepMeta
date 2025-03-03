source("R/Functions_PrognosticMetaAnalysis.R")
# 1 load survivaldataList ####
GBMsurvivaldataList <- get(load("data/Data_used_for_MetaAnalysisFunctions/GBMsurvivaldataList_Zscore.Rdata"))
# 2 Get input for meta-analysis ####
GBMMetaInputList <- list()
for (cohortNM in names(GBMsurvivaldataList)){
  survivaldata <- GBMsurvivaldataList[[cohortNM]]
  expr <- t(survivaldata[,4:ncol(survivaldata)])
  expr <- convertGeneName(expr, "Ensembl")
  temp <- t(expr)
  survivaldata_Ens <- cbind(survivaldata[,1:3],temp)
  GBMMetaInputList[[i]] <- GetCoxData(GBMsurvivaldataList[[i]])
}
save(GBMMetaInputList, file = "data/Data_used_for_MetaAnalysisFunctions/GBMMetaInputList.Rdata")
# 3 Perform meta-analysis ####
GBMMetaInputList <- get(load("data/Data_used_for_MetaAnalysisFunctions/GBMMetaInputList.Rdata"))
metaRes <- PrognosticMetaAnalysis(MetaInputList = GBMMetaInputList,
                                  EffectModel = "random",
                                  Pvalue_Thresholds = 0.05, 
                                  HR_Thresholds = 1, 
                                  I2_Thresholds = 0.5, 
                                  Pvalue.Q_Thresholds = 0.1)
# 4 visualization ####
PlotMetaForest(GBMMetaInputList,Symbol = "MCM3",Colors = c("#155289","#94221F"))
