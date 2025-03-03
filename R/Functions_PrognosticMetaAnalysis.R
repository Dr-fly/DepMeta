#' Univariate cox regression (input for prognostic meta-analysis) for all genes ####
#' 
#' @param survivaldata A data frame containing both gene expressions and clinical features.
#'        The column names should indicate the features: "barcode", "overall_survival", "deceased", and gene names.
#'        Gene names must be in the Ensembl format.
#' @return a data frame containing results for all genes in one cohost
#' @examples
#' load("data/Data_used_for_MetaAnalysisFunctions/GBMsurvivaldataList_Zscore.Rdata")
#' GBMMetaInputList <- list()
#' for (cohortNM in names(GBMsurvivaldataList)){
#'   survivaldata <- GBMsurvivaldataList[[cohortNM]]
#'   expr <- t(survivaldata[,4:ncol(survivaldata)])
#'   expr <- convertGeneName(expr, "Ensembl")
#'   temp <- t(expr)
#'   survivaldata_Ens <- cbind(survivaldata[,1:3],temp)
#'   GBMMetaInputList[[i]] <- GetCoxData(GBMsurvivaldataList[[i]])
#' }
GetCoxData <- function(survivaldata){
  library(survminer)
  library(survival)
  library(RegParallel)
  library(DESeq2)
  library(tibble)
  res <- RegParallel(
    data = survivaldata,
    formula = 'Surv(overall_survival, deceased) ~ [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'breslow',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(survivaldata)[4:ncol(survivaldata)],
    blocksize = 1000,
    cores = 100,
    nestedParallel = F,
    p.adjust = "BH")
  #convert ensmble to symbol
  convertGeneName <- function(expr) {
    #Convert ensembl gene id to gene symbol
    expr <- cbind(sub("\\..*", "", rownames(expr)),expr)
    expr <- expr[!duplicated(expr[,1]),]
    rownames(expr) <- expr[,1]
    expr <- expr[,-1]
    geneNM <- get(load("CommonInput/geneNM.Rdata"))
    temp <-
      merge(
        data.frame("ensembl_gene_id" = rownames(expr)),
        data.frame("ensembl_gene_id" = geneNM$ensembl_gene_id,geneNM$external_gene_name),all.x=T)
    temp <- temp[order(temp$ensembl_gene_id),]
    temp <- temp[!duplicated(temp$geneNM.external_gene_name),]
    temp <- temp[!is.na(temp$geneNM.external_gene_name),]
    
    expr <- expr[order(rownames(expr)),]
    temp1 <- temp$ensembl_gene_id
    expr <- expr[rownames(expr)%in%temp1,]
    rownames(expr) <- temp[,2]
    expr
  }
  res <- as.data.frame(res)
  rownames(res) <- res$Variable
  res <- convertGeneName(res, "Symbol")
  res$Variable <- rownames(res)
  res <- res[!is.na(res$Variable),]
  res <- res[!is.na(res$P),]
  res <- res[order(res$P.adjust),]
  res
}

#' Prognostic meta-analysis for all genes ####
#' 
#' @param MetaInputList A list containing data frames, where each data frame represents the meta-analysis input for a specific cohort (output of GetCoxData function).
#' @param EffectModels Choose from c("random", "fixed").
#'        The effect model to be used in the meta-analysis. 
#'        "random" refers to the random-effects model, which assumes that the true effect sizes vary across studies.
#'        "fixed" refers to the fixed-effects model, which assumes that the true effect sizes are the same across studies.
#' @param Pvalue_Thresholds The threshold for the p-value.
#'        The p-value is a measure of the statistical significance in hypothesis testing. 
#'        In the context of meta-analysis, it represents the p-value resulting from the statistical test
#'        that combines the individual study effect sizes to obtain an overall effect size.
#'        Studies with p-values below this threshold will be considered statistically significant.
#'        Default value: 0.05.
#' @param HR_Thresholds The threshold for the hazard ratio.
#'        The hazard ratio (HR) represents the ratio of the hazard rates between two groups. 
#'        Studies with hazard ratios above this threshold will be considered as having a higher risk or impact, 
#'        while studies with hazard ratios below this threshold will be considered as having a lower risk or impact.
#'        Default value: 1.
#' @param I2_Thresholds The threshold for the I2 statistic.
#'        The I2 statistic measures the percentage of total variation across studies due to heterogeneity rather than chance.
#'        Studies with I2 values below this threshold will be considered as having low heterogeneity.
#'        Default value: 0.5.
#' @param Pvalue.Q_Thresholds The threshold for the p-value of the Q statistic.
#'        The Q statistic tests the null hypothesis of no heterogeneity. 
#'        Studies with p-values of the Q statistic above this threshold will be considered as having low heterogeneity.
#'        Default value: 0.1.
#'
#' @return a list containing total result of meta-analysis, filtered result of meta-analysis and prognostic genes from prognositc meta-analysis
#' @examples
#' MetaInputList <- get(load("data/Data_used_for_MetaAnalysisFunctions/GBMMetaInputList.Rdata"))
#' PrognosticMetaAnalysis(MetaInputList = MetaInputList,
#'                        EffectModel = "random",
#'                        Pvalue_Thresholds = 0.05, 
#'                        HR_Thresholds = 1, 
#'                        I2_Thresholds = 0.5, 
#'                        Pvalue.Q_Thresholds = 0.1)
PrognosticMetaAnalysis <- function(MetaInputList,
                                   EffectModel,
                                   Pvalue_Thresholds = 0.05, 
                                   HR_Thresholds = 1, 
                                   I2_Thresholds = 0.5, 
                                   Pvalue.Q_Thresholds = 0.1){
  library(meta)
  symbols <- rownames(MetaInputList[[1]])
  res <- data.frame(matrix(ncol = 9, nrow = 0))
  colnames(res) <- c("symbol", "HR.random","pval.random","HR.fixed","pval.fixed","I2", "pval.Q", "lower.random", "upper.random")
  x <- 1  
  for (symbol in symbols) {
    tryCatch({
      y <- NULL
      for (i in 1:length(MetaInputList)) {
        df <- MetaInputList[[i]]
        if (symbol %in% rownames(df)) {
          y <- rbind(y, df[symbol, ])
        }
      }
      y$dataset <- names(MetaInputList)
      z <- y[,c("dataset","Variable","HR","HRlower","HRupper")]
      mg0 <- metagen(log(HR),lower=log(HRlower), upper=log(HRupper),data = z, 
                     sm="HR", studlab = paste(z$dataset, z$Variable, sep = "-"))
      mg0$TE.random
      mg0$TE.common
      mg0$TE.fixed
      mg0$I2
      
      res[x,1] <- symbol
      res[x,2] <- exp(mg0$TE.random)
      res[x,3] <- mg0$pval.random
      res[x,4] <- exp(mg0$TE.fixed)
      res[x,5] <- mg0$pval.fixed
      res[x,6] <- mg0$I2  
      res[x,7] <- mg0$pval.Q  
      res[x,8] <- exp(mg0$lower.random)
      res[x,9] <- exp(mg0$upper.random)
      x <- x+1
    }, error = function(e) {
      # Handle error cases in the error function
      message(paste("Error occurred for variable", symbol, ": ", e$message))
    })
  }
  pvaluecol <- paste("pval.",EffectModel,sep="")
  HRcol <- paste("HR.",EffectModel,sep="")
  lowercol <- paste("lower.",EffectModel,sep="")
  uppercol <- paste("upper.",EffectModel,sep="")
  filtered_res <- res[res[,pvaluecol] < Pvalue_Thresholds,]
  filtered_res <- filtered_res[filtered_res[,HRcol] > HR_Thresholds & filtered_res[,lowercol] > 1 |
                                 filtered_res[,HRcol] < HR_Thresholds & filtered_res[,uppercol] > 1,]
  filtered_res <- filtered_res[filtered_res[,"I2"] < I2_Thresholds & filtered_res[,"pval.Q"] > Pvalue.Q_Thresholds,]
  PrognosticMetagene <- filtered_res$symbol
  return(list(res = res, filtered_res = filtered_res, PrognosticMetagene = PrognosticMetagene))
}
#' (Visualization) Forest Plot ####
#' 
#' @param MetaInputList A list containing data frames, where each data frame represents the meta-analysis input for a specific cohort (output of GetCoxData function)
#' @param Symbol gene symbols you want to display (must be genes in the MetaInputList)
#' @param Colors a vector of colors, for example c("#155289","#94221F")
#' @return a forest plot showing the result of prognostic meta-analysis for one gene
#' @examples
#' GBMMetaInputList <- get(load("data/GBMMetaInputList.Rdata"))
#' PlotMetaForest(GBMMetaInputList,Symbol = "MCM3",Colors = c("#155289","#94221F"))
PlotMetaForest <- function(MetaInputList, Symbol, Colors){
  library(meta)
  y <- NULL
  for (i in 1:length(MetaInputList)) {
    df <- MetaInputList[[i]]
    if (Symbol %in% rownames(df)) {
      y <- rbind(y, df[Symbol, ])
    }
  }
  y$dataset <- names(MetaInputList)
  z <- y[,c("dataset","Variable","HR","HRlower","HRupper")]
  
  mg0 <- metagen(log(HR),lower=log(HRlower), upper=log(HRupper),data = z, 
                 sm="HR", studlab = paste(z$dataset, z$Variable, sep = "-"))
  forest(
    mg0,
    layout = "JAMA",
    col.square = Colors[1],
    fontfamily = "serif",
    col.square.lines = Colors[1],
    col.diamond.random = Colors[2],
    col.diamond.lines = Colors[2],
    studlab = T,
    common = F,
    weight.study = "same",
    spacing = 1,
    xlab = "HR (high versus low gene expression)"
  )
}

