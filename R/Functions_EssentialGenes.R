#' Extract DepMap data for concerned cancer type from pan-cancer data ####
#' @param modelNMs ids of cell lines of concerned cancer type in the metadata of DepMap
#'                Manual filtering is required in the OncotreeLineage, OncotreePrimaryDisease and OncotreeSubtype columns of the metadata
#' @param fmetadata metadata for concerned modelNMs
#' @return a list containing Chronos score, metadata, model names, mutation, CNV and expression data for concerned cancer type
#' @examples
#' DepMapList <- get(load("data/Data_used_for_EssentialGenes/DepMapList.Rdata"))
#' metadata <- DepMapList[["metadata"]]
#' fmetadata_GBM <- metadata[metadata$OncotreeLineage=="CNS/Brain",] #164
#' fmetadata_GBM <- fmetadata_GBM[fmetadata$OncotreePrimaryDisease=="Diffuse Glioma",] #132
#' fmetadata_GBM <- fmetadata_GBM[fmetadata$OncotreeSubtype=="Glioblastoma"|fmetadata$OncotreeSubtype=="Glioblastoma Multiforme",] #98
#' modelNMs_GBM <- fmetadata_GBM$ModelID
#' GBMDepMapList<- GetDepMapList(modelNMs_GBM,fmetadata_GBM)
GetDepMapList <- function(modelNMs,fmetadata){
  effect <- DepMapList[["effect"]]
  feffect <- effect[rownames(effect)%in%modelNMs,] 
  feffect <- feffect[, colSums(is.na(feffect)) == 0] #
  
  mut <- DepMapList[["mut"]] 
  fmut <- mut[mut$ModelID%in%modelNMs,] 
  length(unique(fmut$HugoSymbol)) 
  
  CNV <- DepMapList[["CNV"]]
  fCNV <- CNV[rownames(CNV)%in%modelNMs,] 
  
  expr <- DepMapList[["expr"]]
  fexpr <- expr[rownames(expr)%in%modelNMs,] 
  fexpr <- fexpr[,colSums(fexpr == 0) <= nrow(fexpr)/2]
  
  GBMDepMapList <- list(
    effect = feffect,
    metadata = fmetadata,
    modelNM = modelNMs,
    CNV = fCNV,
    mut = fmut,
    expr = fexpr
  )
  GBMDepMapList
}

#' Screening essential genes at multi-omics levels####
#' @param subDepMapList a list containing Chronos score, metadata, model names, mutation, CNV and expression data for concerned cancer type
#'                      output of GetDepMapList function
#' @param omicsNM choose from c("mut","CNV","expr")
#' @return a data frame containing screening results at selected omics level
#'         Genomic levels (mutation by Wilcoxon test and CNV by correlation analysis) 
#'         Transcriptomic level (mRNA by correlation analysis)
#' @examples
#' GBMDepMapList <- get(load("data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata"))
#' GBMDepMapResList <- list(
#' mut = GetEssentialGenes(GBMDepMapList,"mut"),
#' CNV = GetEssentialGenes(GBMDepMapList,"CNV"),
#' expr = GetEssentialGenes(GBMDepMapList,"expr")
#' )
GetEssentialGenes <- function(subDepMapList,omicsNM){
  feffect <- subDepMapList[["effect"]]
  if (omicsNM == "mut"){
    wilcox <- data.frame(matrix(ncol = 7, nrow = 0))
    colnames(wilcox) <- c("symbol", "p.value", 
                          "posi_number","nega_number",
                          "posi_median", "nega_median", "direction")
    fmut <- subDepMapList[["mut"]]
    for (symbol in unique(fmut$HugoSymbol)) {
      mutmodelNM <- fmut[fmut$HugoSymbol==symbol,]$ModelID
      if(symbol%in%colnames(feffect) & any(mutmodelNM%in%rownames(feffect))){
        temp <- data.frame(row.names = rownames(feffect),
                           effect = feffect[,symbol])
        temp$group <- ifelse(rownames(temp)%in%mutmodelNM,1,0)
        
        posi_mid <- median(temp[temp$group == 1, "effect"], na.rm = TRUE)
        nega_mid <- median(temp[temp$group == 0, "effect"], na.rm = TRUE)
        
        posi_number <- nrow(temp[temp$group == 1,])
        nega_number <- nrow(temp[temp$group == 0,])
        
        if (posi_mid < nega_mid) {
          direction <- 1 # 1 means lower Chronos score (more negative scores represent cell lines that are more dependent on the gene for proliferation)
        } else {
          direction <- 0
        }
        result <- wilcox.test(temp$effect ~ temp$group, data = temp)
        wilcox[nrow(wilcox) + 1,] <- c(symbol, result$p.value, 
                                       posi_number,nega_number,
                                       posi_mid, nega_mid, direction)
      }}    
    wilcox$p.value <- as.numeric(wilcox$p.value)
    wilcox$posi_number <- as.numeric(wilcox$posi_number)
    DEResMut <- wilcox
    rownames(DEResMut) <- DEResMut[,1]
    DEResMut
  } else {
    fCNVorExpr <- subDepMapList[[omicsNM]]
    commonrow <- intersect(rownames(fCNVorExpr),rownames(feffect))
    commoncol <- intersect(colnames(fCNVorExpr),colnames(feffect))
    temp1 <- fCNVorExpr[commonrow,commoncol]
    temp2 <- feffect[commonrow,commoncol]
    cor_list <- list()
    for (i in commoncol) {
      # 相关性检验
      cor_res <- cor(x = temp1[,i],y = temp2[,i],method = 'spearman')
      # pvalue
      cor_pval <- cor.test(x = temp1[,i],y = temp2[,i],)$p.value
      # 合并结果
      final_res <- data.frame(gene_name = i,
                              cor_results = cor_res,
                              cor_pvalue = cor_pval)
      # 储存
      cor_list[[i]] <- final_res
    }
    CorRes <- do.call('rbind',cor_list)
    CorRes
  }
}

#' (Visualization) Box plot between mutant cell lines and wild cell lines####
#' @param Symbol gene or genelist from screening at the Genomic(Mutation) level
#' @param subDepMapList a list containing Chronos score, metadata, model names, mutation, CNV and expression data for concerned cancer type
#'                      output of GetDepMapList function
#' @param Colors a vector of colors, for example c("#155289","#94221F")
#' @return box plots showing the results of DepMap screening at the Genomic(Mutation) level
#' @examples
#' GBMDepMapList <- get(load("data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata"))
#' PlotMutBoxFromDepMapList(c("DNMT3A","BCR","PIK3R1","ARHGEF40","GLI2","PROM1"),GBMDepMapList,c("#C25757", "#3A68AE"))
PlotMutBoxFromDepMapList <- function(Symbol,subDepMapList,Colors){
  feffect <- subDepMapList[["effect"]]
  fmut <- subDepMapList[["mut"]]
  input <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(input) <- c("effect","group","symbol")
  for (symbol in geneInmut){
    # symbol <- geneInmut[1]
    mutmodelNM <- fmut[fmut$HugoSymbol%in%symbol,]$ModelID
    temp <- data.frame(row.names = rownames(feffect),
                       effect = feffect[,symbol])
    temp$group <- ifelse(rownames(temp)%in%mutmodelNM,1,0)
    temp$symbol <- rep(symbol,nrow(temp))
    temp[temp == 1] <- "Mut"
    temp[temp == 0] <- "Wild"
    temp$group <- as.factor(temp$group)
    input <- rbind(input,temp)
  }
  library(ggplot2)
  ggplot(input, aes(y = effect, x = group, color = group)) +
    geom_boxplot(aes(x = group, y = effect), notch = FALSE, size = 0.6) +
    geom_jitter(position = position_jitter(0.17), size = 1) +
    theme_bw() +
    theme(strip.background = element_blank(), #把头上的灰框去掉,
          strip.text =element_text(size = 14),
          axis.text = element_text(colour = "black", size = 20),
          axis.title = element_text(size = 20),  # 设置坐标轴标题字体大小
          axis.text.x = element_text(angle = 45, hjust = 1)) +  # 设置x轴文本角度和对齐方式
    labs(x = "Status", y = "CRISPER effect", color = "Status") +
    scale_y_continuous(expand = c(0, 0), limit = c(-1.25,0.8)) +
    ggsignif::geom_signif(aes(x = group, y = effect),
                          comparisons = list(c("Mut", "Wild")),
                          step_increase = 0.1, map_signif_level = TRUE, test = "wilcox.test",
                          textsize = 6) +
    scale_color_manual(values = color) +
    facet_grid(~symbol, drop = TRUE, scale = "free", space = "free_x")
}


#' (Visualization) Correlation scatterplot between copy number/mRNA expression and Chronos score for selected genes####
#' @param Symbol gene or genelist from screening at the Genomic(copy number) level of Transcriptomic (mRNA) level
#' @param subDepMapList a list containing Chronos score, metadata, model names, mutation, CNV and expression data for concerned cancer type
#'                      output of GetDepMapList function
#' @param omicsNM choose from c("mut","CNV","expr")
#' @param Colors a vector of colors, for example c("#155289","#94221F")
#' @return box plots showing the results of DepMap screening at the Genomic(Mutation) level
#' @examples
#' GBMDepMapList <- get(load("data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata"))
#' PlotScatterFromDepMapList(c("LMNB2","E2F3","KIFC1","MCM3"),"CNV",GBMDepMapList,c("#3A68AE","#C25757"))
#' PlotScatterFromDepMapList(c("NFIX","CDK6"),"expr",GBMDepMapList,c("#3A68AE","#C25757"))
PlotScatterFromDepMapList <- function(geneInCNVorExpr,omicsNM,subDepMapList,Colors){
  #'@geneInCNVorExpr "LMNB2" "E2F3"  "KIFC1" "MCM3"
  #'@omicNM "CNV" "expr"
  #'@subDepMapList 
  #'@color color <- c("#94221F",  "#155289")
  feffect <- subDepMapList[["effect"]]
  fdat <- subDepMapList[[omicsNM]]
  commonrow <- intersect(rownames(feffect),rownames(fdat))
  temp1 <- feffect[commonrow,geneInCNVorExpr]
  temp2 <- fdat[commonrow,geneInCNVorExpr]
  index <- match(rownames(temp1),rownames(temp2)) #adjust the order
  temp1 <- temp1[index,]
  
  plotlst <- list()
  for (symbol in geneInCNVorExpr) {
    library(tidyverse)
    library(ggplot2)
    library(ggpubr)
    library(cowplot)
    library(ggExtra)
    input <- as.data.frame(cbind(temp1[,symbol],temp2[,symbol]))
    colnames(input) <- c("CRISPER_effect",symbol)
    p <- ggplot(input,aes_string(x = symbol,y = "CRISPER_effect")) +
      geom_point(size = 6,color = Colors[1],alpha = 0.5) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()
      ) +
      # 添加回归线
      geom_smooth(method = 'lm',se = T,color =  Colors[2],size = 1.5,fill =  Colors[2]) +
      # 添加相关性系数及p值
      stat_cor(method = "spearman",digits = 3,size=6)
    plotlst[[symbol]] <- p
  }
  allplot <- plot_grid(plotlist = plotlst,ncol = length(geneInCNVorExpr) ,align = "hv")
  allplot
}

