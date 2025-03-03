#' Wash clinical features ####
#' 
#' @param clinical a dataframe containing at least three renamed variables including 
#' barcode(id of patients), overall_survival(days) and deceased(TRUE/1 means "Death",FALSE/0 means "Alive")
#' @return a data frame including clinical features used for survival analysis without NA
#' @examples
#' clinical_cleaned <- CliClean(clinical)
CliClean <- function(clinical){
  clinical <- clinical[!duplicated(clinical$barcode) & !is.na(clinical$barcode),]
  clinical <- clinical[!is.na(clinical$overall_survival),]
  rownames(clinical) <- clinical$barcode
  clinical$overall_survival <- as.numeric(clinical$overall_survival)
  clinical
}


#' Merge expression matrix and clinical features ####
#' 
#' @param expr Gene Expression Matrix, rows are genes, columns are samples
#' @param clinical a dataframe containing at least three variables including 
#' barcode(id of patients), overall_survival(days) and deceased(TRUE/1 means "Death",FALSE/0 means "Alive")
#' @return a data frame (survivaldata) containing both gene expressions and clinical features, colnames are features 
#' (barcode, overall_survival, deceased, and genes), rownames are patients(barcodes)
#' @examples
#' survivaldata <- MergeCliExp(expr, clinical)
MergeCliExp <- function (expr, clinical) {
  clinical <- clinical[,c("barcode", "overall_survival", "deceased")]
  expr <- as.data.frame(t(expr))
  rownames(expr) <- gsub("\\.","-",rownames(expr)) 
  common <- intersect(rownames(expr),rownames(clinical))
  expr <- expr[common,]
  clinical <- clinical[common,]
  temp <- merge(clinical, expr, by = 0, all.x = T) #colnames(expr)=barcode, rownames(clinical)=barcode
  rownames(temp) <- temp[, 1]
  temp <- temp[, -1]
  temp
}

#' Convert name of genes from Ensembl(ENSG00000141510) to symbol(TP53), or conversely####
#' 
#' @param expr Gene Expression Matrix, rownames are genes, colnames are samples
#' @param toWhichType choose from "Symbol" or "Ensembl"
#' @return expression matrix with expected gene names
#' @examples
#' expr_Symbol <- ConvertGeneName(expr_Ensembl,"Symbol")
ConvertGeneName <- function(expr,toWhichType) {
  expr <- cbind(sub("\\..*", "", rownames(expr)),expr) # "ENSG00000000003.15" to ENSG00000000003
  expr <- expr[!duplicated(expr[,1]),]
  rownames(expr) <- expr[,1] 
  expr <- expr[,-1]
  geneNM <- get(load("data/Data_used_for_DataClean/geneNM.Rdata")) #data storing different types of gene names
  geneNM <- geneNM[!duplicated(geneNM[,toWhichType]) & !is.na(geneNM[,toWhichType]),]
  if (toWhichType == "Symbol"){
    geneNM <- geneNM[!duplicated(geneNM[,"Ensembl"]) & !is.na(geneNM[,"Ensembl"]),]
    rownames(geneNM) <- geneNM[,"Ensembl"]
    temp <- merge(expr, geneNM, by=0, all.x = T)
    temp <- temp[!duplicated(temp[,"Ensembl"]) & !is.na(temp[,"Ensembl"]),]
  } else {
    geneNM <- geneNM[!duplicated(geneNM[,"Symbol"]) & !is.na(geneNM[,"Symbol"]),]
    rownames(geneNM) <- geneNM[,"Symbol"]
    temp <- merge(expr, geneNM, by=0, all.x = T)
    temp <- temp[!duplicated(temp[,"Symbol"]) & !is.na(temp[,"Symbol"]),]
  }
  temp <- temp[!duplicated(temp[,toWhichType]) & !is.na(temp[,toWhichType]),]
  rownames(temp) <- temp[,toWhichType]
  expr <- temp[,!colnames(temp)%in%c("Row.names","Symbol","Ensembl","Biotype")]
  expr
}


#' Convert RNAseq data format from counTheoph#' Convert RNAseq data format from count to TPM ####
#'
#' @param expr Gene Expression Matrix, rownames are genes, colnames are samples, without log-trans, gene names are symbols
#' @return expression matrix of TPM format
#' @examples
#' expr_TPM <- CounttoTPM(expr_count)
CounttoTPM <- function(expr){
  length <- get(load("data/Data_used_for_DataClean/GeneLength.Rdata")) # data storing length of genes
  common <- intersect(rownames(length),rownames(expr))
  a <- expr[common, ]
  b <- length[common,]
  b <- as.numeric(b[,2])
  tpm <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums(x)))
  }
  c <- as.data.frame(tpm(a,b))
}


#' Perform tsne analysis for query expression matrix and ref expression matrix ####
#' 
#'@param datarRef expression matrix as reference, rownames are symbol, colnames are samples
#'@param dataset expression matrix as query, rownames are symbol, colnames are samples
#'@return a tsne plot to show the batch between two expression matrix
#'@examples
#'Mytsne(dataset,dataRef)
Mytsne <- function(dataset, dataRef) {
  require(M3C)
  genesInCommon <- intersect(rownames(dataset), rownames(dataRef))
  dataset <- dataset[rownames(dataset) %in% genesInCommon, ]
  dataRef <- dataRef[rownames(dataRef) %in% genesInCommon, ]
  batch <- c(rep(1, length(dataset[1, ])), rep(2, length(dataRef[1, ])))
  dataset <- cbind(dataset, dataRef)
  dataset <- unique(dataset)
  dataset <- dataset[complete.cases(dataset),]
  tsne(dataset,
       dotsize = 1.5,
       labels = as.factor(batch))
}


#' Remove batch effects of dataset with the reference of dataRef ####
#' 
#'@param datarRef expression matrix as reference, rownames are symbol, colnames are samples
#'@param dataset expression matrix as query, rownames are symbol, colnames are samples
#'@return a dataset (gene expression matrix) without batch effect compare with dataRef
#'@examples
#'dataset_RemovedBatchSize <- batchEffectRemoval(dataset_WithBatchEffect, dataRef)
BatchEffectRemoval <- function(dataset, dataRef) {
  require(sva)
  genesInCommon <-intersect(rownames(dataset), rownames(dataRef))
  dataset <- dataset[rownames(dataset) %in% genesInCommon, ]
  dataRef <- dataRef[rownames(dataRef) %in% genesInCommon, ]
  
  ref_rows <- as.character(rownames(dataRef))
  query_rows <- as.character(rownames(dataset))
  query_match <- match(query_rows, ref_rows)
  dataset <- dataset[order(query_match), ]
  
  batch <- c(rep(1, length(dataset[1, ])), rep(2, length(dataRef[1, ])))
  temp = ComBat(cbind(dataset, dataRef),
                         batch = batch,
                         ref.batch = 2)
  dataset <- temp[, 1:length(dataset[1, ])]
  dataset <- as.data.frame(dataset)
  dataset
}

#' z-score tranlation for gene expression matrix ####
#' @param expr Gene Expression Matrix, rownames are genes, colnames are samples, without log-trans, gene names are symbols
#' @return z-score translated expression matrix
#' @examples
#' expr_zscore <- get_z_score(expr_log2TPM)
Get_z_score <- function(expr){
  cal_z_score <- function(x) {
    #z score conversion
    x <- (x - mean(x)) / sd(x)
  }
  expr <- expr[rowSums(expr) != 0,] #delete all-zere rows
  expr <- expr[apply(expr, 1, sd) != 0, ] #delete rows whose sd=0
  expr <- as.matrix(expr)
  expr <- apply(expr,1,cal_z_score)
  expr <- as.data.frame(t(expr))
  expr
}

# Calculate overall survival form TCGA clinical data downloaded from TCGAbiolinks ####
#' 
#' @param expr Gene Expression Matrix, rownames are genes, colnames are samples, without log-trans, gene names are symbols
#' @return z-score translated expression matrix
#' @examples
#' expr_zscore <- get_z_score(expr_log2TPM)
CalculateOSfromTCGA <- function(Cli_GBM_fromTCGAbiolinks) {
  Cli_GBM_fromTCGAbiolinks$deceased <- Cli_GBM_fromTCGAbiolinks$vital_status == "Dead"
  Cli_GBM_fromTCGAbiolinks$overall_survival <- ifelse(Cli_GBM_fromTCGAbiolinks$deceased,
                                                      Cli_GBM_fromTCGAbiolinks$days_to_death,
                                     Cli_GBM_fromTCGAbiolinks$days_to_last_follow_up)
  Cli_GBM_fromTCGAbiolinks
}

# Unify the number of columns or rows between different cohorts in a list ####
#' 
#' @param unintersectList a list with different columns or rows between cohorts
#' @param  RowOrCol choose from "Row" and "Col", choose "Col" means unifying columns between cohorts
#' @return a list with same columns or rows between cohorts
#' @examples
#' intersectList <- GetIntersectList(unintersectList, "Col")
GetIntersectList <- function(unintersectList, RowOrCol){
  if (RowOrCol == "Row"){
    common_rows <- rownames(unintersectList[[1]])
    for (i in 2:length(unintersectList)) {
      common_rows <- intersect(common_rows, rownames(unintersectList[[i]]))
    } 
    filter_list <- list()
    for (i in 1:length(unintersectList)) {
      filter_list[[i]] <- unintersectList[[i]][common_rows, ]
      names(filter_list)[i] <- names(unintersectList)[i]
    } 
    filter_list
  } else {
    common_cols <- colnames(unintersectList[[1]])
    for (i in 2:length(unintersectList)) {
      common_cols <- intersect(common_cols, colnames(unintersectList[[i]]))
    } 
    filter_list <- list() 
    for (i in 1:length(unintersectList)) {
      filter_list[[i]] <- unintersectList[[i]][,common_cols]
      names(filter_list)[i] <- names(unintersectList)[i]
    } 
    filter_list
  }
  
}

