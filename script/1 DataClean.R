source("R/Functions_DataClean.R")
# 1 Create reference data (TCGA-GBM) ########
clinical <- read.table("data/Data_used_for_DataClean/TCGA_survival.tsv",
                       row.names = 1,header = T)
colnames(clinical)[2] <- "barcode"
colnames(clinical)[1] <- "deceased"
colnames(clinical)[3] <- "overall_survival"
clinical <- clinical[!duplicated(clinical$barcode),]
clinical <- CliClean(clinical)

expr <- as.data.frame(get(load("data/Data_used_for_DataClean/TCGA_TPM.Rdata")))
expr <- ConvertGeneName(expr,"Symbol")
expr <- log2(expr+1)
colnames(expr) <- gsub("\\.","-",substr(colnames(expr),1,12))
# expr <- Get_z_score(expr)

survivaldataTCGA <- MergeCliExp(expr, clinical)
dataRef <- t(survivaldataTCGA[,4:ncol(survivaldataTCGA)]) #for removing batch

# 2 Data Cleaning (Mol Cell 2021 GSE147352)################################################
clinical <- as.data.frame(read.csv("data/Data_used_for_DataClean/GSE147352_survival.csv"))
colnames(clinical)[1] <- "barcode"
colnames(clinical)[6] <- "deceased"
colnames(clinical)[5] <- "overall_survival"
clinical <- CliClean(clinical)
clinical$overall_survival <- as.numeric(clinical$overall_survival)*365


expr <- get(load("data/Data_used_for_DataClean/GSE147352_Count.Rdata"))
expr <- ConvertGeneName(expr,"Symbol")
expr <- CounttoTPM(expr)
expr <- log2(expr+1)
# expr <- Get_z_score(expr)

survivaldataGSE147352 <- MergeCliExp(expr, clinical)

#remove batch
dataset <- t(survivaldataGSE147352[,4:ncol(survivaldataGSE147352)])
Mytsne(dataset, dataRef) # with batch effect
dataset <- BatchEffectRemoval(dataset, dataRef)
Mytsne(dataset, dataRef) # without batch effect

survivaldataGSE147352<- mergeCliExp(dataset, clinical)



