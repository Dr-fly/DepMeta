source("R/Functions_EssentialGenes.R")
# 1 Extract DepMap data for GBM from pan-cancer data
DepMapList <- get(load("data/Data_used_for_EssentialGenes/DepMapList.Rdata"))
metadata <- DepMapList[["metadata"]]
fmetadata_GBM <- metadata[metadata$OncotreeLineage=="CNS/Brain",] #164
fmetadata_GBM <- fmetadata_GBM[fmetadata$OncotreePrimaryDisease=="Diffuse Glioma",] #132
fmetadata_GBM <- fmetadata_GBM[fmetadata$OncotreeSubtype=="Glioblastoma"|fmetadata$OncotreeSubtype=="Glioblastoma Multiforme",] #98
modelNMs_GBM <- fmetadata_GBM$ModelID
GBMDepMapList<- GetDepMapList(modelNMs_GBM,fmetadata_GBM)
save(GBMDepMapList, file = "data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata")
# 2 Screening essential genes at multi-omics levels
GBMDepMapList <- get(load("data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata"))
GBMDepMapResList <- list(
  mut = GetEssentialGenes(GBMDepMapList,"mut"),
  CNV = GetEssentialGenes(GBMDepMapList,"CNV"),
  expr = GetEssentialGenes(GBMDepMapList,"expr")
  )
save(GBMDepMapResList, file = "data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata")
# 3 Visualization
GBMDepMapList <- get(load("data/Data_used_for_EssentialGenes/GBMDepMapList.Rdata"))
PlotMutBoxFromDepMapList(c("DNMT3A","BCR","PIK3R1","ARHGEF40","GLI2","PROM1"),GBMDepMapList,c("#C25757", "#3A68AE"))
PlotScatterFromDepMapList(c("LMNB2","E2F3","KIFC1","MCM3"),"CNV",GBMDepMapList,c("#3A68AE","#C25757"))
PlotScatterFromDepMapList(c("NFIX","CDK6"),"expr",GBMDepMapList,c("#3A68AE","#C25757"))
