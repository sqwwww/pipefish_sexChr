
library(Seurat)
library(tidyverse)

Mbr_dataT = readRDS("Mbr_afterAnno1_zbfconvert_v5NormFinal_allCells_rmBias.rds")
published_data <- readRDS("zbf.allCell.rds")
Ha_dataT <- readRDS("Ha_afterAnno_convert_zbfID.v7Norm_rmBias.rds")
Cint_dataT <- readRDS("Cin_afterAnno1_zbfconvert_v9.rds")
published_data$species <- "Dre"
Ha_dataT$species <- "Hab"
Cint_dataT$species <- "Cin"
Mbr_dataT$species <- "Mbr"


seurat_list <- list(published_data, Ha_dataT, Cint_dataT, Mbr_dataT)

for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}



anchors <- FindIntegrationAnchors(object.list = seurat_list,dims = 1:30, k.anchor = 3)
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:30)

data.combined <- integrated_data
#data.combined <- NormalizeData(data.combined)
data.combined <- FindVariableFeatures(data.combined)
data.combined <- ScaleData(data.combined)
data.combined <- RunPCA(data.combined)
data.combined <- FindNeighbors(data.combined,dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 1)
data.combined <- RunUMAP(data.combined, dims = 1:10)

saveRDS(data.combined, "combined_HaCinMbr_Zbfmo5_13000.v2.rds")
