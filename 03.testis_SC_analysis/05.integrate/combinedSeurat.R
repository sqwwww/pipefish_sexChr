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
saveRDS(integrated_data, "combined_HaCinMbr13000_Zbfmo5_13000.v4.rds")
