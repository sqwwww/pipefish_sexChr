library(Seurat)
library(tidyverse)

cdata <- readRDS("mid_data/Mbr_testis/output/v2/middleFiles/rds/combined_HaCinMbr_Zbfmo5_13000.v2_250603.rds")
DimPlot(cdata)


zbf_dataT <- readRDS("mid_data/published/zbf_Frontier_2022/publishedAnno/zebrafish_testis_10XscRNASeq_with_cell_type_annotation_and_umap.modifiedAnno.rds")
Ha_dataT <- readRDS("mid_data/Ha_testis/v9_removeAlleleBias/rds/Ha_afterAnno_convert_zbfID.v6Norm_rmBias_annotated.rds")
Cin_dataT <- readRDS("mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v9_annotated.rds")
Mbr_dataT <- readRDS("mid_data/Mbr_testis/output/v3_removeAlleleBias/rds/Mbr_afterAnno1_zbfconvert_v6NormFinal_allCells_rmBias.annotated.rds")

data.combined <- cdata

data.combined$cell_type3 <- NA

add_cluster <- function(obj, seu, cluster_col) {
  x <- as.character(seu[[cluster_col]][, 1])
  names(x) <- colnames(seu)
  cells <- intersect(names(x), colnames(obj))
  obj$cell_type3[cells] <- x[cells]
  return(obj)
}

data.combined <- add_cluster(data.combined, zbf_dataT, "zbf_cluster")
data.combined <- add_cluster(data.combined, Ha_dataT,  "Hab_cluster")
data.combined <- add_cluster(data.combined, Cin_dataT, "Cin_cluster")
data.combined <- add_cluster(data.combined, Mbr_dataT, "Mbr_cluster")

Idents(data.combined) <- "cell_type3"
table(Idents(data.combined), useNA = "ifany")

cells.use <- colnames(data.combined)[!is.na(Idents(data.combined))]

data.combined2 <- subset(data.combined, cells = cells.use)

Idents(data.combined2) <- data.combined2$cell_type3
Idents(data.combined) %>% table()
my_cols <- c(
  'spermatogonia'='#0067a5',
  'spermatocyte'='#f38400',
  'early spermatid' = "#8db600",
  'late spermatid' = '#008856',
  'sertoli' = '#e68fac',
  'leydig' = '#b3446c',
  'immune' = '#654522',
  'erythroid' = "#f3c300",
  'other' = '#bebebe')

combined_umap <- DimPlot(data.combined2, label = F,cols = my_cols)+
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )
combined_umap

combined_umap <- DimPlot(data.combined2, label = F,cols = my_cols)+
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

p3 <- DimPlot(data.combined2, reduction = "umap",split.by = "species",
              cols = my_cols,ncol = 2)+
  coord_fixed() +  # 保持 x/y 比例为1:1，使图形为正方形
  theme_minimal() + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加方框
    aspect.ratio = 1,  # 控制为正方形比例（部分 ggplot2 版本可能不生效，可配合 coord_fixed）
    strip.background = element_rect(fill = "gray90", color = "black"),  # 分面标签背景加边框
    panel.grid = element_blank(),  # 可选：去除网格
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
library(patchwork)
all_p <- combined_umap | p3
all_p

ggsave(filename = "plots/combined/final/allSpecies.v4.umap.combined.pdf", plot = all_p, width = 12, height = 5)


# output pseudobulk
zbf_dataT2 <- readRDS(file = "mid_data/published/zbf_Frontier_2022/Lawson_anno/annotation/zbf_mo5_ALLcells_Pei.v5_final.rds")
head(zbf_dataT2)
zbf_dataT2$zbf_clusterFinal <- zbf_dataT$zbf_cluster


pseudo_Mbr <- AverageExpression(Mbr_dataT, assays = "RNA",group.by = c("Mbr_cluster"))
pseudo_Mbr_data <- pseudo_Mbr$RNA
write.table(x = pseudo_Mbr_data, file = "mid_data/combined/pseudobulk/Mbr_rmBias_pseudobulk.v2.tsv", quote = F, row.names = T, col.names = T, sep = "\t")
pseudo_Hab <- AverageExpression(Ha_dataT, assays = "RNA",group.by = c("Hab_cluster"))
pseudo_Hab_data <- pseudo_Hab$RNA
write.table(x = pseudo_Hab_data, file = "mid_data/combined/pseudobulk/Hab_rmBias_pseudobulk.v2.tsv", quote = F, row.names = T, col.names = T, sep = "\t")
pseudo_Cin <- AverageExpression(Cin_dataT, assays = "RNA",group.by = c("Cin_cluster"))
pseudo_Cin_data <- pseudo_Cin$RNA
write.table(x = pseudo_Cin_data, file = "mid_data/combined/pseudobulk/Cin_rmBias_pseudobulk.v3.tsv", quote = F, row.names = T, col.names = T, sep = "\t")
pseudo_zbf <- AverageExpression(zbf_dataT2, assays = "RNA",group.by = c("zbf_clusterFinal"))
pseudo_zbf_data <- pseudo_zbf$RNA
write.table(x = pseudo_zbf_data, file = "mid_data/combined/pseudobulk/zbf_rmBias_pseudobulk.v2.tsv", quote = F, row.names = T, col.names = T, sep = "\t")
