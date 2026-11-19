library(Seurat)
library(tidyverse)


# public annotation
reported <- readRDS(file = "mid_data/published/zbf_Frontier_2022/publishedAnno/zebrafish_testis_10XscRNASeq_with_cell_type_annotation_and_umap.rds")
Idents(reported) <- reported$cell_type
DimPlot(reported, label = TRUE)

write.table(reported@meta.data, file = "mid_data/published/zbf_Frontier_2022/publishedAnno/zebrafish_testis_10XscRNASeq_with_cell_type_annotation_and_umap.rds_meta.table", quote = F, sep = "\t")

reported$zbf_cluster <- case_when(
  as.character(reported$cell_type) %in% c("early round spermatid", "middle round spermatid") ~ "early spermatid",
  as.character(reported$cell_type) %in% c("late round spermatid", "elongated spermatid") ~ "late spermatid",
  as.character(reported$cell_type) == "erythrocyte" ~ "erythroid",
  as.character(reported$cell_type) %in% c("SPG1", "SPG2") ~ "spermatogonia",
  TRUE ~ as.character(reported$cell_type)
)

table(reported$zbf_cluster)
Idents(reported) <- reported$zbf_cluster
DimPlot(reported, label = TRUE)

# sub-cluster of sertoli
data.sc <- reported
data.sc[["main.cluster"]] <- Idents(data.sc);
clusterToSubset <- "sertoli";
data.sub <- subset(data.sc, subset = main.cluster == clusterToSubset);
data.sub <- SCTransform(data.sub, vars.to.regress = "percent.mt", verbose = FALSE)
data.sub <- RunPCA(data.sub, features = VariableFeatures(object = data.sub))
data.sub <- RunUMAP(data.sub, dims = 1:15)
data.sub <- FindNeighbors(data.sub, dims = 1:15)
data.sub <- FindClusters(data.sub, resolution =1.2)
head(Idents(data.sub), 5)
DimPlot(data.sub, reduction = "umap",label = T)
FeaturePlot(data.sub, features = c("dmrt1"))
FeaturePlot(data.sub, features = c("gsdf"))
FeaturePlot(data.sub, features = c("pcna"))
FeaturePlot(data.sub, features = c("ptprc"))
#sub_allMarkers <- FindAllMarkers(data.sub, min.pct = 0.5, only.pos = T, logfc.threshold = 1)
new.cluster.ids.subData <- c("sertoli","sertoli","immune","sertoli","sertoli")
names(new.cluster.ids.subData) <- levels(data.sub@meta.data$seurat_clusters)
data.sub <- RenameIdents(data.sub, new.cluster.ids.subData)
data.sub$newName <- Idents(data.sub)
DimPlot(data.sub, reduction = "umap",label = T)
data.sc$new_cluster <- data.sc$main.cluster %>% as.character()
data.sc$new_cluster[data.sc$main.cluster == "sertoli"] <- data.sub$newName %>% as.character()
Idents(data.sc) <- data.sc$new_cluster
p1 <- DimPlot(data.sc, reduction = "umap",label = T, cols = my_cols)+
  theme(panel.border = element_rect(fill = NA, colour = "black"),legend.position = "none",
        aspect.ratio = 1,
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
p1

Idents(reported) <- data.sc$new_cluster
reported$zbf_cluster <- data.sc$new_cluster
DimPlot(reported)

FeaturePlot(reported, features = c("ptprc"), max.cutoff = 4)
FeaturePlot(reported, features = c("krt18a.1"), max.cutoff = 10)
FeaturePlot(reported, features = c("star"))
FeaturePlot(reported, features = c("gsdf"), max.cutoff = 10)


write.table(reported@meta.data, file = "mid_data/published/zbf_Frontier_2022/publishedAnno/zebrafish_testis_10XscRNASeq_with_cell_type_annotation_and_umap.rds_meta.v2.table", quote = F, sep = "\t")
saveRDS(reported, "mid_data/published/zbf_Frontier_2022/publishedAnno/zebrafish_testis_10XscRNASeq_with_cell_type_annotation_and_umap.modifiedAnno.rds")





# library(SeuratDisk)
# 
# 
# 
# DefaultAssay(reported) <- "RNA"
# 
# 
# 
# reported@meta.data[] <- lapply(reported@meta.data, function(x) {
#   if (is.factor(x)) as.character(x) else x
# })
# 
# SaveH5Seurat(
#   reported,
#   filename = "mid_data/published/zbf_Frontier_2022/publishedAnno/reported.h5Seurat",
#   overwrite = TRUE
# )
# 
# Convert(
#   "mid_data/published/zbf_Frontier_2022/publishedAnno/reported.h5Seurat",
#   dest = "h5ad",
#   overwrite = TRUE
# )
