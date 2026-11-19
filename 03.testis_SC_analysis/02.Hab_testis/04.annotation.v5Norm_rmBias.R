library(Seurat)
library(tidyverse)

dataT.data <- Read10X(data.dir = "mid_data/Ha_testis/v9_removeAlleleBias/cellRanger/")
dataT <- CreateSeuratObject(counts = dataT.data, project = "Ha_testis", min.cells = 3, min.features = 200)
colnames(dataT) <- paste('Hab', colnames(x = dataT), sep = '_')
dataT[["percent.mt"]] <- PercentageFeatureSet(dataT, pattern = "^MT-")
Hab_dataT <- dataT
Hab_dataT <- subset(dataT, subset = nCount_RNA >  5000 )
Hab_dataT <- subset(Hab_dataT, subset = nFeature_RNA >  200 & nFeature_RNA <  quantile(nFeature_RNA, probs = 0.9) & percent.mt < 5)


# output filtered data ----
library(Matrix)
counts <- Hab_dataT@assays$RNA$counts
barcodes <- colnames(Hab_dataT)
gene_names <- rownames(Hab_dataT)
filtered_dir = 'mid_data/Ha_testis/v9_removeAlleleBias/cellRanger/filteredCounts3000Cells/'

if (!dir.exists(filtered_dir)) {
  dir.create(filtered_dir, recursive = TRUE)
}
## Output counts matrix
writeMM(counts, paste0(filtered_dir, 'matrix.mtx'))

# Output cell barcodes
barcodes <- colnames(counts)
write.table(as.data.frame(barcodes), paste0(filtered_dir, 'barcodes.tsv'),
            col.names = FALSE, row.names = FALSE, sep = "\t")
# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,type = "Gene Expression")
write.table(as.data.frame(features), sep = "\t", paste0(filtered_dir, 'genes.tsv'),
            col.names = FALSE, row.names = FALSE)



# Normalizing the data
Hab_dataT <- NormalizeData(Hab_dataT)
Hab_dataT <- FindVariableFeatures(Hab_dataT, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Hab_dataT)
Hab_dataT <- ScaleData(Hab_dataT, features = all.genes)
Hab_dataT <- RunPCA(Hab_dataT, verbose = FALSE)
ElbowPlot(Hab_dataT)
Hab_dataT <- FindNeighbors(Hab_dataT, dims = 1:15, verbose = FALSE)
Hab_dataT <- FindClusters(Hab_dataT, verbose = FALSE,resolution = 0.7)
Hab_dataT <- RunUMAP(Hab_dataT, dims = 1:15, verbose = FALSE)

DimPlot(Hab_dataT, label = TRUE)
FeaturePlot(object = Hab_dataT, features = c("dnd1"))
FeaturePlot(object = Hab_dataT, features = c("nanos3"))
FeaturePlot(object = Hab_dataT, features = c("dazl"))
FeaturePlot(object = Hab_dataT, features = c("piwil1"))

FeaturePlot(object = Hab_dataT, features = c("spag8"))
FeaturePlot(object = Hab_dataT, features = c("spag6"))
FeaturePlot(object = Hab_dataT, features = c("larp1b"))
FeaturePlot(object = Hab_dataT, features = c("armc2"))

FeaturePlot(object = Hab_dataT, features = c("star"))
FeaturePlot(object = Hab_dataT, features = c("krt18a.1"))
FeaturePlot(object = Hab_dataT, features = c("gsdf"))
FeaturePlot(object = Hab_dataT, features = c("dmrt1"))
FeaturePlot(object = Hab_dataT, features = c("ptprc"))
FeaturePlot(object = Hab_dataT, features = c("dcn"))
FeaturePlot(object = Hab_dataT, features = c("hbbe2"))

# my_cols = c('unknown'='#bdbdbd',"spermatogonia" = "#1f77b4",
#             'spermatocyte'='#f38400',
#             'sertoli'='#604e97','setoli'='#604e97',
#             'leydig'='#e68fac',
#             'early spermatid'='#8db600','late spermatid' ='#008856',
#             "spermatozoa" = "#2b3d26",
#             "immune" = "#654522", 
#             "erythroid"="#f3c300")
# 
# # check umap in original annotation
# ori_HabData <- readRDS("mid_data/Ha_testis/v7anno/rds/Ha_afterAnno_convert_zbfID.v3.rds")
# table(colnames(Hab_dataT) %in% colnames(ori_HabData))
# Hab_dataT$Hab_cellType <- ori_HabData$main.cluster
# Idents(Hab_dataT) <- Hab_dataT$Hab_cellType
# DimPlot(Hab_dataT, label = TRUE, cols = my_cols)
# p1 <- DimPlot(Hab_dataT, reduction = "umap",label = T, cols = my_cols)+
#   theme(panel.border = element_rect(fill = NA, colour = "black"),legend.position = "none",
#         aspect.ratio = 1,
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank())
# p1

saveRDS(Hab_dataT, "mid_data/Ha_testis/v9_removeAlleleBias/rds/Ha_afterAnno_convert_zbfID.v6Norm_rmBias.rds")
write.table(Hab_dataT@meta.data, file = "mid_data/Ha_testis/v9_removeAlleleBias/rds/Ha_afterAnno_convert_zbfID.v6Norm_rmBias.rds.seuratCluster.meta.table", quote = F, sep = "\t")
# ggsave(filename = "plots/testis_Ha/v7_anno/annotated_umap.v5.pdf", width = 6, height = 6, plot = p1)

# annotation
new.cluster.ids <- c("spermatogonia", #0
                     "early spermatid", #1
                     "spermatogonia", #2
                     "spermatogonia", #3
                     "sertoli", #4
                     "spermatogonia", #5
                     "spermatocyte", #6
                     "spermatogonia", #7
                     'spermatogonia', #8
                     "spermatocyte", #9
                     "leydig", #10
                     'early spermatid', #11
                     'sertoli', #12
                     'leydig', #13
                     'late spermatid' #14
)
Idents(Hab_dataT) <- Hab_dataT$seurat_clusters
names(new.cluster.ids) <- levels(Hab_dataT@meta.data$seurat_clusters)
Hab_dataT <- RenameIdents(Hab_dataT, new.cluster.ids)
p1 <- DimPlot(Hab_dataT, reduction = "umap", label = TRUE)
Hab_dataT$Hab_cluster <- Idents(Hab_dataT)
head(Hab_dataT)
ggsave(filename = "plots/testis_Ha/v9_removeAlleleBias/umap.pdf", plot = p1, width = 9, height = 7)

write.table(Hab_dataT@meta.data, file = "mid_data/Ha_testis/v9_removeAlleleBias/rds/Ha_afterAnno_convert_zbfID.v6Norm_rmBias.rds.cellType.meta.table", quote = F, sep = "\t")
saveRDS(Hab_dataT, "mid_data/Ha_testis/v9_removeAlleleBias/rds/Ha_afterAnno_convert_zbfID.v6Norm_rmBias_annotated.rds")
