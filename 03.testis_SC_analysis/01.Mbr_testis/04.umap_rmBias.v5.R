library(Seurat)
library(tidyverse)

dataT.data <- Read10X(data.dir = "mid_data/Mbr_testis/output/v3_removeAlleleBias/",gene.column=2)
dataT <- CreateSeuratObject(counts = dataT.data, project = "Mbr_testis", min.cells = 3, min.features = 200)
colnames(dataT) <- paste('Mbr', colnames(x = dataT), sep = '_')
dataT
dataT[["percent.mt"]] <- PercentageFeatureSet(dataT, pattern = "^MT-")
# Mbr_dataT <- subset(dataT, subset = nCount_RNA >  4500 )
# Visualize QC metrics as a violin plot
Mbr_dataT <- dataT
VlnPlot(Mbr_dataT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Mbr_dataT <- subset(Mbr_dataT, subset = nFeature_RNA >  200 & nFeature_RNA <  quantile(nFeature_RNA, probs = 0.9) & percent.mt < 5)


# output filtered data ----
library(Matrix)
counts <- Mbr_dataT@assays$RNA$counts
barcodes <- colnames(Mbr_dataT)
gene_names <- rownames(Mbr_dataT)
filtered_dir = 'mid_data/Mbr_testis/output/v3_removeAlleleBias/filtered_count/'

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


# Normalizing the data ---
Mbr_dataT <- NormalizeData(Mbr_dataT)
Mbr_dataT <- FindVariableFeatures(Mbr_dataT, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Mbr_dataT)
Mbr_dataT <- ScaleData(Mbr_dataT, features = all.genes)
Mbr_dataT <- RunPCA(Mbr_dataT, verbose = FALSE)
ElbowPlot(Mbr_dataT)
Mbr_dataT <- FindNeighbors(Mbr_dataT, dims = 1:15, verbose = FALSE)
Mbr_dataT <- FindClusters(Mbr_dataT, verbose = FALSE,resolution = 0.7)
Mbr_dataT <- RunUMAP(Mbr_dataT, dims = 1:15, verbose = FALSE)

DimPlot(Mbr_dataT, label = TRUE)

FeaturePlot(Mbr_dataT, features = c("ptprc"))
FeaturePlot(Mbr_dataT, features = c("hbbe2"), max.cutoff = 1)
FeaturePlot(Mbr_dataT, features = c("alas2"))
FeaturePlot(Mbr_dataT, features = c("klf1"))
FeaturePlot(Mbr_dataT, features = c("sycp1"))
FeaturePlot(Mbr_dataT, features = c("sycp2"))
FeaturePlot(Mbr_dataT, features = c("meiob"))
FeaturePlot(Mbr_dataT, features = c("dmc1"))
FeaturePlot(Mbr_dataT, features = c("efcab2"))
FeaturePlot(Mbr_dataT, features = c("spag8"), min.cutoff = 1)
FeaturePlot(Mbr_dataT, features = c("spag6"))
FeaturePlot(Mbr_dataT, features = c("larp1b"))
FeaturePlot(Mbr_dataT, features = c("cfap70"))
FeaturePlot(Mbr_dataT, features = c("akap4"))
FeaturePlot(Mbr_dataT, features = c("star"))
FeaturePlot(Mbr_dataT, features = c("cyp11a1"))
FeaturePlot(Mbr_dataT, features = c("krt18a.1"))
FeaturePlot(Mbr_dataT, features = c("gsdf"))
FeaturePlot(Mbr_dataT, features = c("dmrt1"))
FeaturePlot(Mbr_dataT, features = c("cyp17a1"))
FeaturePlot(Mbr_dataT, features = c("itgb2"))

# # check umap in original annotation
# ori_MbrData <- readRDS("mid_data/Mbr_testis/output/v2/middleFiles/rds/Mbr_afterAnno1_zbfconvert_v2_all9000cells.rds")
# table(colnames(Mbr_dataT) %in% colnames(ori_MbrData))
# Mbr_dataT$Mbr_cellType <- ori_MbrData$new_cluster
# Idents(Mbr_dataT) <- Mbr_dataT$Mbr_cellType
# DimPlot(Mbr_dataT, label = TRUE, cols = my_cols)
# my_cols = c('unknown'='#bdbdbd',"spermatogonia" = "#1f77b4",
#             'spermatocyte'='#f38400',
#             'sertoli'='#604e97','setoli'='#604e97',
#             'leydig'='#e68fac',
#             'early spermatid'='#8db600','late spermatid' ='#008856',
#             "spermatozoa" = "#2b3d26",
#             "immune" = "#654522", 
#             "erythroid"="#f3c300")
# p1 <- DimPlot(Mbr_dataT, reduction = "umap",label = T, cols = my_cols)+
#   theme(panel.border = element_rect(fill = NA, colour = "black"),legend.position = "none",
#         aspect.ratio = 1,
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank())
# p1
saveRDS(Mbr_dataT, "mid_data/Mbr_testis/output/v3_removeAlleleBias/rds/Mbr_afterAnno1_zbfconvert_v6NormFinal_allCells_rmBias.rds")
write.table(Mbr_dataT@meta.data, file = "mid_data/Mbr_testis/output/v3_removeAlleleBias/rds/Mbr_afterAnno1_zbfconvert_v6NormFinal_allCells_rmBias.rds_meta.table", quote = F, sep = "\t")

# ggsave(filename = "plots/testis_Mbr/v2/annotated_umap.v3_all9000cells.pdf", width = 6, height = 6, plot = p1)

Idents(Mbr_dataT) <- Mbr_dataT$seurat_clusters
allMarkers <- FindAllMarkers(object = Mbr_dataT, logfc.threshold = 1, only.pos = T)
m_15 <- allMarkers %>%filter(cluster == 15)
m_15$gene %>% head(100) %>% paste0(collapse = ",")
m_19 <- allMarkers %>%filter(cluster == 19)
m_19$gene %>% head(100) %>% paste0(collapse = ",")
FeaturePlot(object = Mbr_dataT, features = c("scube3"))

DimPlot(Mbr_dataT,label = TRUE)
new.cluster.ids <- c("spermatogonia", #0
                     "early spermatid", #1
                     "early spermatid", #2
                     "spermatocyte", #3
                     "early spermatid", #4
                     "early spermatid", #5
                     "spermatocyte", #6
                     "early spermatid", #7 
                     'early spermatid', #8
                     "spermatogonia", #9
                     "spermatocyte", #10
                     'late spermatid', #11
                     'early spermatid', #12
                     'sertoli', #13
                     'late spermatid', #14
                     'other', #15
                     'immune', #16
                     'late spermatid', #17
                     'spermatogonia', #18
                     "other"
)
names(new.cluster.ids) <- levels(Mbr_dataT)
Mbr_dataT <- RenameIdents(Mbr_dataT, new.cluster.ids)
#DimPlot(Mbr_dataT, reduction = "umap", label = TRUE, pt.size = 0.5, cols = my_cols) + NoLegend()
Mbr_dataT$Mbr_cluster <- Idents(Mbr_dataT)
write.table(Mbr_dataT@meta.data, file = "mid_data/Mbr_testis/output/v3_removeAlleleBias/rds/Mbr_afterAnno1_zbfconvert_v6NormFinal_allCells_rmBias.rds.cellType.meta.table", quote = F, sep = "\t")
saveRDS(Mbr_dataT, "mid_data/Mbr_testis/output/v3_removeAlleleBias/rds/Mbr_afterAnno1_zbfconvert_v6NormFinal_allCells_rmBias.annotated.rds")


# tmp <- readRDS(file = "mid_data/Mbr_testis/output/v3_removeAlleleBias/rds/integrated_anno_Mbr13000.final.2.rds")
# Mbr_dataT$tmp <- tmp$combinedAnno
# Idents(Mbr_dataT) <- Mbr_dataT$tmp
# DimPlot(Mbr_dataT)
