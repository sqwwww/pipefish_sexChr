library(Seurat)
library(tidyverse)

# options(Seurat.object.assay.version = "v5")

dataT.data <- Read10X(data.dir = "mid_data/Cint_testis/v3/")
dataT <- CreateSeuratObject(counts = dataT.data, project = "Cint_testis", min.cells = 3, min.features = 200)
colnames(dataT) <- paste('Cin', colnames(x = dataT), sep = '_')

dataT[["percent.mt"]] <- PercentageFeatureSet(dataT, pattern = "^MT-")
dataT <- subset(dataT, percent.mt < 5)
dataT <- subset(dataT, subset = nFeature_RNA > 200 & nCount_RNA < quantile(nCount_RNA, probs = 0.90))
Cin_dataT <- dataT

# output filtered data ----
library(Matrix)
counts <- Cin_dataT@assays$RNA$counts
barcodes <- colnames(Cin_dataT)
gene_names <- rownames(Cin_dataT)
filtered_dir = 'mid_data/Cint_testis/v3/filteredCounts/filteredCounts_3000Cells/'

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
# run sctransform
#Cin_dataT <- SCTransform(Cin_dataT, vars.to.regress = "percent.mt", verbose = FALSE)

Cin_dataT <- NormalizeData(Cin_dataT)
Cin_dataT <- FindVariableFeatures(Cin_dataT, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Cin_dataT)
Cin_dataT <- ScaleData(Cin_dataT, features = all.genes)
Cin_dataT <- RunPCA(Cin_dataT, verbose = FALSE)
ElbowPlot(Cin_dataT)
Cin_dataT <- FindNeighbors(Cin_dataT, dims = 1:10, verbose = FALSE)
Cin_dataT <- FindClusters(Cin_dataT, verbose = FALSE,resolution = 0.7)
Cin_dataT <- RunUMAP(Cin_dataT, dims = 1:10, verbose = FALSE)

DimPlot(Cin_dataT, label = TRUE)
write.table(Cin_dataT@meta.data, file = "mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v8.rds.seuratCluster.meta.table", quote = F, sep = "\t")

# check umap in original annotation
# ori_CinData <- readRDS("mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v7_final.rds")
# table(colnames(Cin_dataT) %in% colnames(ori_CinData))
# Cin_dataT$tmp <- ori_CinData$new_cluster
# Idents(Cin_dataT) <- Cin_dataT$Cin_cellType
# my_cols = c('unknown'='#bdbdbd',"spermatogonia A" = "#1f77b4",
#             "spermatogonia B" = "#1f77b4", 'spermatogonia'='#1f77b4','spermatocytes'='#f38400','spermatocyte'='#f38400',
#             'sertoli'='#604e97','setoli'='#604e97',
#             'leydig'='#e68fac',
#             'early spermatids'='#8db600','late spermatids' ='#008856',"spermatid" = "#008856",'early spermatid'='#8db600','late spermatid' ='#008856',
#             "spermatoza" = "#2b3d26","spermatozoa" = "#2b3d26",
#             "leukocytes" = "#654522","immune" = "#654522",
#             "dmrt1+ cells"="#604e97","erythroid"="#f3c300")
# DimPlot(Cin_dataT, label = TRUE, cols = my_cols)
# p1 <- DimPlot(Cin_dataT, reduction = "umap",label = T, cols = my_cols)+
#   theme(panel.border = element_rect(fill = NA, colour = "black"),legend.position = "none",
#         aspect.ratio = 1,
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank())
# p1
Idents(Cin_dataT) <- Cin_dataT$tmp
DimPlot(Cin_dataT, label = T)
saveRDS(Cin_dataT, "mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v8.rds")

# ggsave(filename = "plots/testis_Cint/v3/annotated_umap.v5.pdf", width = 6, height = 6, plot = p1)



Cin_dataT <- readRDS("mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v8.rds")
my_cols = c('unknown'='#bdbdbd',"spermatogonia" = "#1f77b4",
            'spermatocyte'='#f38400',
            'sertoli'='#604e97','setoli'='#604e97',
            'leydig'='#e68fac',
            'early spermatid'='#8db600','late spermatid' ='#008856',
            "spermatozoa" = "#2b3d26",
            "immune" = "#654522",
            "erythroid"="#f3c300")






# spermatogonia, spermatocyte, early spermatid, late spermatid, spermatozoa,sertoli,leydig, immune

grep(pattern = "cyp11", x = Features(Cin_dataT), value = T)

FeaturePlot(Cin_dataT, features = c("ptprc"))
FeaturePlot(Cin_dataT, features = c("star"))
FeaturePlot(Cin_dataT, features = c("krt18a.1"))
FeaturePlot(Cin_dataT, features = c("gsdf"))
FeaturePlot(Cin_dataT, features = c("dmrt1"))
FeaturePlot(Cin_dataT, features = c("cyp11c1"))
FeaturePlot(Cin_dataT, features = c("hbbe2"))

Idents(Cin_dataT) <- Cin_dataT$seurat_clusters
DimPlot(Cin_dataT,label = TRUE)
new.cluster.ids <- c("spermatogonia", #0
                     "immune", #1
                     "late spermatid", #2
                     "spermatogonia", #3
                     "spermatocyte", #4
                     "spermatogonia", #5
                     "sertoli", #6
                     "immune", #7 
                     'sertoli', #8
                     "early spermatid", #9
                     "early spermatid", #10
                     'late spermatid', #11
                     'leydig', #12
                     'spermatocyte', #13
                     'immune', #14
                     'sertoli' #15
)
names(new.cluster.ids) <- levels(Cin_dataT)
Cin_dataT <- RenameIdents(Cin_dataT, new.cluster.ids)
DimPlot(Cin_dataT, reduction = "umap", label = TRUE, pt.size = 0.5, cols = my_cols) + NoLegend()



# # sub-cluster of sertoli
# data.sc <- Cin_dataT
# data.sc[["main.cluster"]] <- Idents(data.sc);
# clusterToSubset <- "0";
# data.sub <- subset(data.sc, subset = main.cluster == clusterToSubset);
# data.sub <- SCTransform(data.sub, vars.to.regress = "percent.mt", verbose = FALSE)
# data.sub <- RunPCA(data.sub, features = VariableFeatures(object = data.sub))
# data.sub <- RunUMAP(data.sub, dims = 1:15)
# data.sub <- FindNeighbors(data.sub, dims = 1:15)
# data.sub <- FindClusters(data.sub, resolution =1.2)
# head(Idents(data.sub), 5)
# DimPlot(data.sub, reduction = "umap",label = T)
# FeaturePlot(data.sub, features = c("star"))
# FeaturePlot(data.sub, features = c("pcna"))
# FeaturePlot(data.sub, features = c("hbbe2"))
# #sub_allMarkers <- FindAllMarkers(data.sub, min.pct = 0.5, only.pos = T, logfc.threshold = 1)
# new.cluster.ids.subData <- c("spermatogonia","spermatogonia","spermatogonia","spermatogonia","spermatogonia","spermatogonia","spermatogonia","erythroid")
# names(new.cluster.ids.subData) <- levels(data.sub@meta.data$seurat_clusters)
# data.sub <- RenameIdents(data.sub, new.cluster.ids.subData)
# data.sub$newName <- Idents(data.sub)
# DimPlot(data.sub, reduction = "umap",label = T)
# data.sc$new_cluster <- data.sc$main.cluster %>% as.character()
# data.sc$new_cluster[data.sc$main.cluster == "0"] <- data.sub$newName %>% as.character()
# Idents(data.sc) <- data.sc$new_cluster
# p1 <- DimPlot(data.sc, reduction = "umap",label = T, cols = my_cols)+
#   theme(panel.border = element_rect(fill = NA, colour = "black"),legend.position = "none",
#         aspect.ratio = 1,
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank())
# p1

# ggsave(filename = "plots/testis_Cint/v3/annotated_umap.v6.pdf", plot = p1, width = 9, height = 7)
#Idents(Cin_dataT) <- data.sc$new_cluster
Cin_dataT$Cin_cluster <- Idents(Cin_dataT)
DimPlot(Cin_dataT)
write.table(Cin_dataT@meta.data, file = "mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v9.rds.cellType.meta.table", quote = F, sep = "\t")

saveRDS(Cin_dataT, "mid_data/Cint_testis/v3/rds/Cin_afterAnno1_zbfconvert_v9_annotated.rds")



# # sertoli cells
# FeaturePlot(Cin_dataT, features = c("dmrt1"))
# FeaturePlot(Cin_dataT, features = c("wt1b"))
# FeaturePlot(Cin_dataT, features = c("wt1a"))
# # krt18a.1
# FeaturePlot(Cin_dataT, features = c("LOC130912002"))
# FeaturePlot(Cin_dataT, features = c("gsdf"))
# VlnPlot(Cin_dataT, features = c("gsdf"))
# 
# # leydig cell
# FeaturePlot(Cin_dataT, features = c("star"))
# FeaturePlot(dataT, features = c("gata4"))
# # cyp17a2
# FeaturePlot(dataT, features = c("LOC130928348"))
# # CYP11A1
# FeaturePlot(dataT, features = c("LOC130916105"))
# 
# # amh, two copy
# FeaturePlot(Cin_dataT, features = c("LOC130912665"))
# FeaturePlot(Cin_dataT, features = c("amh"))
# FeaturePlot(Cin_dataT, features = c("LOC130922043")) # amhr2
# 
# 
# # spermatogonia
# FeaturePlot(Cin_dataT, features = c("piwil1"))
# FeaturePlot(Cin_dataT, features = c("piwil2"))
# FeaturePlot(Cin_dataT, features = c("dnd1"))
# VlnPlot(Cin_dataT, features = c("dnd1"))
# FeaturePlot(Cin_dataT, features = c("pcna"))
# FeaturePlot(Cin_dataT, features = c("dazl"))
# # spermatocytes
# FeaturePlot(Cin_dataT, features = c("dmc1"))
# FeaturePlot(Cin_dataT, features = c("syce2"))
# FeaturePlot(Cin_dataT, features = c("syce3"))
# VlnPlot(Cin_dataT, features = c("sycp1"))
# FeaturePlot(Cin_dataT, features = c("sycp1"))
# FeaturePlot(Cin_dataT, features = c("meiob"))
# FeaturePlot(Cin_dataT, features = c("spo11"))
# VlnPlot(Cin_dataT, features = c("spo11"))
# FeaturePlot(Cin_dataT, features = c("ccnb3"))
# FeaturePlot(Cin_dataT, features = c("sycp2"))
# FeaturePlot(Cin_dataT, features = c("dmc1"))
# # test
# # efcab2 reported in spermatocytes and elongate spermatids
# FeaturePlot(dataT, features = c("efcab2"))
# FeaturePlot(dataT, features = c("esr2a"))
# FeaturePlot(dataT, features = c("gata4"))
# FeaturePlot(dataT, features = c("LOC130908590"))
# FeaturePlot(dataT, features = c("podxl2"))
# 
# 
# # spermtids
# 
# FeaturePlot(dataT, features = c("spag8"))
# FeaturePlot(dataT, features = c("cfap206"))
# FeaturePlot(dataT, features = c("cfap251"))
# FeaturePlot(dataT, features = c("armc2"))
# VlnPlot(dataT, features = c("armc2"))
# 
# 
# # test
# FeaturePlot(dataT, features = c("cfap70"))
# FeaturePlot(dataT, features = c("LOC130930638")) # cfp99
# 
# cfap206
# FeaturePlot(dataT, features = c("dync2li1"))
# 
# # SMC -like
# FeaturePlot(dataT, features = c("tagln"))
# FeaturePlot(dataT, features = c("acta2"))
# FeaturePlot(dataT, features = c("dcn"))
# 
# FeaturePlot(Cin_dataT, features = c("hbbe2"))
# 
# # leukocytes 
# FeaturePlot(Cin_dataT, features = c("ptprc"))
# FeaturePlot(Cin_dataT, features = c("ptprb"))
# 
# FeaturePlot(dataT, features = c("cxcr4b"))
# FeaturePlot(dataT, features = c("itgb2"))
# FeaturePlot(dataT, features = c("lcp1"))
# 
# # leukocytes detail
# FeaturePlot(dataT, features = c("runx3"))
# FeaturePlot(dataT, features = c("mpx"))
# FeaturePlot(dataT, features = c("pax5"))
# FeaturePlot(dataT, features = c("pax5"))
# 
# 
# 
# 
