# Download fails
# brain_data1<-read.csv("C:/Users/rodic/OneDrive/Документы/tissue_positions_list_1.csv",header=TRUE)

library (Lahman)
library(tidyverse)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
brain_data2<-read_rds("C:/Users/rodic/OneDrive/Документы/Human_Breast_Cancer_Whole_Transcriptome_Analysis_10xvisium.rds")


# Extract the dimension of the active assay using dim, nrow, or ncol
dim(x = brain_data2) # the number of features (genes) by samples (spots)
nrow(x = brain_data2) # the number of features
ncol(x = brain_data2) # the number of samples


# Extract the feature and sample names using rownames or colnames.
head(x = rownames(brain_data2), n = 10)
tail(x = colnames(brain_data2), n = 10)

# # Первый столбец как название 
# row.names(brain_data) <- make.names(brain_data[,1],TRUE)

# # in SeuratProject
 library(Seurat)
# CreateSeuratObject(
#   brain_data,
#   project = "SeuratProject",
#   assay = "RNA",
#   names.field = 1,
#   names.delim = "_",
#   meta.data = NULL,
# )
class(brain_data2[])

 
# # download Image
# InputImage <- Read10X_Image(
#   "C:/Users/rodic/OneDrive/Документы/spatial",
#   image.name = "tissue_hires_image.png",
#   filter.matrix = FALSE
# )

# # download filtered matrix
# H2_1_Seurat <- Load10X_Spatial(
#   "C:/Users/rodic/OneDrive/Документы",
#   filename = "filtered_feature_bc_matrix.h5",
#   assay = "Spatial",
#   image = InputImage
# )
# summary(brain_data)
#                          meta.data


 
# # Переименовать столбец
# colnames(brain_data[])
# names(brain_data)[names(brain_data)=="X3910"] <- "nCount_Spatial"
# names(brain_data)[names(brain_data)=="X3159"] <- 
# "nFeature_Spatial"
# install.packages("installr")
# library(installr)
# updateR()
# head(brain_data@meta.data)
colnames(brain_data2[]) # automatically calculated while creating the Seurat object

# BiocManager::install('Seurat', update = F, ask = F)
library(Seurat)
# 
# brain_data_1 <- as.data.frame(t(brain_data1))
# names(brain_data_1) <- brain_data_1[1,]
# brain_data_1 <- brain_data_1[-1,]
# pbmc<- CreateSeuratObject(counts = brain_data_1, project = "RNA")
 
brain_data <-brain_data2@meta.data
head(brain_data)
brain_data$nCount_Spatial[1:3]

# nFeature_Spatial: the number of unique genes in each sample
sum(brain_data$nFeature_Spatial ==  colSums(brain_data@assays$Spatial@counts > 0))
# summary(brain_data$nFeature_RNA)

# nCount_Spatial: the total number of detected molecules in each sample
# sum(brain_data$nCount_RNA ==  colSums(brain_data@assays$Spatial@counts))
sum(brain_data$nFeature_Spatial ==  colSums(brain_data@assays$Spatial@counts > 0))
# summary(brain_data$nCount_RNA)

# A vector of names of associated objects can be had with the names function
# These can be passed to the double [[ extract operator to pull them from the Seurat object
names(x = brain_data2)


brain_data2[['Spatial']] # equivalent to: brain_data@assays$Spatial

brain_data2[['image']] # equivalent to: brain_data@images$slice1
brain_data2@assays$Spatial@counts[5:10, 1:3] 

brain_data[['Spatial']]@meta.features
head(brain_data[['Spatial']][[]])

# The number of unique genes detected in each sample (nFeature_Spatial).
# The total number of molecules detected within a sample (nCount_Spatial).
# The percentage of reads that map to the mitochondrial genome.
brain_data2[["percent.mt"]] <- PercentageFeatureSet(brain_data2, pattern = "^mt-")

brain_data2[["percent.rb"]] <- PercentageFeatureSet(brain_data2, pattern = "^RP[SL]")

# Способ 1
VlnPlot(
  brain_data2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
  pt.size = 0.1, ncol = 3) &
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Способ 2
VlnPlot(brain_data2, features = c("nFeature_Spatial","nCount_Spatial","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))


# Jointly (rather than separately) consider the QC metrics when filtering
plot1 <- FeatureScatter(
  brain_data2, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(
  brain_data2, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
  NoLegend()
plot1 + plot2


 brain_subset <- subset(
  brain_data2, 
  subset = nFeature_Spatial < 8000 & nFeature_Spatial > 1000 &     nCount_Spatial < 50000 & percent.mt < 30)
 

# Normalization
SpatialFeaturePlot(
  brain_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt") &
  theme(legend.position = "bottom"))  

library (Seurat)
brain_norm <- SCTransform(brain_subset, assay = "Spatial", verbose = FALSE)
# (много весит и долго)

names(brain_norm)

dim(brain_norm@assays$SCT@counts) 
dim(brain_norm@assays$SCT@data) 
dim(brain_norm@assays$SCT@scale.data) 


# Downstream tasks (dimension reduction, data visualization, cluster annotation, differential expression)

brain_obj <- RunPCA(brain_norm, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
brain_obj <- FindNeighbors(brain_obj, reduction = "pca", dims = 1:30)
# Leiden algorithm for community detection
brain_obj <- FindClusters(brain_obj, verbose = FALSE)
# PCA result is the default UMAP input, use dimensions 1:30 as input features
brain_obj <- RunUMAP(brain_obj, reduction = "pca", dims = 1:30)

plot3 <- DimPlot(brain_obj, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(brain_obj, label = TRUE, label.size = 3) + NoLegend()
plot3 + plot4


# identity class of each sample
table(brain_obj@active.ident)

# find all markers of cluster 10
cluster10_markers <- FindMarkers(brain_obj, ident.1 = 9, min.pct = 0.25)
head(cluster10_markers, n = 5)


VlnPlot(brain_obj, features = c("SLC7A21", "CPB1", "PHGR1"))

SpatialFeaturePlot(object = brain_obj, 
                   features = rownames(cluster10_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)


# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
# this code chunk is not evaluated for now because of time constraints
brain_obj_markers <- FindAllMarkers(brain_obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
brain_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
top10 <- brain_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(brain_obj, features = top10$gene) + NoLegend()

# Moran’s I

brain_moransi <- FindSpatiallyVariableFeatures(
  brain_obj, assay = "SCT", 
  features = VariableFeatures(brain_obj)[1:10],
  selection.method = "moransi") 


moransi_output_df <- brain_moransi@assays$SCT@meta.features %>%
  na.exclude
head(moransi_output_df[order(moransi_output_df$MoransI_observed, decreasing = T), ])

library(ggplot2)
library(forecast)
top_features_moransi <- head(
  SpatiallyVariableFeatures(brain_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(brain_moransi, 
                   features = top_features_moransi, ncol = 3, alpha = c(0.1, 1)) + plot_annotation(
    title = "Top 3 genes with the largest Moran's I",
    subtitle = "among 10 top variable genes for illustration purposes")
 


bottom_features_moransi <- tail(
  SpatiallyVariableFeatures(brain_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(brain_moransi, 
                   features = bottom_features_moransi, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "Bottom 3 genes with the smallest Moran's I",
    subtitle = "among 10 top variable genes for illustration purposes")


# marker-variogram


brain_variogram <- FindSpatiallyVariableFeatures(
  brain_obj, assay = "SCT", 
  features = VariableFeatures(brain_obj)[1:10],
  selection.method = "markvariogram")  
variogram_output_df <- brain_variogram@assays$SCT@meta.features %>%
  na.exclude # there are NA rows b/c we only calculated the variogram for 10 genes
head(variogram_output_df[order(variogram_output_df$r.metric.5), ])



top_features_variogram <- head(
  SpatiallyVariableFeatures(brain_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain_variogram, 
                   features = top_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "3 genes with the top spatially variable rank (by mark-variogram)",
    subtitle = "among 10 top variable genes for illustration purposes")


bottom_features_variogram <- tail(
  SpatiallyVariableFeatures(brain_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(brain_variogram, 
                   features = bottom_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
    title = "3 genes with the bottom spatially variale rank (by mark-variogram)",
    subtitle = "among 10 top variable genes for illustration purposes")



# Convert data formats between R and Python
# SeuratDisk vignette to convert between Seurat and AnnData: https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
# sceasy to convert single-cell datasets across different formats: https://github.com/cellgeni/sceasy
# reticulate vignette to call Python from R and convert data formats: https://github.com/cellgeni/sceasy
# 
# 


