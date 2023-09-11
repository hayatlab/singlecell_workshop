library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

seurat <- readRDS("Visium_Heart_deconvoluted.rds")

#Show the seurat object
seurat

#Display the meta data
head(seurat@meta.data)

#Plot some QC and stats
VlnPlot(seurat, features = c("n_genes_by_counts", "total_counts"), pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(seurat, features = c("total_counts"))
SpatialFeaturePlot(seurat, features = c("n_genes_by_counts"))
SpatialFeaturePlot(seurat, features = c("ACTC1","DCN"))
SpatialFeaturePlot(seurat, features = c("Fibroblast","Cardiomyocytes"))

#Normalize the data
seurat <- SCTransform(seurat, assay = "Spatial", verbose = FALSE)

seurat <- RunPCA(seurat, assay = "SCT", verbose = FALSE)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
seurat <- FindClusters(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)

DimPlot(seurat, reduction = "umap", label = TRUE)
SpatialDimPlot(seurat, label = TRUE, label.size = 3)

