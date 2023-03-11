#####################################################################
#### this script is to pre-processing scRNA-seq data ################
#### main package: seurat ###########################################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################
####ref: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html #

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


combined_filter <- NormalizeData(combined_filter, normalization.method = "LogNormalize", scale.factor = 10000)
#combined_filter <- FindVariableFeatures(combined_filter, selection.method = "vst", nfeatures = 5000)
combined_filter <- FindVariableFeatures(combined_filter, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(combined_filter)
variable_gene<-combined_filter@assays[["RNA"]]@var.features


combined_filter<-ScaleData(combined_filter,vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt") )     
#combined_filter<-ScaleData(combined_filter)
combined_filter  <- RunPCA(object = combined_filter,features = variable_gene)

## Determine the ‘dimensionality’ of the dataset 
combined_filter <- JackStraw(combined_filter, num.replicate = 100,dims = 50)
combined_filter <- ScoreJackStraw(combined_filter, dims = 1:50)
JackStrawPlot(combined_filter, dims = 1:50)
ElbowPlot(combined_filter,ndims = 50)



combined_filter <- RunUMAP(combined_filter,dims  = 1:10)
combined_filter <- RunUMAP(combined_filter,features  = all.genes)
combined_filter <- RunTSNE(combined_filter,dims  = 1:10)

combined_filter <- FindNeighbors(combined_filter,dims = 1:10)
combined_filter <- FindClusters(combined_filter, resolution = 0.1)


##UMAP & TSNE plot
DimPlot(combined_filter, label = TRUE,  pt.size = 1, reduction = 'umap',group.by = "Group")
DimPlot(combined_filter, label = TRUE,  pt.size = 1, reduction = 'tsne',group.by = "Group")

#save(combined_filter,file="combined_filter.Rdata")





