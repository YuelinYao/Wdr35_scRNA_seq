#####################################################################
#### this script is to read 10x scRNA-seq data from cellranger output 
#### packages: seurat ###############################################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################
#####################################################################

rm(list=ls())
library(Seurat)

sample="K19Wdr3512m_200"
{
  barcodes<-read.table(paste0(sample,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),header = F)
  barcodes
  
  pre<- Read10X(data.dir = paste0(sample,"/filtered_feature_bc_matrix/"))
  pre <- CreateSeuratObject(counts = pre, project = sample)
  
  id<-gsub("K19Wdr3512m_",replacement = "",sample)
  id
  table(colnames(pre)==barcodes$V1)
  pre$barcodes<-barcodes$V1
  pre@meta.data
  pre<-RenameCells(object = pre,new.names =paste0(id,"_",colnames(pre)) )
  pre$sample<-id
  
  pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^mt-")
  pre$log10GenesPerUMI <- log10(pre$nFeature_RNA) / log10(pre$nCount_RNA)
  
  
  pre <- subset(pre, cells = sample(Cells(pre), 1500,replace = F))  ###downsample to 1500
  
  counts <- GetAssayData(object = pre, slot = "counts")
  nonzero <- counts > 0
  index<-Matrix::rowSums(nonzero) >= 3  ##keep genes that are present in more than 3 cells in each sample.
  keep_genes <- rownames(nonzero)[index]
  
  
}

keep_gene_final<-keep_genes
Merge<-pre

for (sample in c("K19Wdr3512m_186","K19Wdr3512m_177","K19Wdr3512m_176","K19Wdr3512m_168","K19Wdr3512m_158","K19Wdr3512m_155"))
{
  print(sample)
  barcodes<-read.table(paste0(sample,"/filtered_feature_bc_matrix/barcodes.tsv.gz"),header = F)
  barcodes
  
  pre<- Read10X(data.dir = paste0(sample,"/filtered_feature_bc_matrix/"))
  pre <- CreateSeuratObject(counts = pre, project = sample)
  
  id<-gsub("K19Wdr3512m_",replacement = "",sample)
  id
  table(colnames(pre)==barcodes$V1)
  pre$barcodes<-barcodes$V1
  pre@meta.data
  pre<-RenameCells(object = pre,new.names =paste0(id,"_",colnames(pre)) )
  pre$sample<-id
  
  pre[["percent.mt"]] <- PercentageFeatureSet(pre, pattern = "^mt-")
  pre$log10GenesPerUMI <- log10(pre$nFeature_RNA) / log10(pre$nCount_RNA)

  Merge<-merge(x=pre,y=Merge) 
  
  counts <- GetAssayData(object = pre, slot = "counts")
  nonzero <- counts > 0
  index<-Matrix::rowSums(nonzero) >= 3
  keep_genes <- rownames(nonzero)[index]

  keep_gene_final<-intersect(keep_genes,keep_gene_final)
  
}



Merge<-Merge[keep_gene_final,]
Merge

Merge$sample
Merge$Group<-NA
Merge$Group[Merge$sample%in%c("155", "158", "176","200")] <- "Wild type"
Merge$Group[Merge$sample%in%c("168", "177","186")] <- "Mutant"



table(Merge$Group)
table(Merge$sample)
Merge <-subset(Merge, subset= nFeature_RNA > 595 & percent.mt < 5 & nFeature_RNA < 7000 & nCount_RNA < 50000 &Epcam>0 & Spp1>0)
#save(Merge,file="Merge_pre.Rdata")
