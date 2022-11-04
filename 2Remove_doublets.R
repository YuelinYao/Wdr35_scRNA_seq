#####################################################################
#### this script is to remove doublet ############################### 
#### packages: seurat  doubletfinder#################################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################
#####################################################################
#### reference: https://github.com/chris-mcginnis-ucsf/DoubletFinder#
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder) 
library(Seurat)
#### load previously merged dataset
load("Merge_pre.Rdata")

Singlet_final<-NULL
sample=155
for (sample in (unique(Merge$sample))){
  print(sample)
  sample_subset<-Merge[,Merge$sample%in%sample]
  dim(sample_subset)
  
  ##pre-processing
  sample_subset <- NormalizeData(sample_subset)
  sample_subset <- FindVariableFeatures(sample_subset, selection.method = "vst", nfeatures = 2000)
  sample_subset <- ScaleData(sample_subset)
  sample_subset <- RunPCA(sample_subset)
  sample_subset <- RunUMAP(sample_subset, dims = 1:30)
  
  sample_subset<- FindNeighbors(sample_subset, dims = 1:30)
  sample_subset <- FindClusters(sample_subset, resolution = 0.1)
  sample_subset@meta.data
  
  ## pK Identification (no ground-truth) --------------------------------------------
  ## pK ~ This defines the PC neighborhood size used to compute pANN, expressed as a proportion of the merged real-artificial data. No default is set, as pK should be adjusted for each scRNA-seq dataset. 
  ## Optimal pK values should be estimated using the strategy described below.
  sweep.res.list <- paramSweep_v3(sample_subset, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn<- find.pK(sweep.stats)
  bcmvn<-bcmvn[order(bcmvn$MeanBC,decreasing = T),]
  pk<-as.numeric(bcmvn$pK[1])
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------
  nExp_poi <- round(0.1*nrow(sample_subset@meta.data)) #Assuming 10% doublet formation rate
  
  annotations<-sample_subset@meta.data[["seurat_clusters"]]
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies -----------------------
  ## pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data. Default is set to 25%
  ## DoubletFinder performance is largely invariant of pN value selection
  ## nExp ~ This defines the pANN threshold used to make final doublet/singlet predictions.
  ## This value can best be estimated from cell loading densities into the 10X/Drop-Seq device, and adjusted according to the estimated proportion of homotypic doublets.
  
  sample_subset<- doubletFinder_v3(sample_subset, PCs = 1:10, pN = 0.25, pK =  pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  reuse.pANN<-colnames(sample_subset@meta.data)[grepl("pANN_0.25",colnames(sample_subset@meta.data))]
  sample_subset <- doubletFinder_v3(sample_subset, PCs = 1:10, pN = 0.25, pK =  pk, nExp = nExp_poi.adj,reuse.pANN=reuse.pANN,sct=FALSE)
  
  Singlet<-rownames(sample_subset@meta.data)[sample_subset@meta.data[,13]=="Singlet"]
  Singlet_final<-c(Singlet_final,Singlet)
}




combined_filter<-Merge[,colnames(Merge)%in%Singlet_final]
#save(combined_filter,file="combined_filter.Rdata")
                        
