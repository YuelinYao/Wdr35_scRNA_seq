#####################################################################
#### this script is to perform hierachical clustering  ##############
#### and computer p-value for clusters ##############################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################
#### ref: clusterpval https://www.lucylgao.com/clusterpval/ #########
#### devtools::install_github("lucylgao/clusterpval")
library(dplyr)
library(Seurat)
library(ggplot2)
library(clusterpval)
library(fastcluster)
library(DropletUtils)
library(patchwork)
library(pheatmap)
source("util.R") #github("lucylgao/clusterpval")
load("combined_filter.Rdata")

pca<-combined_filter@reductions[["pca"]]@cell.embeddings

X<-as.matrix(pca)


cluster_number<-c(3,4,5,6)
dir.create(paste0("~/clustering_stability/"))
for (k in cluster_number){
  dir.create(paste0("~/clustering_stability/",k))
  pc_ndim<-c(5,10,15,20,25,30,35,40,45,50)
  for ( i in pc_ndim){
    X_var<-X[,1:i]
    print("starting hclstering...")
    hcl <- hclust(dist(X_var)^2, method="ward.D")
    print("done")
    combined_filter@meta.data[["hierarchical_clustering"]]<-cutree(hcl, k)
    clu.out<-cutree(hcl, k)
    save(clu.out,file=paste0("~/clustering_stability/",k,"/pc_",i,".Rdata"))
    pdf(file=paste0("~/clustering_stability/",k,"/pc_",i,".pdf"),width = 6,height = 5)
    p<-DimPlot(combined_filter, label = TRUE,  pt.size = 1, reduction = 'umap',group.by  = "hierarchical_clustering")+ggtitle(" ")
    print(p)
    dev.off()
    p_value<-array(data = NA,dim = c(k,k))
    p_value
    print("starting testing p value")
    for (p in 1:k){
      for (j in 1:k){
        if(j==p) {
          p_value[p,j]=1
        } else{
          result<-test_hier_clusters_exact(X_var, link="ward.D", K=k, k1=p, k2=j, hcl=hcl)
          p_value[p,j]<-result$pval
        }
      }}
    save(p_value,file=paste0("~/clustering_stability/",k,"/pc_p_value_",i,".Rdata"))
    print("done p value")
  }
}

#####################################
####### plot heatmap for p-value ####
#####################################

for (i in c(3:5)){
  setwd(paste0("~/clustering_stability/",i))
  file=dir()
  file<-file[grepl("p_value_pc", file)]
  file
  f="p_value_pc_10.Rdata"
  for (f in file){
    
    load(f)
    p_value
    MyFDR<-p.adjust(p_value,method = "BH")
    MyFDR_sig<-MyFDR
    MyFDR_sig[MyFDR_sig<=0.05] <- "*"
    MyFDR_sig[MyFDR_sig>0.05] <- " "
    MyFDR_sig <- matrix(MyFDR_sig,nrow = dim(p_value)[1],byrow = F)
    MyFDR_m<-matrix(MyFDR,nrow = dim(p_value)[1],byrow = F)
    MyFDR_m
    MyFDR<--log10(MyFDR_m)
    n=dim(p_value)[1]
    n
    rownames(MyFDR)<-c(1:n)
    colnames(MyFDR)<-c(1:n)
    library(pheatmap)
    MyFDR
    breaksList = seq(0, 10, by = 1)
    p_value
    f
    #display_numbers = MyFDR_sig
    name<-gsub("Rdata"," ",f)
    name
    pdf(file = paste0("~/clustering_stability/",i,"/","hc_",name,"pdf"))
    pheatmap(MyFDR,display_numbers = round(MyFDR,digits = 2),fontsize_number = 15,fontsize = 18,color = colorRampPalette(c("white","firebrick3"))(10),breaks = breaksList,cluster_rows = F,cluster_cols = F)
    
    dev.off()
    
  }
}







