#####################################################################
#### This script is to identify marker for hc clustering  ###########
#### and computer p-value for clusters ##############################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################


library(Seurat)
library(sctransform)
library(xlsx)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
load("combined_filter.Rdata")


for (hc in c(3,4,5,6)){
  setwd(paste0("~/clustering_stability/",hc))
  print(hc)
  file=dir()
  file<-file[grepl(pattern = "^pc_10.Rdata", file)]
  file
  f="pc_10.Rdata"
  dir.create(paste0("~/cluster_marker/",hc))
  
  for (f in file){
    
    load(f)
    
    print(f)
    combined_filter@meta.data[["hierarchical_clustering"]]<-clu.out
    
    
    table(combined_filter$hierarchical_clustering)
    
    Idents(combined_filter)<-combined_filter$hierarchical_clustering
    for (i in unique(combined_filter$hierarchical_clustering))
      
    {
      print("========================")
      print(i)
      cluster_marker<-FindMarkers(combined_filter,ident.1 = i,min.diff.pct = 0.1)
      cluster_marker<-cluster_marker[cluster_marker$p_val_adj<0.05,]
      cluster_marker$abs<-abs(cluster_marker$avg_log2FC)
      cluster_marker<-cluster_marker[cluster_marker$abs>0.5,]
      cluster_marker$gene<-rownames(cluster_marker)
      cluster_marker$change<-NA
      cluster_marker$change[cluster_marker$avg_log2FC>0]<-"UP"
      cluster_marker$change[cluster_marker$avg_log2FC<0]<-"DOWN"
      
      write.csv(cluster_marker,file=paste0("~/cluster_marker/",hc,"/cluster_",i,".csv"))
      
      
      
      print("========================")
      
      
      
      for (cluster in unique(cluster_marker$change)){
        
        print(cluster)
        marker<-cluster_marker$gene[cluster_marker$change==cluster]
        marker
        gene<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                    filter="external_gene_name", values=marker, uniqueRows=TRUE)
        
        group_kegg <- enrichKEGG(gene = gene$entrezgene_id,organism = "mouse",
                                 pAdjustMethod = "BH",
                                 minGSSize = 1,
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.05)
        
        
        
        if (dim(group_kegg)[1]==0) {
          print("no result")
        } else {
          write.csv(group_kegg,file=paste0("~/cluster_marker/",hc,"/cluster_",i,"_",cluster,"_KEGG.csv"))
        }
        
        
        group_GO <- enrichGO(gene = gene$entrezgene_id,
                             OrgDb= org.Mm.eg.db,
                             ont = "BP",
                             pAdjustMethod = "BH",
                             minGSSize = 1,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE)
        
        if (dim(group_GO)[1]==0) {
          print("no result")
        } else {
          write.csv(group_GO,file=paste0("~/cluster_marker/",hc,"/cluster_",i,"_",cluster,"_BP_GO.csv"))
        }
        
        
        
        group_GO <- enrichGO(gene = gene$entrezgene_id,
                             OrgDb= org.Mm.eg.db,
                             ont = "CC",
                             pAdjustMethod = "BH",
                             minGSSize = 1,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE)
        
        if (dim(group_GO)[1]==0) {
          print("no result")
        } else {
          write.csv(group_GO,file=paste0("~/cluster_marker/",hc,"/cluster_",i,"_",cluster,"_CC_GO.csv"))
        }
        
        
        
        
        group_GO <- enrichGO(gene = gene$entrezgene_id,
                             OrgDb= org.Mm.eg.db,
                             ont = "MF",
                             pAdjustMethod = "BH",
                             minGSSize = 1,
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.05,
                             readable = TRUE)
        if (dim(group_GO)[1]==0) {
          print("no result")
        } else {
          write.csv(group_GO,file=paste0("~/cluster_marker/",hc,"/cluster_",i,"_",cluster,"_MF_GO.csv"))
        }
        
        
      }
      
    }
  }
}
