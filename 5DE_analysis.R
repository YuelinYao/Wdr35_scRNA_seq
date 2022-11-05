#####################################################################
#### This script is to perform DE analysis between cluster1/2 #######
#### main package: seurat ###########################################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################
#### ref: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#
#### hierarchical_clustering: 4 clusters, pc10 ######################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(xlsx)
load("combined_filter.Rdata")
load("background_gene.Rdata")


#######cluster1_cluster2
cluster1<-all_cell[combined_filter@meta.data[["hierarchical_clustering"]]=="1"]
cluster2<-all_cell[combined_filter@meta.data[["hierarchical_clustering"]]=="2"]

cluster1_cluster2_deg<-FindMarkers(combined_filter,ident.1 = cluster1,ident.2 = cluster2, min.pct = 0.25, logfc.threshold = 0.25)

logFC_cutoff=0.25
cluster1_cluster2_deg$DEG = as.factor(ifelse(cluster1_cluster2_deg$p_val_adj < 0.05 & abs(cluster1_cluster2_deg$avg_log2FC) > logFC_cutoff,
                                             ifelse(cluster1_cluster2_deg$avg_log2FC > logFC_cutoff ,'UP','DOWN'),'NOT'))



cluster1_cluster2_deg<-cluster1_cluster2_deg[!cluster1_cluster2_deg$DEG=="NOT",]

cluster1_cluster2_deg<-cluster1_cluster2_deg[order(cluster1_cluster2_deg$avg_log2FC),]


down_top<-cluster1_cluster2_deg$gene[1:20]
up_top<-tail(cluster1_cluster2_deg$gene,20)


pdf("~/heatmap_cluster1_cluster2_new.pdf",width=12,height=10)
p<-DoHeatmap(combined_filter,features =c(down_top,up_top),cells = c(cluster1,cluster2),group.by = "class")
print(p)
dev.off()


up_regulated_mut<-cluster1_cluster2_deg$gene[cluster1_cluster2_deg$DEG=="DOWN"]
up_regulated_wild<-cluster1_cluster2_deg$gene[cluster1_cluster2_deg$DEG=="UP"]


write.csv(up_regulated_mut,file = "~/up_regulated_gene_in_cluster2.csv")
write.csv(up_regulated_wild,file = "~/up_regulated_gene_in_cluster1.csv")


#### UMAP/tsne for DEGs
for (i in 1:length(up_regulated_mut)){
  if(up_regulated_mut[i]%in%rownames(combined_filter)){
    pdf(file = paste0("~/up_regulated_cluster2/umap/",up_regulated_mut[i],".pdf"),width=4,height=3)
    a<-FeaturePlot(combined_filter, features =up_regulated_mut[i],reduction = 'umap')
    print(a)
    dev.off()} else{
      print(paste0(up_regulated_mut[i]," are not founded"))
    }
}


for (i in 1:length(up_regulated_mut)){
  if(up_regulated_mut[i]%in%rownames(combined_filter)){
    pdf(file = paste0("~/up_regulated_cluster2/tsne/",up_regulated_mut[i],".pdf"),width=4,height=3)
    a<-FeaturePlot(combined_filter, features =up_regulated_mut[i],reduction = 'tsne')
    print(a)
    dev.off()} else{
      print(paste0(up_regulated_mut[i]," are not founded"))
    }
}


for (i in 1:length(up_regulated_wild)){
  if(up_regulated_wild[i]%in%rownames(combined_filter)){
    pdf(file = paste0("~/up_regulated_cluster1/umap/",up_regulated_wild[i],".pdf"),width=4,height=3)
    a<-FeaturePlot(combined_filter, features =up_regulated_wild[i],reduction = 'umap')
    print(a)
    dev.off()} else{
      print(paste0(up_regulated_wild[i]," are not founded"))
    }
}


for (i in 1:length(up_regulated_wild)){
  if(up_regulated_wild[i]%in%rownames(combined_filter)){
    pdf(file = paste0("~/up_regulated_cluster1/tsne/",up_regulated_wild[i],".pdf"),width=4,height=3)
    a<-FeaturePlot(combined_filter, features =up_regulated_wild[i],reduction = 'tsne')
    print(a)
    dev.off()} else{
      print(paste0(up_regulated_wild[i]," are not founded"))
    }
}




############### GO/KEGG analysis

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)

up_wt_mut_marker_genes<-up_regulated_wild
background_gene<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
                       filter="external_gene_name", values=background_gene, uniqueRows=TRUE)

gene<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
            filter="external_gene_name", values=up_wt_mut_marker_genes, uniqueRows=TRUE)

group_kegg <- enrichKEGG(gene = gene$entrezgene_id,universe = as.character(background_gene$entrezgene_id),organism = "mouse",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)

write.xlsx(group_kegg,file=paste0("~/wt_mut/up_wt_mut_KEGG.xlsx"),showNA=TRUE)



group_GO <- enrichGO(gene = gene$entrezgene_id,universe = as.character(background_gene$entrezgene_id),
                     OrgDb= org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)

write.xlsx(group_GO,file=paste0("~/wt_mut/up_wt_mut_GO.xlsx"),showNA=TRUE)

#####
up_mut_wt_marker_genes<-up_regulated_mut


gene<-getBM(mart=mart, attributes=c("external_gene_name","entrezgene_id"),
            filter="external_gene_name", values=up_mut_wt_marker_genes, uniqueRows=TRUE)


group_kegg <- enrichKEGG(gene = gene$entrezgene_id,universe = as.character(background_gene$entrezgene_id),organism = "mouse",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
                         
write.xlsx(group_kegg,file=paste0("~/wt_mut/up_mut_wt_KEGG.xlsx"),showNA=TRUE)

group_GO <- enrichGO(gene = gene$entrezgene_id,universe = as.character(background_gene$entrezgene_id),
                     OrgDb= org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
                     
write.xlsx(group_GO,file=paste0("~/wt_mut/up_mut_wt_GO.xlsx"),showNA=TRUE)







rm(list=ls())

list<-read.xlsx("~/result/Selected_GO.xlsx",sheetIndex = 1,header = F)
list<-list$X1
enrichment<-c(list,enrichment)
enrichment<-unique(enrichment)
enrichment
file<-dir("~/wt_mut",full.names = T)
file


condition<-substr(file,nchar("~/wt_mut/")+1,nchar(file)-nchar(".xlsx"))
condition
Group<-substr(condition,start = 1,stop = 10)
Group
Method<-substr(condition,start = 12,stop = 15)
Method

all_function<-data.frame(ID=NA,Description=NA,BgRatio=NA,GeneRatio=NA,p.adjust=NA,Group=NA,Method=NA,Count=NA)


for( i in 1:length(file)){
  in_go<-read.xlsx(file[i],header = T,stringsAsFactors = F,sheetIndex = 1)
  head(in_go)
  in_go<-in_go[in_go$ID%in%enrichment,]
  in_go<-in_go[,c(c(2,3,4,5,7,10))]
  in_go$Group<-Group[i]
  in_go$Method<-Method[i]
  in_go$Description=paste0(in_go$Method,":",in_go$Description)
  all_function<-rbind(all_function,in_go)
}


table(all_function$Group)
#bp + facet_grid(. ~ supp)
all_function<-all_function[-1,]
#######################################################
# PART 2: plot the representative GO terms by condition
# or list of DEGs and compare the enrichment with 
##########################################

# Prepare dataframe
#------------------
# Import the table containing the enriched GO terms by groups
all_function

# List objects and their structure contained in the dataframe 'GO_gp'
ls.str(all_function)

# Transform the column 'Gene_number' into a numeric variable
all_function$Gene_number <- all_function$Count

# Replace all the "_" by a space in the column containing the GO terms
#GO_gp$GO_biological_process <- chartr("_", " ", GO_gp$GO_biological_process)

# Transform the column 'GO_biological_process' into factors
#GO_gp$GO_biological_process<-as.factor(GO_gp$GO_biological_process)
all_function$Description<-as.factor(all_function$Description)
# Transform FDR values by -log10('FDR values')
all_function$'|log10(FDR)|' <- -(log10(all_function$p.adjust))

# Change factor order
condition
all_function$Group<- factor(all_function$Group,levels = c("up_wt_mut","up_mut_wt"))
all_function$Description<-factor(all_function$Description,levels=rev(levels(all_function$Description)))

# Create a vector with new names for groups to use in the plot
# Replace the terms by your own (\n allow to start a new line)
group.labs <- c(`up_wt_mut` = "up-regulated \nWt vs. Mut",
                `up_mut_wt` = "up-regulated \nMut vs. Wt")
group.labs

# Draw the plot with ggplot2
# to represent -log10(FDR) and Number of genes 
# of each GO biological process per group 
#---------------------------------------------------
ggplot(all_function, aes(x = Description, y = Group)) +
  geom_point(data=all_function,aes(x=Description, y=Group,size = Gene_number, colour = `|log10(FDR)|`), alpha=.7)+
  scale_y_discrete(labels =group.labs)+
  scale_color_gradient(low = "blue", high = "orange", limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        axis.title.x=element_blank())+
  xlab("Functional annotation")+
  labs(color="-log10(FDR)", size="Number\nof genes")




######################################
######extract cell id for RNA velocity
######################################
all_cells<-Cells(combined_filter)
all_id<-as.data.frame(t(as.data.frame(strsplit(all_cells, split = "-" ))))
all_id<-as.data.frame(t(as.data.frame(strsplit(all_id$V1, split = "_" ))))
all_id$sample<-NA
all_id$sample[all_id$V1=="155"]<-"K19Wdr3512m_155:"
all_id$sample[all_id$V1=="158"]<-"K19Wdr3512m_158:"
all_id$sample[all_id$V1=="186"]<-"K19Wdr3512m_186:"
all_id$sample[all_id$V1=="200"]<-"K19Wdr3512m_200:"
all_id$sample[all_id$V1=="168"]<-"K19Wdr3512m_168:"
all_id$sample[all_id$V1=="176"]<-"K19Wdr3512m_176:"
all_id$sample[all_id$V1=="177"]<-"K19Wdr3512m_177:"
all_id$sample<-paste0(all_id$sample,all_id$V2,"x")

all_id$sample

write.csv(all_id$sample,file="all_id.csv",row.names = FALSE)
write.csv(Embeddings(combined_filter, reduction = "umap"), file = "cell_embeddings.csv")
write.csv(combined_filter@meta.data[["hierarchical_clustering"]], file = "hierarchical_clustering_4.csv")

#color<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
#names(color)<-c(1,2,3,4)


