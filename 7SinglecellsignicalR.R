#####################################################################
#### this script is to perform singlecellsignalR  ###################
#### wild-type group: 155, 158, 176, 200 ############################
#### mutated group: 168, 177, 186 ###################################
###########Ligand/Receptor analysis using SingleCellSignalR

rm(list=ls())
setwd("~/result/")
load("combined_filter.Rdata")
library(Seurat)
#BiocManager::install("SingleCellSignalR")
library(SingleCellSignalR)

load("combined_filter.Rdata")
all_cell<-colnames(combined_filter)
wildtype<-all_cell[combined_filter@meta.data[["class"]]=="Wild type"]
mutanttype<-all_cell[combined_filter@meta.data[["class"]]=="Mutant type"]
all.genes <- rownames(combined_filter)


#within mutant animals only
mutant_expression<-combined_filter[,mutanttype]
data = data.frame(mutant_expression[["RNA"]]@data)

cluster_all= combined_filter@meta.data[["hierarchical_clustering"]]
names(cluster_all)<-all_cell
cluster<-cluster_all[mutanttype]
table(cluster)



signal = cell_signaling(data=data,genes=all.genes,cluster=cluster,species = "mus musculus",int.type = "autocrine")
inter.net <- inter_network(data = data, signal = signal, genes = all.genes, cluster =cluster, write = F)



cluster1_2<-signal[["cluster 1-cluster 2"]]


pdf("～/cluster1_cluster2.pdf",width = 6,height = 8)
visualize_interactions(signal=signal,show.in = 2,limit = 200)
dev.off()

pdf("～/cluster2_cluster1.pdf",width = 5,height = 8)
visualize_interactions(signal=signal,show.in = 5,limit = 200)
dev.off()


###a
cluster2_cluster1<-inter.net[["individual-networks"]][["cluster 2-cluster 1"]]


write.xlsx(cluster2_cluster1,"~/Within_mutant_samples_cluster2_cluster1.xlsx",row.names = F)


cluster1_cluster2<-inter.net[["individual-networks"]][["cluster 1-cluster 2"]]
cluster1_cluster2

write.xlsx(cluster1_cluster2,"~/Within_mutant_samples_cluster1_cluster2.xlsx",row.names = F)



cluster2<-read.csv("~/result/up_regulated_gene_in_cluster2.csv")
cluster1<-read.csv("~/result/up_regulated_gene_in_cluster1.csv")

###b
DEG<-c(cluster2$x,cluster1$x)
DEG
cluster2_cluster1

deg_cluster2_ligand<-cluster2_cluster1[cluster2_cluster1$ligand.name%in%DEG,]
deg_cluster2_ligand


write.xlsx(deg_cluster2_ligand,"~/Within_mutant_samples_cluster2(DEG_ligand)_cluster1_receptor.xlsx",row.names = F)

cluster1_cluster2_receptor<-cluster1_cluster2[cluster1_cluster2$receptor.name%in%DEG,]
dim(cluster1_cluster2_receptor)
write.xlsx(cluster1_cluster2_receptor,"~/Within_mutant_samples_cluster2(DEG_receptor)_cluster1_ligand.xlsx",row.names = F)

deg_cluster2_ligand<-deg_cluster2_ligand[,c(3,4,7,8)]
cluster1_cluster2_receptor<-cluster1_cluster2_receptor[,c(3,4,7,8)]

colnames(deg_cluster2_ligand)<-c("cluster 2","cluster 1", "interaction type","LRscore")
colnames(cluster1_cluster2_receptor)<-c("cluster 1","cluster 2", "interaction type","LRscore")


list_a<-list(deg_cluster2_ligand,cluster1_cluster2_receptor)
list_a
pdf("~/cluster 2 DEGs (as ligand) vs all genes for cluster 1 (as receptor).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_a,show.in = 1,limit = 200)
dev.off()

pdf("~/cluster 2 DEGs (as receptor) vs all genes for cluster 1 (as ligand).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_a,show.in = 2,limit = 200)
dev.off()



###c
deg_cluster1_ligand<-cluster1_cluster2[cluster1_cluster2$ligand.name%in%DEG,]
write.xlsx(deg_cluster1_ligand,"~/Within_mutant_samples_cluster1(DEG_ligand)_cluster2_receptor.xlsx",row.names = F)


deg_cluster1_receptor<-cluster2_cluster1[cluster2_cluster1$receptor.name%in%DEG,]
deg_cluster1_receptor
write.xlsx(deg_cluster1_receptor,"~/Within_mutant_samples_cluster1(DEG_receptor)_cluster2_ligand.xlsx",row.names = F)



deg_cluster1_receptor<-deg_cluster1_receptor[,c(3,4,7,8)]
deg_cluster1_ligand<-deg_cluster1_ligand[,c(3,4,7,8)]

colnames(deg_cluster1_receptor)<-c("cluster 2","cluster 1", "interaction type","LRscore")
colnames(deg_cluster1_ligand)<-c("cluster 1","cluster 2", "interaction type","LRscore")


list_c<-list(deg_cluster1_receptor,deg_cluster1_ligand)
list_c
pdf("~/Cluster 1 DEG (as receptor) vs all genes in cluster 2 (as ligand).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_c,show.in = 1,limit = 50)
dev.off()

pdf("~/cluster 1 DEGs (as ligand) vs all genes in cluster 2 (as receptor).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_c,show.in = 2,limit = 50)
dev.off()


###d

degcluster1_ligand_degcluster2_receptor<-cluster1_cluster2[cluster1_cluster2$receptor.name%in%DEG & cluster1_cluster2$ligand.name%in%DEG,]
degcluster1_ligand_degcluster2_receptor

write.xlsx(degcluster1_ligand_degcluster2_receptor,"~/Within_mutant_samples_cluster1(DEG_ligand)_cluster2(DEG_receptor).xlsx",row.names = F)


degcluster2_ligand_degcluster1_receptor<-cluster2_cluster1[cluster2_cluster1$receptor.name%in%DEG & cluster2_cluster1$ligand.name%in%DEG,]
degcluster2_ligand_degcluster1_receptor

write.xlsx(degcluster2_ligand_degcluster1_receptor,"~/Within_mutant_samples_cluster1(DEG_receptor)_cluster2(DEG_ligand).xlsx",row.names = F)




degcluster1_ligand_degcluster2_receptor<-degcluster1_ligand_degcluster2_receptor[,c(3,4,7,8)]
degcluster2_ligand_degcluster1_receptor<-degcluster2_ligand_degcluster1_receptor[,c(3,4,7,8)]

colnames(degcluster1_ligand_degcluster2_receptor)<-c("cluster 1","cluster 2", "interaction type","LRscore")
colnames(degcluster2_ligand_degcluster1_receptor)<-c("cluster 2","cluster 1", "interaction type","LRscore")


list_d<-list(degcluster1_ligand_degcluster2_receptor,degcluster2_ligand_degcluster1_receptor)
list_d
pdf("~/cluster 2 DEGs (as receptor) vs cluster 1 DEGs (as ligand).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_d,show.in = 1,limit = 50)
dev.off()

pdf("~/cluster 2 DEGs (as ligand) vs cluster 1 DEGs (as receptor).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_d,show.in = 2,limit = 50)
dev.off()




##e

cluster2_cluster2<-inter.net[["individual-networks"]][["cluster 2-cluster 2"]]
cluster2_cluster2

pdf("~/cluster2_cluster2_top125.pdf",width = 7,height = 5.5)
visualize_interactions(signal=signal,show.in = 6,limit = 125)
dev.off()





write.xlsx(cluster2_cluster2,"~/Within_mutant_samples_cluster2_cluster2.xlsx",row.names = F)


cluster2_DEG<-cluster2_cluster2[cluster2_cluster2$ligand.name%in%DEG & cluster2_cluster2$receptor.name%in%DEG,]
cluster2_DEG



cluster2_DEG<-cluster2_DEG[,c(3,4,7,8)]
colnames(cluster2_DEG)<-c("cluster 2","cluster 2", "interaction type","LRscore")

list_e<-list(cluster2_DEG,cluster2_DEG)
list_e
pdf("~/cluster 2 DEGs (as receptor) vs cluster 2 DEGs (as ligand).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list_e,show.in = 1,limit = 50)
dev.off()



write.xlsx(cluster2_DEG,"~/Within_mutant_samples_cluster2(DEG_receptor)_cluster2(DEG_ligand).xlsx",row.names = F)



####wild-type
wildtype_expression<-combined_filter[,wildtype]
wildtype_expression[1:10,1:10]
data = data.frame(wildtype_expression[["RNA"]]@data)

cluster_all= combined_filter@meta.data[["hierarchical_clustering"]]
names(cluster_all)<-all_cell
cluster<-cluster_all[wildtype]
table(cluster)
head(data[1:10,1:10])


signal = cell_signaling(data=data,genes=all.genes,cluster=cluster,species = "mus musculus",int.type = "autocrine")
signal

inter.net <- inter_network(data = data, signal = signal, genes = all.genes, cluster =cluster, write = F)
inter.net

cluster1_cluster1<-inter.net[["individual-networks"]][["cluster 1-cluster 1"]]
cluster1_cluster1

write.xlsx(cluster1_cluster1,"~/Within_wild_type_samples_cluster1_cluster1.xlsx",row.names = F)
cluster1_cluster1


pdf("~/cluster1_cluster1_top75.pdf",width = 5,height = 5.5)
visualize_interactions(signal=signal,show.in = 1,limit = 75)
dev.off()




cluster1_DEG<-cluster1_cluster1[cluster1_cluster1$ligand.name%in%DEG & cluster1_cluster1$receptor.name%in%DEG,]
cluster1_DEG
cluster1_DEG<-cluster1_DEG[,c(3,4,7,8)]


colnames(cluster1_DEG)<-c("cluster 1","cluster 1", "interaction type","LRscore")
list<-list(cluster1_DEG,cluster1_DEG)
list
pdf("~/All cluster 1 DEGs (autocrine).pdf",width = 5,height = 5.5)
visualize_interactions(signal=list,show.in = 1,limit = 50)
dev.off()




write.xlsx(cluster1_DEG,"~/Within_wild_type_samples_cluster1(DEG_receptor)_cluster1(DEG_ligand).xlsx",row.names = F)



