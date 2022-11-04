library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dbplyr)


meta_data<-Merge@meta.data
tapply(meta_data$nCount_RNA, meta_data$sample, summary)
tapply(meta_data$nFeature_RNA, meta_data$sample, summary)
tapply(meta_data$log10GenesPerUMI, meta_data$sample, summary)
tapply(meta_data$percent.mt, meta_data$sample, summary)



p<-VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","log10GenesPerUMI"), ncol = 4,group.by = "sample")
p& theme(axis.title.x = element_blank())




p<-meta_data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 50000,col="red")
p



p<-meta_data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + ylab("log10 Cell density")+
  geom_vline(xintercept = 600,col="red")+geom_vline(xintercept = 7000,col="red")
p


p<-meta_data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 650,col="red") +geom_vline(xintercept = 50000,col="red")+ 
  geom_hline(yintercept = 600,col="red") +geom_hline(yintercept = 7000,col="red") +
  facet_wrap(~sample)
p

p<-meta_data %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +ylab("log10 Cell density")+
  geom_vline(xintercept = 0.8,col="red")
p



p<-meta_data %>% 
  ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 5,col="red")

p



###Cell-level filtering
Merge <-subset(Merge, subset= nFeature_RNA > 600 & percent.mt < 5 & nFeature_RNA < 7000 & nCount_RNA < 50000 &Epcam>0 & Spp1>0)






###table for after QC
meta_data<-filtered_data@meta.data
samples<-unique(meta_data$sample)
samples

results<-NULL

i=6
for (i in 1:length(samples)){
  s=samples[i]
  nCount_RNA<-summary(meta_data$nCount_RNA[meta_data$sample==s])
  
  nFeature_RNA<-summary(meta_data$nFeature_RNA[meta_data$sample==s])
  nFeature_RNA
  
  
  cell_count<-table(meta_data$sample)[s]
  result<-NULL
  result$SampleID<-s
  result$Median.Genes.per.Cell<-round(nFeature_RNA[3],digits = 0)
  result$Median.UMI.Counts.per.cell<-as.integer(nCount_RNA[3])
  result$Number.of.Cells<-round(cell_count,digits = 0)
  results<-rbind(result,results)
  
}

results


