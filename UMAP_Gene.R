load("combined_filter.Rdata") # load seruat object


genelist<-c("Smad1","Smad2","Smad3","Smad4","Smad5","Smad6","Smad7") # list of genes

# loop for plot umap for all the genes in genelist
for (g in 1:length(genelist)){
  if(genelist[g]%in%rownames(combined_filter)){
    pdf(file = paste0("~/SMAD/",genelist[g],".pdf"),width=4,height=3)
    a<-FeaturePlot(combined_filter, features =genelist[g] ,reduction = 'umap')
    print(a)
    dev.off()} else{
      print(genelist[g])
    } 
  
}
