rm(list=ls())

### this script is to extract background genes for go/kegg analysis
### Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123133
setwd("~/result/")
vitro<-read.table("GSE123133_TPM_RNAseq_vitro.txt",sep="\t",header = T)
vitro<-vitro[!duplicated(vitro$X),]
rownames(vitro)<-vitro$X

vitro<-vitro[,-1]
vitro_bulk_bkg<-rowSums(vitro>0)
vitro_tpm<-vitro[!vitro_bulk_bkg==0,]



vivo<-read.table("GSE123133_TPM_RNAseq_vivo.txt",sep="\t",header = T)

vivo<-vivo[!duplicated(vivo$ext_gene),]
rownames(vivo)<-vivo$ext_gene
vivo<-vivo[,-1]
vivo_bulk_bkg <- rowSums(vivo>0)
vivo_tpm<-vivo[!vivo_bulk_bkg==0,]

background_gene<-intersect(rownames(vivo_tpm),rownames(vitro_tpm))

#save(background_gene,file="background_gene.Rdata")
