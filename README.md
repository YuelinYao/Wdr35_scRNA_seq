# ScRNA_seq_K19WDR
Codes and scripts for scRNA-seq analysis. This project includes the following steps:


## 1. Processing the raw sequencing data using Cell ranger pipeline (version 5.0.0). 
(1) FASTQ files was generated using ‘cellranger mkfastq’ from raw base call files. This step was performed by the sequencing company and they provided a summary report: SummaryReport.pdf

(2) We then aligned FASTQ files to a pre-build mouse reference genome (GRCm38/mm10) with cellranger count to generate single cell feature counts. The reference genome was downloaded from 10x genomes website (https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). 


(3) We also additionally calculated the unspliced and spliced count matrices for RNA velocity analysis using velocyto (version 0.17). 


## 2. scRNA-seq data analysis, mainly with Seurat (version 4.0.6). 

(1) We first created seurat objest, merging all the samples. 

(2) We then perform QC:
  
  Gene level: only keep genes which are expressed in 3 or more cells in each sample.
  
  Cell level: nFeature_RNA > 600, percent.mt < 5, nFeature_RNA < 7000, nCount_RNA < 50000, Epcam>0, Spp1>0, and remove doublets used DoubletFinder (V2.0).

(3) Preprocessing steps: normalisation, finding variable features, scale, PCA, UMAP, TSNE. 

(4) Clustering: we applied a **selective inference** framework from lucy gao (https://www.lucylgao.com/clusterpval/index.html) to identify "true" 
clusters. It tests for the difference in means after any type of clustering and is especially efficient when applied to hierarchical clustering. 

(5) Differential expression analysis between cluster 1 and cluster 2. For functional enrichment analysis, We obtained the background genes selected from a previous public dataset (https://www.nature.com/articles/s41556-019-0402-6). The background genes is a set of genes known to expressed in the EpCam sorted cholangiocytes from mouse.

(6) Perform GO enrichment analysis for the identified clusters. 

(7) Receptor-ligand analysis was performed using R package SingleCellSingalR. Ref: https://www.bioconductor.org/packages/release/bioc/vignettes/SingleCellSignalR/inst/doc/UsersGuide.html

(8) RNA velocity analysis. Ref: https://scvelo.readthedocs.io/getting_started/
