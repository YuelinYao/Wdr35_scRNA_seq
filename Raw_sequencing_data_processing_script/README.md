# Raw sequencing_data processing

Data was processed on Eddie (High-performance computing, University of Edinburgh, http://www.ecdf.ed.ac.uk/).

SW_NS2K_190321, SW_NS2K_190321 are the output from cellranger mkfastq, provdied by the sequencing company (See SummaryReport.pdf).

We have 7 samples: 

  **Wild-type group**: 155, 158, 176, 200
  
  **Mutated group**: 168, 177, 186
  
See Yuelin.yml for conda enviroment

## Download reference genome
```
########1. Download the reference genome####
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -xzvf refdata-gex-mm10-2020-A.tar.gz

#ref: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input
#ref: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/glossary#gem_well
```

## Cell ranger count
```
###################################################
###################################################
############Cell ranger count######################
###################################################
qsub 155.sh
qsub 158.sh
qsub 168.sh
qsub 176.sh
qsub 177.sh
qsub 186.sh
qsub 200.sh
```

## Velocyto

```
qsub velocyto10x.sh
```
