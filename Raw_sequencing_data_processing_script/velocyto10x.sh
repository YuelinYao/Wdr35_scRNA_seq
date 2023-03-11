########

#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N velocyto10x              
#$ -cwd                  
#$ -l h_rt=48:00:00 
#$ -l h_vmem=16G
#$ -pe sharedmem 12
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load Python
module add anaconda
source activate Yuelin

# Run the program
ncores=30

mask=/exports/igmm/eddie/boulter-lab/sw_singlecell/velocyto/RNA_velocity/mm10_rmsk.gtf
#from UCSC genome browser 
ref=/exports/eddie/scratch/s1914230/refdata-cellranger-mm10-1.2.0/genes/genes.gtf
# from 10x genomics

for i in $(ls /exports/eddie/scratch/s1914230/Single_Cell/result/)

do
velocyto run10x -m $mask /exports/eddie/scratch/s1914230/Single_Cell/result/$i  $ref

done
