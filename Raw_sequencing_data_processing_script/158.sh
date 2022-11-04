####
#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N sc158         
#$ -cwd                  
#$ -l h_rt=48:00:00 
#$ -l h_vmem=16G
#$ -pe sharedmem 12
#$ -M s1914230@ed.ac.uk
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem
# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load Python
module add anaconda
module load igmm/apps/cellranger/5.0.0
# Run the program


tx=/exports/eddie/scratch/s1914230/refdata-gex-mm10-2020-A
run1=/exports/eddie/scratch/s1914230/Single_Cell/SW_NS2K_190321/outs/fastq_path
run2=/exports/eddie/scratch/s1914230/Single_Cell/SW_260321/outs/fastq_path

cd /exports/eddie/scratch/s1914230/Single_Cell/result

s=K19Wdr3512m_158


cellranger count --id=$s \
--transcriptome=$tx \
--fastqs=$run1,$run2 \
--sample=$s
