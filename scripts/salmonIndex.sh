#!/bin/bash -login
#PBS -l walltime=04:00:00,nodes=1:ppn=1,mem=64Gb
#mdiag -A ged
#PBS -m abe
#PBS -N salmonIndex

cd $PBS_O_WORKDIR
module load GNU/4.8.2
# The following have been reloaded with a version change:
# 1) Boost/1.47.0 => Boost/1.55.0     2) GNU/4.4.5 => GNU/4.8.2     3) OpenMPI/1.4.3 => OpenMPI/1.6.5     4) R/2.15.1 => R/3.1.0
module load salmon/0.5.0
# The following have been reloaded with a version change:
# 1) CMake/2.8.5 => CMake/3.1.0
salmon index --index "$index" --transcripts "$transcriptome" --type quasi

