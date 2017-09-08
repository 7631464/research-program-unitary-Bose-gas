#!/bin/sh

#PBS -N BoseGS
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00

module load devel
cd $PBS_O_WORKDIR

./beccalc < parameter7.inp >output.dat
