#!/bin/bash
#PBS -q D 
#PBS -l select=1:ncpus=96:mpiprocs=48:ompthreads=2
#PBS -N swpc-N

module load PrgEnv-intel; module load cray-pals; module load cray-pmi
module load HDF5/1.14.5/intel/2024.2.1 NetCDF/4.9.2/intel/2024.2.1

export PALS_CPU_BIND=core

cd $PBS_O_WORKDIR
export DATASET=${HOME}/dataset

mpirun -d 2 ./bin/swpc_3d.x -i example/input_NJapan.inf
