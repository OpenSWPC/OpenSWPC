#!/bin/bash

#PBS -q regular-g
#PBS -l select=4:mpiprocs=1
#PBS -l walltime=01:00:00
#PBS -W group_list=${GROUP}  #<-- 自分の所属グループに変更した後このコメントを削除
#PBS -j oe
#PBJ -N miyabi_job
#PBS -o miyabi_job.out

cd ${PBS_O_WORKDIR}
module load nvidia nv-hpcx netcdf hdf5 netcdf-fortran
mpirun ./bin/swpc_3d.x -i example/input.inf