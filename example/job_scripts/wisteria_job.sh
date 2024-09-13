#!/bin/bash

#PJM -L    rscgrp=regular-o
#PJM -L    node=2x2:mesh
#PJM --mpi proc=4
#PJM -L    elapse=00:30:00
#PJM -g    ${GROUP}  #<-- 自分の所属グループに変更した後コメントを削除
#PJM --omp thread=48
#PJM -N    bdec-001         
#PJM -o    bdec-001.out
#PJM -j 
# ---------- 

# 計算に必要なモジュールのロード．OpenSWPCはNetCDFを利用するため，それと関連のモジュールをロードしている
module load fj fjmpi netcdf hdf5 netcdf-fortran

# プログラムの実行．mpiコードの実行コマンドは mpiexec
mpiexec ./bin/swpc_3d.x -i in/input.inf
