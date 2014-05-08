#! /bin/bash

set -x
. ~/.bash_profile
#module load bgq-gnu/4.4.6
#module load compilers/bgclang.feelpp
module load fermi/bgqgcc472.feelpp

export CC=$WORK/local/bin/mpicc
export CXX=$WORK/local/bin/mpic++
prefix=$WORK/local/petsc/3.4.4-bgq
export MPI_ROOT=/bgsys/drivers/ppcfloor/comm/
export PETSC_DIR=`pwd`
./configure  --prefix=$prefix\
    --with-shared-libraries=0 \
    --with-debugging=0 \
    --download-f-blas-lapack=1 \
    --with-mpi=1 \
    --with-mpi-include=$MPI_ROOT/include \
    --with-mpi-compilers=1 \
    --with-mpi-lib=[$MPI_ROOT/lib/libmpichcxx-gcc.a,$MPI_ROOT/lib/libmpich-gcc.a] \
    --download-mumps=0 \
    --download-pastix=0 \
    --download-ptscotch=1 \
    --download-scalapack=1 \
    --download-ml=1 \
    --download-superlu=1 \
    --download-hypre=1

make PETSC_DIR=/gpfs/work/LI03s_DDFD/sources/petsc/petsc-3.4.4 PETSC_ARCH=arch-linux2-c-opt all
make PETSC_DIR=/gpfs/work/LI03s_DDFD/sources/petsc/petsc-3.4.4 PETSC_ARCH=arch-linux2-c-opt install



#    --with-blas-lib=[/cineca/prod/libraries/blas/2007/bgq-xl--1.0/lib/libblas.a] \
#    --with-lapack-lib=[/cineca/prod/libraries/lapack/3.4.1/bgq-xl--1.0/lib/liblapack.a ] \
