#! /bin/bash

set -x
. ~/.bash_profile
#module load bgq-gnu/4.4.6
module load compilers/bgclang.feelpp
#module load fermi/bgqgcc472.feelpp

export CC=$WORK/bgclang/mpi/bgclang/bin/mpicc
export CXX=$WORK/bgclang/mpi/bgclang/bin/mpic++11
export FC=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-gfortran
#export FC=$WORK/local/bin/mpif90
#export CXX=$WORK/local/bin/mpic++
#export CC=/bgsys/drivers/V1R2M1/ppc64/comm/bin/gcc/mpicc
#export FC=/bgsys/drivers/V1R2M1/ppc64/comm/bin/gcc/mpif90

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
    --with-mpi-lib=[/bgsys/drivers/V1R2M1/ppc64/comm/lib/libmpich-gcc.a,/bgsys/drivers/V1R2M1/ppc64/comm/lib/libopa-gcc.a,/bgsys/drivers/V1R2M1/ppc64/comm/lib/libmpl-gcc.a,/bgsys/drivers/V1R2M1/ppc64/comm/lib/libpami-gcc.ay,/bgsys/drivers/V1R2M1/ppc64/spi/lib/libSPI.a,/bgsys/drivers/V1R2M1/ppc64/spi/lib/libSPI_cnk.a] \
    --download-mumps=0 \
    --download-pastix=0 \
    --download-ptscotch=1 \
    --download-scalapack=1 \
    --download-ml=0 \
    --download-superlu=1 --CXX=$CXX --CC=$CC --FC=$FC

#    --download-hypre=1

make PETSC_DIR=/gpfs/work/LI03s_DDFD/sources/petsc/petsc-3.4.4 PETSC_ARCH=arch-linux2-c-opt all
make PETSC_DIR=/gpfs/work/LI03s_DDFD/sources/petsc/petsc-3.4.4 PETSC_ARCH=arch-linux2-c-opt install



#    --with-blas-lib=[/cineca/prod/libraries/blas/2007/bgq-xl--1.0/lib/libblas.a] \
#    --with-lapack-lib=[/cineca/prod/libraries/lapack/3.4.1/bgq-xl--1.0/lib/liblapack.a ] \
