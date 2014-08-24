#! /bin/bash

set -x
. ~/.bash_profile

export CXX=powerpc64-bgq-linux-clang++
export CC=powerpc64-bgq-linux-clang

prefix=$WORKDIR/local/petsc/3.5.1-bgq
export MPI_ROOT=/bgsys/drivers/ppcfloor/comm/
export PETSC_DIR=`pwd`
./configure  --prefix=$prefix\
--known-level1-dcache-size=32768 --known-level1-dcache-linesize=32 --known-level1-dcache-assoc=0 --known-memcmp-ok=1 --known-sizeof-char=1 --known-sizeof-void-p=8 \
--known-sizeof-short=2 --known-sizeof-int=4 --known-sizeof-long=8 --known-sizeof-long-long=8 --known-sizeof-float=4 --known-sizeof-double=8 --known-sizeof-size_t=8 \
--known-bits-per-byte=8 --known-sizeof-MPI_Comm=4 --known-sizeof-MPI_Fint=4 --known-mpi-long-double=1 --known-mpi-c-double-complex=1\
    --with-shared-libraries=0 \
    --with-debugging=0 \
    --with-mpi=1 \
    --with-mpi-compilers=1 \
    --with-batch \
    --known-mpi-shared-libraries=0 \
    --with-blas-lapack-lib="[/bglocal/cn/pub/LAPACK/3.3.1/lib/liblapack.a,/bgsys/ibm_essl/prod/opt/ibmmath/essl/5.1/lib64/libesslbg.a]" \
    --CXX=mpicxx --CC=mpicc --FC=mpif90 --CXXFLAGS=-cxx=$CXX --CFLAGS=-cc=$CC


make PETSC_DIR=$WORKDIR/petsc-3.5.1 PETSC_ARCH=arch-linux2-c-opt all
#make PETSC_DIR=/gpfs/work/LI03s_DDFD/sources/petsc/petsc-3.4.4 PETSC_ARCH=arch-linux2-c-opt install



#    --with-blas-lib=[/cineca/prod/libraries/blas/2007/bgq-xl--1.0/lib/libblas.a] \
#    --with-lapack-lib=[/cineca/prod/libraries/lapack/3.4.1/bgq-xl--1.0/lib/liblapack.a ] \
