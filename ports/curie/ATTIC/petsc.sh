#! /bin/bash

. ~/.bash_profile
module load c++/gnu/4.5.1   
prefix=$WORKDIR/local-gcc45
export PETSC_DIR=`pwd` 
./configure --prefix=$prefix \
    --with-shared-libraries=1 \
    --with-fortran=0 \
    --with-cc=/usr/local/gcc-4.5.1/bin/gcc \
    --with-cxx=/usr/local/gcc-4.5.1/bin/g++ \
    --with-mpi-include=$MPI_ROOT/include \
    --with-mpi-lib=[$MPI_ROOT/lib/libmpi_cxx.so,$MPI_ROOT/lib/libmpi.so] \
    --with-blas-lib=$prefix/lib/libblas.so \
    --with-lapack-lib=$prefix/lib/liblapack.so

