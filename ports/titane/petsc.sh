#!/bin/sh

# PETSC
# tested on Fri, Dec 14, 2012

# file petsc-lite-3.3-p5.tar.gz
# size 6741006
# md5  f7a240400dcfcc3eb877a047efb20d31
# sha1 ede11c0459ba37de1fb44218230689c61ac1eccd

export PREFIX=$CCCWORKDIR/local-titane
export MPI_ROOT=/applications/openmpi-1.4.2_gnu
export PETSC_DIR=$PWD

module unload intel/11.1.056
module unload bullxmpi/1.0.2
module load cmake/2.8.0
module load gcc/4.6.3
module load openmpi/1.4.2_gnu

./configure \
  --prefix=$PREFIX \
  --with-shared-libraries=1 \
  --with-fortran=0 \
  --with-mpi-include=$MPI_ROOT/include \
  --with-mpi-lib=[$MPI_ROOT/lib/libmpi_cxx.so,$MPI_ROOT/lib/libmpi.so] \
  --with-blas-lib=$PREFIX/lib/libblas.so \
  --with-lapack-lib=$PREFIX/lib/liblapack.so

make -j8 PETSC_DIR=$PWD PETSC_ARCH=linux all install
