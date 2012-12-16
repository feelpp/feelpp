#!/bin/bash

# GMSH
# tested on Fri, Dec 14, 2012

# file gmsh-2.6.1-source.tgz
# size 6193189
# md5  815511cad883db20b966ba0f773ab339
# sha1 6d470e2286ebf9534e1af8b407df59a86315e15e

export PREFIX=$CCCWORKDIR/local-titane
export MPI_ROOT=/applications/openmpi-1.4.2_gnu

module unload intel/11.1.056
module unload bullxmpi/1.0.2
module load cmake/2.8.0
module load gcc/4.6.3
module load openmpi/1.4.2_gnu

mkdir build
cd build

cmake \
  -DCMAKE_SHARED_LINKER_FLAGS="-L$PREFIX/lib;-L$MPI_ROOT/lib" \
  -DCMAKE_Fortran_COMPILER=/applications/gcc-4.6.3/bin/gfortran \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DENABLE_PETSC=OFF \
  ..

make -j8 all lib shared install 
