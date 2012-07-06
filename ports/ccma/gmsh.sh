#!/bin/sh

module load cmake/2.8.7   
module load gcc/4.6.3  
module load openmpi/gcc/1.4.5  

export MPI=/softs/openmpi/gcc/1.4.5/
export prefix=/softs/cemracs/gmsh/2.6.0

cmake \
 -DCMAKE_SHARED_LINKER_FLAGS="-L$MPI/lib" \
 -DCMAKE_Fortran_COMPILER="$MPI/bin/mpif90" \
 -DCMAKE_INSTALL_PREFIX=$prefix \
 -DENABLE_MPI=ON \
 ..

#make all lib shared
