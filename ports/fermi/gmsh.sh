#! /bin/bash

. ~/.bash_profile
#module load fermi/bgclang.profile
module load fermi/bgqgcc472.feelpp

cmake \
 -DCMAKE_INSTALL_PREFIX=$WORK/local/gmsh/2.8.4-bgq \
 -DENABLE_PETSC=OFF \
 -DENABLE_BUILD_LIB=ON\
 -DGMP_INC=$GMP_INC\
 -DGMP_LIB=$GMP_LIB\
 ..

#make all lib shared
