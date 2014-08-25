#! /bin/bash

. ~/.bash_profile
#module load fermi/bgclang.profile
module load fermi/bgqgcc472.feelpp

cmake \
 -DCMAKE_CXX_FLAGS="-Bdynamic -fPIC -std=c++11 -std=c++0x"\
 -DCMAKE_EXE_LINKER_FLAGS="-dynamic -std=c++11 -std=c++0x"\
 -DCMAKE_INSTALL_PREFIX=$WORK/local/gmsh/2.8.4-bgq \
 -DENABLE_PETSC=OFF \
 -DENABLE_KBIPACK=OFF\
 -DENABLE_METIS=OFF\
 -DENABLE_GMM=OFF\
 -DENABLE_BAMG=OFF\
 -DENABLE_BUILD_SHARED=ON\
 ..
 #-DGMP_INC=$GMP_INC\
 #-DGMP_LIB=$GMP_LIB\
 #..

#make all lib shared
