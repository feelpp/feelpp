#! /bin/bash

#. ~/.bash_profile
#module load fermi/bgclang.profile
#module load fermi/bgqgcc472.feelpp

cmake \
 -DCMAKE_CXX_FLAGS="-Bdynamic -fPIC  -std=c++0x"\
 -DCMAKE_EXE_LINKER_FLAGS="-dynamic  -std=c++0x"\
 -DCMAKE_INSTALL_PREFIX=$WORKDIR/local/gmsh/2.8.5-bgq \
 -DENABLE_PETSC=OFF \
 -DENABLE_KBIPACK=OFF\
 -DENABLE_METIS=OFF\
 -DENABLE_GMM=OFF\
 -DENABLE_BAMG=OFF\
 -DENABLE_BUILD_SHARED=OFF\
 -DENABLE_BUILD_LIB=ON \
 ..
 #-DGMP_INC=$GMP_INC\
 #-DGMP_LIB=$GMP_LIB\
 #..

#make all lib shared
