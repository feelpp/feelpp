#! /bin/bash


# DO NOT FORGET to modify the CMakeLists.txt 
# - has to take into account the $SLEPC_DIR variable !
# - Exclude all petsc installation on system if custom install



cmake \
  -DCMAKE_CXX_COMPILER=`which g++` \
  -DCMAKE_C_COMPILER=`which gcc` \
  -DCMAKE_INSTALL_PREFIX=/data/software/install/gmsh-2.9.4 \
  -DCMAKE_BUILD_TYPE=release \
  -DENABLE_BUILD_LIB=ON \
  -DENABLE_BUILD_SHARED=ON \
  -DENABLE_BUILD_DYNAMIC=ON \
  -DENABLE_MPI=ON \
  -DENABLE_MUMPS=ON \
  -DENABLE_OPENMP=ON  \
  ..
