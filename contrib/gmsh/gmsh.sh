#! /bin/bash

wget http://geuz.org/gmsh/src/gmsh-svn-source.tgz

export GMSH_DIR=$1

cmake .. -DCMAKE_INSTALL_PREFIX=$1 -DCMAKE_BUILD_TYPE=release -DENABLE_MPI=ON -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 -DENABLE_SLEPC=0

make -j${2} lib
make -j${2} shared
make -j${2} install

