#!/bin/bash

wget -c http://geuz.org/gmsh/src/gmsh-svn-source.tgz
tar xzf gmsh-svn-source.tgz

cd `tar -tvf gmsh-svn-source.tgz | head -n 1 | sed "s/\// /g" | awk '{print $7}'`
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=${1} -DCMAKE_BUILD_TYPE=release -DENABLE_MPI=ON -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 -DENABLE_SLEPC=0

make -j${2} lib
make -j${2} shared
make -j${2} install

