#!/bin/bash

set -e

export CXX=clang++-3.7
export CC=clang-3.7
export FEELPP_DEP_INSTALL_PREFIX=/usr/local
export LD_LIBRARY_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib:${FEELPP_DEP_INSTALL_PREFIX}/lib/paraview-4.4:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH
export PYTHONPATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/python2.7/site-packages:${FEELPP_DEP_INSTALL_PREFIX}/lib/paraview-4.4/site-packages:$PYTHONPATH
export MANPATH=${FEELPP_DEP_INSTALL_PREFIX}/share/man:$MANPATH

echo '--- build directory'
if [ -d build ]; then rm -rf build; fi
mkdir -p build
cd build

echo '--- configure -r'
../configure -r  --cmakeflags="-DFEELPP_ENABLE_VTK_INSITU=ON -DCMAKE_INSTALL_PREFIX=${FEELPP_HOME}"

echo '--- make feelpp library'
make -j4 feelpp
