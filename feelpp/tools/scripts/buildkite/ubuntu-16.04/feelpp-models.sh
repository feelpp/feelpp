#!/bin/bash

set -e

#export CXX=clang++-3.8
#export CC=clang-3.8
export FEELPP_DEP_INSTALL_PREFIX=/usr/local
#export LD_LIBRARY_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib:${FEELPP_DEP_INSTALL_PREFIX}/lib/paraview-4.4:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH
export PYTHONPATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/python2.7/site-packages:/usr/lib/paraview/site-packages:$PYTHONPATH
export MANPATH=${FEELPP_DEP_INSTALL_PREFIX}/share/man:$MANPATH

cd build

echo '--- make feelpp model library'
make -j10 feelpp-models
