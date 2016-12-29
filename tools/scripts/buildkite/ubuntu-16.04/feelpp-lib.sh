#!/bin/bash

set -e

#export CXX=clang++-3.8
#export CC=clang-3.8
# export FEELPP_DEP_INSTALL_PREFIX=/usr/local
# #export LD_LIBRARY_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib:${FEELPP_DEP_INSTALL_PREFIX}/lib/paraview-4.4:$LD_LIBRARY_PATH
# export PKG_CONFIG_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH
# export PYTHONPATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/python2.7/site-packages:/usr/lib/paraview/site-packages:$PYTHONPATH
# export MANPATH=${FEELPP_DEP_INSTALL_PREFIX}/share/man:$MANPATH

echo '--- build directory'
if [ -d build ]; then rm -rf build; fi
mkdir -p build
cd build
echo 'FEELPP_HOME: $FEELPP_HOME'
echo '--- configure -r'
#../configure -r  --cmakeflags="-DFEELPP_ENABLE_VTK_INSITU=ON -DCMAKE_INSTALL_PREFIX=${FEELPP_HOME}"
../configure -r  --cmakeflags="-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}"

echo '--- make feelpp library'
make -j10 feelpp
