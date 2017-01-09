#!/bin/bash

set -e

#export CXX=clang++-3.8
#export CC=clang-3.8
# export FEELPP_DEP_INSTALL_PREFIX=/usr/local
# #export LD_LIBRARY_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib:${FEELPP_DEP_INSTALL_PREFIX}/lib/paraview-4.4:$LD_LIBRARY_PATH
# export PKG_CONFIG_PATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/pkgconfig:$PKG_CONFIG_PATH
# export PYTHONPATH=${FEELPP_DEP_INSTALL_PREFIX}/lib/python2.7/site-packages:/usr/lib/paraview/site-packages:$PYTHONPATH
# export MANPATH=${FEELPP_DEP_INSTALL_PREFIX}/share/man:$MANPATH


#echo '--- get docker'
#sudo apt-get update
#sudo apt-get install -y -q -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" docker-engine

echo '--- get feelpp/docker'
if [ ! -d docker]; then git clone --depth=1 https://github.com/feelpp/docker; fi

echo '--- building feelpp-libs'
cd docker/dockerfiles/feelpp-libs && bash mkimg.sh -f ${TARGET} --jobs=${JOBS} --branch=${BUILDKITE_BRANCH} --cxx="${CXX}" --cc="${CC}" --

