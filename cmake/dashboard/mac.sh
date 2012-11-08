#!/bin/sh
# $1 provides the path to the toplevel source Feel++ directory
# e.g. $HOME/Devel/FEEL/feelpp.git
# $2 provide the ctest mode: e.g. Nightly or Experimental

export PETSC_DIR=/opt/local/lib/petsc
export SLEPC_DIR=/opt/local/lib/petsc
export PATH=/opt/local/bin/:/opt/local/lib/openmpi/bin:$PATH

COMMON="/opt/local/bin/ctest -VV -S $1/cmake/dashboard/testsuite.cmake,FEELPP_CTEST_CONFIG=$1/cmake/dashboard/feelpp.site.`hostname -s`.cmake,FEELPP_MODE=$2"

if [ -x /opt/local/bin/g++-mp-4.6 ]; then
    $COMMON,FEELPP_CXXNAME=gcc-4.6.2,FEELPP_CXX=/opt/local/bin/g++-mp-4.6,FEELPP_EXPLICIT_VECTORIZATION=SSE2
fi
if [ -x /opt/local/bin/g++-mp-4.7 ]; then
    $COMMON,FEELPP_CXXNAME=gcc-4.7.2,FEELPP_CXX=/opt/local/bin/g++-mp-4.7,FEELPP_EXPLICIT_VECTORIZATION=SSE2
fi

if [ -x /opt/local/bin/clang-mp-3.1 ]; then
    $COMMON,FEELPP_CXXNAME=clang-3.1,FEELPP_CXX=/opt/local/bin/clang-mp-3.1,FEELPP_EXPLICIT_VECTORIZATION=SSE2
fi
