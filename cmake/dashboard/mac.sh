#!/bin/sh
# $1 provides the path to the toplevel source Feel++ directory
# e.g. $HOME/Devel/FEEL/feelpp.git
# $2 provide the ctest mode: e.g. Nightly or Experimental

COMMON="/usr/local/bin/ctest -VV -S $1/cmake/dashboard/testsuite.cmake,FEELPP_CTEST_CONFIG=$1/cmake/dashboard/feelpp.site.`hostname -s`.cmake,FEELPP_MODE=$2"
 
if [ -x /usr/bin/clang ]; then
    $COMMON,FEELPP_CXXNAME=clang-3.4,FEELPP_CXX=/usr/bin/clang,FEELPP_EXPLICIT_VECTORIZATION=SSE2
fi
