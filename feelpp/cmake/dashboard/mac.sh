#!/bin/sh
# $1 provides the path to the toplevel source Feel++ directory
# e.g. $HOME/Devel/FEEL/feelpp.git
# $2 provide the ctest mode: e.g. Nightly or Experimental

export PATH=/usr/local/bin:$PATH
COMMON="/usr/local/bin/ctest -VV -S $1/cmake/dashboard/testsuite.cmake,FEELPP_CTEST_CONFIG=$1/cmake/dashboard/feelpp.site.`hostname -s`.cmake,FEELPP_MODE=$2"

if [ -x /usr/bin/clang++ ]; then
    export FEELPP_WORKDIR=/tmp/feel-clang
    rm -rf $FEELPP_WORKDIR
    clang_version=`echo | clang -dM -E - | grep clang_version | awk '{print $3}' | sed "s/\"//g"`
    $COMMON,FEELPP_CXXNAME=clang-$clang_version,FEELPP_CXX=/usr/bin/clang++
    rm -rf $FEELPP_WORKDIR
fi

# if [ -x /usr/local/bin/g++-4.8 ]; then
#     export FEELPP_WORKDIR=/tmp/feel-clang
#     rm -rf $FEELPP_WORKDIR
#     $COMMON,FEELPP_CXXNAME=gcc-4.8,FEELPP_CXX=/usr/local/bin/g++-4.8
#     rm -rf $FEELPP_WORKDIR
# fi
