#!/bin/bash
#set -x


# $1 provides the path to the toplevel source Feel++ directory
# e.g. $HOME/Devel/FEEL/feelpp.git
COMMON="ctest -VV -S $1/cmake/dashboard/testsuite.cmake,FEELPP_CTEST_CONFIG=$1/cmake/dashboard/feelpp.site.`hostname -s`.cmake,FEELPP_MODE=$2"

#To make available specific compilation instead of "whatever is available"
compiler_list=${3:-"gcc47,gcc48,gcc49,clang"}
#do_gcc46=`echo $compiler_list | grep gcc46`
do_gcc47=`echo $compiler_list | grep gcc47`
do_gcc48=`echo $compiler_list | grep gcc48`
do_gcc49=`echo $compiler_list | grep gcc49`
do_clang=`echo $compiler_list | grep clang`

export FEELPP_SCRATCHDIR=/tmp/feel-logs/
rm -rf $FEELPP_SCRATCHDIR
#if [ ! -z "$do_gcc46" -a -x /usr/bin/g++-4.6 ]; then
#    export FEELPP_WORKDIR=/tmp/feel-gcc46
#    rm -rf $FEELPP_WORKDIR 
#    $COMMON,FEELPP_CXXNAME=gcc-4.6,FEELPP_CXX=/usr/bin/g++-4.6 
#    rm -rf $FEELPP_WORKDIR 
#fi
if [ ! -z "$do_gcc47" -a -x /usr/bin/g++-4.7 ]; then
    export FEELPP_WORKDIR=/tmp/feel-gcc47
    rm -rf $FEELPP_WORKDIR 
    $COMMON,FEELPP_CXXNAME=gcc-4.7,FEELPP_CXX=/usr/bin/g++-4.7
    rm -rf $FEELPP_WORKDIR 
fi
if [ ! -z "$do_gcc48" -a -x /usr/bin/g++-4.8 ]; then
    export FEELPP_WORKDIR=/tmp/feel-gcc48
    rm -rf $FEELPP_WORKDIR 
    $COMMON,FEELPP_CXXNAME=gcc-4.8,FEELPP_CXX=/usr/bin/g++-4.8 
    rm -rf $FEELPP_WORKDIR 
fi
if [ ! -z "$do_gcc49" -a -x /usr/bin/g++-4.9 ]; then
    export FEELPP_WORKDIR=/tmp/feel-gcc49
    rm -rf $FEELPP_WORKDIR 
    $COMMON,FEELPP_CXXNAME=gcc-4.9,FEELPP_CXX=/usr/bin/g++-4.9 
    rm -rf $FEELPP_WORKDIR 
fi
if [ ! -z "$do_clang" -a -x /usr/bin/clang++-3.5 ]; then
    export FEELPP_WORKDIR=/tmp/feel-clang
    rm -rf $FEELPP_WORKDIR 
    clang_version=`echo | clang -dM -E - | grep clang_version | awk '{print $3}' | sed "s/\"//g"`
    $COMMON,FEELPP_CXXNAME=clang-$clang_version,FEELPP_CXX=/usr/bin/clang++-3.5
    rm -rf $FEELPP_WORKDIR 
fi
rm -rf $FEELPP_SCRATCHDIR
