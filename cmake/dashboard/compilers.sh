#!/bin/bash
ARCH=`uname -m`
SITE=`hostname`
#OS_VERSION=ubuntu-natty
OS_VERSION=`uname -s`
#WORK_DIR=/home/prudhomm/sources/
WORK_DIR=$HOME/Devel/FEEL
MAKE_ARGS="-j6"
#mkdir -p $WORK_DIR
# get the last version of the script
#svn export svn://scm.forge.imag.fr/var/lib/gforge/chroot/scmrepos/svn/life/trunk/life/trunk/cmake/dashboard/testsuite.cmake $WORK_DIR/testsuite.cmake
COMMON="ctest -VV -S $WORK_DIR/feel/cmake/dashboard/testsuite.cmake,FEEL_WORK_DIR=$WORK_DIR,FEEL_MAKE_ARGS=$MAKE_ARGS,FEEL_SITE=$SITE,FEEL_MODE=$1,FEEL_BUILD_STRING=$OS_VERSION-$ARCH"
#$COMMON-gcc-4.4.6,FEEL_CXX=g++-4.4,FEEL_EXPLICIT_VECTORIZATION=novec
#$COMMON-gcc-4.4.6,FEEL_CXX=g++-4.4,FEEL_EXPLICIT_VECTORIZATION=SSE2
#$COMMON-gcc-4.5.2,FEEL_CXX=g++-4.5,FEEL_EXPLICIT_VECTORIZATION=novec
$COMMON-gcc-4.5.2,FEEL_CXX=g++-4.5,FEEL_EXPLICIT_VECTORIZATION=SSE2
