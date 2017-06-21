#! /bin/sh

BUILDDIR=$PWD
SCRIPTPATH=$( cd $(dirname $0) ; pwd -P )
cd $BUILDDIR

. $SCRIPTPATH/environment
. $SCRIPTPATH/modules

PROF=gcc630
FEELPP_SOURCES=$LUSTRE/devel/feelpp.git
export CMAKE_BUILD_TYPE=RELEASE
export CLANG_DIR=/opt/cesga/easybuild/software/Clang/3.9.1-foss-2017a
export GMSH_DIR=${INSDIR}/$PROF/gmsh/2.16.0
export PETSC_DIR=${INSDIR}/$PROF/petsc/3.7.3
export CLN_DIR=${INSDIR}/$PROF/cln/1.3.4
export PETSC_ARCH=
export CC=${CLANG_DIR}/bin/clang
export CXX=${CLANG_DIR}/bin/clang++
#export CC=`which gcc` \
    #export CXX=`which g++` \

${INSDIR}/$PROF/cmake/3.8.0/bin/cmake \
    -DFEELPP_PETSC_ENABLE_TESTS=OFF \
    $FEELPP_SOURCES
