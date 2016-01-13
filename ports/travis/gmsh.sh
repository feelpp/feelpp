#! /bin/bash


# DO NOT FORGET to modify the CMakeLists.txt 
# - has to take into account the $SLEPC_DIR variable !
# - Exclude all petsc installation on system if custom install

VERSION=2.10.1
basedir=$HOME
#basedir=/tmp

set -e
# check to see if protobuf folder is empty
if [ ! -d "$basedir/software/install/gmsh" ]; then
    wget http://gmsh.info/src/gmsh-${VERSION}-source.tgz
    tar xvzf gmsh-$VERSION-source.tgz
    cd gmsh-$VERSION-source
    mkdir build
    cd build 
    cmake \
        -DCMAKE_CXX_COMPILER=`which g++` \
        -DCMAKE_C_COMPILER=`which gcc` \
        -DCMAKE_INSTALL_PREFIX=$basedir/software/install/gmsh \
        -DCMAKE_BUILD_TYPE=release \
        -DENABLE_BUILD_LIB=ON \
        -DENABLE_BUILD_SHARED=ON \
        -DENABLE_BUILD_DYNAMIC=ON \
        -DENABLE_MPI=OFF \
        -DENABLE_MUMPS=OFF \
        -DENABLE_OPENMP=ON  \
        ..
    make -j$NPROCS
    make install
    export GMSH_DIR=$basedir/software/install/gmsh
else
    echo 'Using cached directory $basedir/software/install/gmsh';
    echo "Cached version: ${VERSION}";
    export GMSH_DIR=$basedir/software/install/gmsh
fi
