#! /bin/bash


# DO NOT FORGET to modify the CMakeLists.txt 
# - has to take into account the $SLEPC_DIR variable !
# - Exclude all petsc installation on system if custom install

VERSION=1.3.4
basedir=$HOME
#basedir=/tmp

set -e
# check to see if protobuf folder is empty
if [ ! -d "$basedir/software/install/cln" ]; then
    wget http://www.ginac.de/CLN/cln-$VERSION.tar.bz2
    tar xvjf cln-$VERSION.tar.bz2
    cd cln-$VERSION
    mkdir build
    cd build 
    ../configure --prefix=$basedir/software/install/cln --disable-dependency-tracking
    make -j$NPROCS install
    ls -R $basedir/software/install/cln
else
    echo "Using cached directory $basedir/software/install/cln";
fi
