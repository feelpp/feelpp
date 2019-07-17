#! /bin/sh

. ./environment

################################################################################
NAME=gmsh
VERSION=2.16.0
CMAKE_BUILD_TYPE=RELEASE
# Files are installed in PREFIX = $INSDIR/$prof (via script args)
################################################################################

export $CMAKE_BUILD_TYPE

#-------------------------------------------------------------------------------
# OPTIONS (same for all installs)
#-------------------------------------------------------------------------------

HELP="\
Usage: $0 [opt] profilename

Install ${0:2:${#0}} software

Options:
    -r  --rm    remove files
    -h  --help  print help

"

opts=${@:1:$#-1}
for prof; do true; done

if (( $# < 1 )); then
    printf "$HELP"
    exit 0
fi
if [ "${prof:0:1}" == "-" ]; then
    printf "$HELP"
    exit 0
fi

PROFDIR=$INSDIR/$prof
PREFIX=$PROFDIR/$NAME/$VERSION
:
for key in ${opts}
do
    case $key in
        -r|--rm)
            printf "Remove file in $SRCDIR\n"
            rm -rvfI $SRCDIR
            printf "Remove file in $PREFIX\n"
            rm -rvfI $PREFIX
            exit 0
            ;;
        *)
            # unknown option
            ;;
    esac
done

#-------------------------------------------------------------------------------
# INSTALLATION
#-------------------------------------------------------------------------------

set -e

# Avoid overding files.
if [ ! -d "$PREFIX" ]; then
    . ./modules
    printf "installation directory: ${PREFIX}\n"
    mkdir -p $ARCDIR
    mkdir -p $PREFIX
    mkdir -p $SRCDIR
    cd $SRCDIR



    wget http://gmsh.info/src/gmsh-${VERSION}-source.tgz
    tar -xvzf gmsh-$VERSION-source.tgz
    cd gmsh-$VERSION-source
    mkdir build
    cd build
    cmake \
        -DCMAKE_CXX_COMPILER=`which g++` \
        -DCMAKE_C_COMPILER=`which gcc` \
        -DCMAKE_INSTALL_PREFIX=$PREFIX \
        -DCMAKE_BUILD_TYPE=release \
        -DENABLE_BUILD_LIB=ON \
        -DENABLE_BUILD_SHARED=ON \
        -DENABLE_BUILD_DYNAMIC=ON \
        -DENABLE_MPI=ON \
        -DENABLE_MUMPS=OFF \
        -DENABLE_OPENMP=ON  \
        ..
    make -j$NPROCS
    make install

    export GMSH_DIR=$PREFIX


    # Create an archive.
    cd $PROFDIR
    ARCFILE=$NAME-$VERSION-$prof.tar
    tar -cvf $ARCFILE $SRCDIR/*
    mv $ARCFILE $ARCDIR
    rm -rvfI $SRCDIR
else
    printf "Install or cache file already exist for this profile. You should
    remove them before going further.\n$SRCDIR\n$PREFIX\n";
fi
