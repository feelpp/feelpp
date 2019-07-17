#! /bin/sh

. ./environment

################################################################################
PROF=gcc630
NAME=feelpp
VERSION=0.102.0
NPROCS=4
CMAKE_BUILD_TYPE=RELEASE
FEELPP_SOURCES=$LUSTRE/devel/feelpp.git
# Files are installed in PREFIX = $INSDIR/$prof (via script args)
################################################################################

export $CMAKE_BUILD_TYPE
export CLANG_DIR=/opt/cesga/easybuild/software/Clang/3.9.1-foss-2017a
export GMSH_DIR=${INSDIR}/$PROF/gmsh/2.16.0
export PETSC_DIR=${INSDIR}/$PROF/petsc/3.7.3
export CLN_DIR=${INSDIR}/$PROF/cln/1.3.4
export PETSC_ARCH=
export CC=${CLANG_DIR}/bin/clang
export CXX=${CLANG_DIR}/bin/clang++
#export CC=`which gcc` \
    #export CXX=`which g++` \

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



    ${INSDIR}/$PROF/cmake/3.8.0/bin/cmake \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DFEELPP_PETSC_ENABLE_TESTS=OFF \
    $FEELPP_SOURCES

    make install-feelpp-lib -j${NPROCS}




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
