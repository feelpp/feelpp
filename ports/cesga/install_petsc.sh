#!/bin/sh

. ./environment

################################################################################
NAME=petsc
VERSION=3.7.3
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

    wget -c http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${VERSION}.tar.gz

    tar -xvzf petsc-lite-${VERSION}.tar.gz
    cd petsc-${VERSION}

    ./configure --with-shared-libraries=1 \
        --with-debugging=0 \
        COPTFLAGS='-O3' FOPTFLAGS='-O3' \
        --prefix=$PREFIX \
        --download-suitesparse=1 \
        --download-ml \
        --download-blacs \
        --download-scalapack \
        --download-fblaslapack \
        --download-mumps \
        --download-pastix \
        --download-ptscotch
    #--with-mpi-dir=$MPI_ROOT \
    #--download-metis \
        #--download-parmetis \
    make all
    #make test
    make install

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
