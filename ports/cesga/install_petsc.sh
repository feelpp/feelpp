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

    ./configure  COPTFLAGS='-O3' FOPTFLAGS='-O3' \
        --with-shared-libraries=1 \
        --with-debugging=0 \
        --prefix=$PREFIX \
        --with-cc=mpicc \
        --with-cxx=mpicxx \
        --with-fc=mpifort \
        --with-cxx-dialect=C++11 \
        --CFLAGS="-mtune=haswell -O3 -mfma -malign-data=cacheline -finline-functions -fopenmp -lc" \
        --CXXFLAGS="-march=haswell -O3 -mfma -malign-data=cacheline -finline-functions -std=c11 -fopenmp -lc" \
        --FFLAGS="-march=haswell -O3 -mfma -malign-data=cacheline -finline-functions -fopenmp -lc" \
        --download-ml \
        --download-mumps \
        --download-pastix \
        --download-ptscotch \
        --download-suitesparse=1 \
        --with-blas-lapack-lib="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl" \
        --with-scalapack-lib="-L/opt/cesga/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl-gnu"


    make all
#    make test
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


