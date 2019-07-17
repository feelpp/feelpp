#!/bin/sh

. ./environment

################################################################################
NAME=boost
VERSION=1.62.0
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



    VERSION2=`echo "$VERSION" | tr . _`
    NPROCS=${2:-4}
    wget http://sourceforge.net/projects/boost/files/boost/${VERSION}/boost_${VERSION2}.tar.bz2/download -O boost_${VERSION}.tar.bz2;

    tar -xvf boost_${VERSION}.tar.bz2

    cd boost_${VERSION2}

    echo "using mpi ;" >> user-config.jam
    echo "" >> user-config.jam
    ./bootstrap.sh
    ./bjam cxxflags="-std=c++14"  -j$NPROCS install \
        --layout=tagged \
        --prefix=$PREFIX \
        --user-config=user-config.jam \
        variant=release \
        threading=multi \
        link=shared




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
