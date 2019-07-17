#!/bin/sh

set +x
VERSION=1_62_0
VERSION2=`echo "$VERSION" | tr _ .`
basedir=${1:-$HOME}
NPROCS=${2:-4}
#basedir=/tmp

set -e
# check to see if protobuf folder is empty
if [ ! -d "$basedir/software/install/boost-${VERSION}" ]; then

#    if [ ! -f "boost_${VERSION}.tar.bz2"]; then 
        wget http://sourceforge.net/projects/boost/files/boost/${VERSION2}/boost_${VERSION}.tar.bz2/download -O boost_${VERSION}.tar.bz2;
#    fi

    tar xjf boost_${VERSION}.tar.bz2
    cd boost_${VERSION}

    echo "using mpi ;" >> user-config.jam
    echo "" >> user-config.jam
    ./bootstrap.sh
    ./bjam cxxflags="-std=c++14"  -j$NPROCS install \
           --layout=tagged \
           --prefix=$basedir/software/install/boost-${VERSION} \
           --user-config=user-config.jam \
           variant=release \
           threading=multi \
           link=shared
else
    echo 'Using cached directory $basedir/software/install/boost-${VERSION}';
    echo "Cached version: ${VERSION}"
fi
