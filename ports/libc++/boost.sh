#!/bin/sh

set +x
VERSION=1_61_0
VERSION2=`echo "$VERSION" | tr _ .`
basedir=${1:-$HOME}
NPROCS=${2:-4}
#basedir=/tmp
export CXX=clang-libc++
INSTALL_DIR=$basedir/software/install/boost-${VERSION}
set -e
# check to see if protobuf folder is empty
if [ ! -d "$basedir/software/install/boost-${VERSION}" ]; then

#    if [ ! -f "boost_${VERSION}.tar.bz2"]; then 
        #wget http://sourceforge.net/projects/boost/files/boost/${VERSION2}/boost_${VERSION}.tar.bz2/download -O boost_${VERSION}.tar.bz2;
#    fi

    tar xjf boost_${VERSION}.tar.bz2
    cd boost_${VERSION}

    echo "using mpi ;" >> user-config.jam
    echo "" >> user-config.jam
    ./bootstrap.sh --with-toolset=clang --prefix=$INSTALL_DIR
    ./b2 toolset=clang cxxflags="-std=c++11 -stdlib=libc++ -I/usr/include/libcxxabi" linkflags="-stdlib=libc++" \
           -j$NPROCS install \
           --layout=tagged \
           --prefix=$INSTALL_DIR\
           --user-config=user-config.jam \
           variant=release \
           threading=multi \
           link=shared
else
    echo 'Using cached directory $basedir/software/install/boost-${VERSION}';
    echo "Cached version: ${VERSION}"
fi
