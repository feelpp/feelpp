#!/bin/sh

VERSION=1_59_0
basedir=$HOME
#basedir=/tmp

set -e
# check to see if protobuf folder is empty
if [ ! -d "$basedir/software/install/boost" ]; then

wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.bz2/download -O boost_1_59_0.tar.bz2

tar xjf boost_1_59_0.tar.bz2
cd boost_1_59_0

echo "using mpi ;" >> user-config.jam
echo "" >> user-config.jam
./bootstrap.sh
./bjam -j$NPROCS install \
      --layout=tagged \
      --prefix=$basedir/software/install/boost \
      --user-config=user-config.jam \
      variant=release \
      threading=multi \
      link=shared
else
    echo 'Using cached directory $basedir/software/install/boost';
    echo "Cached version: ${VERSION}"
fi
