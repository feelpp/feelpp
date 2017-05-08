#!/bin/bash

# On Ubuntu 14.04.2:
# apt-get install libopenmpi-dev libpapi-dev opari

VERSION=1.4.2

if [[ -f "./scorep-${VERSION}.tar.gz" ]]; then 
    echo "Score-p source code already downloaded"
fi
if [[ -d "./scorep-${VERSION}" ]]; then 
    echo "Score-p seems to be already decompressed/configured"
fi
if [[ -d "/data/software/install/scorep-${VERSION}" ]]; then 
    echo "Score-p seems to be already installed"
    exit 1
fi

wget http://www.vi-hps.org/upload/packages/scorep/scorep-${VERSION}.tar.gz

tar zxvf scorep-${VERSION}.tar.gz

cd scorep-${VERSION}/
#CXXFLAGS="-fpermissive" ./configure --prefix=/data/software/install/scorep-${VERSION}/gcc-4.9.0/openmpi-1.8.5 --enable-shared
# permissive flag not required
./configure --prefix=/data/software/install/scorep-${VERSION}/gcc-4.9.0/openmpi-1.8.5 --enable-shared
# Directly do a make install, some intermediate libraries are installed through the process and needed for the following steps
make install 
