#!/bin/sh

VERSION=$1
basedir=$HOME
#basedir=/tmp

set -e
# check to see if protobuf folder is empty
if [ ! -d "$basedir/software/install/petsc-${VERSION}" ]; then

    wget -c http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${VERSION}.tar.gz

    tar xzf petsc-lite-${VERSION}.tar.gz
    cd petsc-${VERSION}
    
    ./configure --with-shared-libraries=1 \
        --with-debugging=0 \
        COPTFLAGS='-O3' FOPTFLAGS='-O3' \
        --prefix=$HOME/software/install/petsc-${VERSION} \
        --download-suitesparse=1 \
        --download-ml \
        --download-metis \
        --download-parmetis \
        --download-blacs \
        --download-scalapack \
        --download-fblaslapack \
        --download-mumps \
        --download-pastix \
        --download-ptscotch
make all
#make test
make install

else
    echo 'Using cached directory $basedir/software/install/petsc-${VERSION}';
    echo "Cached version: ${VERSION}";
fi
