#! /bin/bash

VERSION=${1:-2.4.2}
PREFIX=${2:-/opt/software/install/singularity/$VERSION}
echo "Install singularity $VERSION in $PREFIX"
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
tar xvf singularity-$VERSION.tar.gz
cd singularity-$VERSION
./configure --prefix=$PREFIX
make
sudo make install
