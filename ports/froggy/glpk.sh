#! /bin/bash

# Optionnaly set version and nbProc
: ${Version:="4.54"}
: ${gccVersion:="4.8.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de glpk-$Version"

 
tar zxf ~/Downloads/glpk_${Version}.orig.tar.gz 

cd glpk-${Version}/
mkdir -p build
cd build/
../configure --prefix=/home/trophime/modules/libs/glpk-${Version}_gcc_${gccVersion}
make -j ${nbProc} 
make install
 
