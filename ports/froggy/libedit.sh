#! /bin/bash

# Optionnaly set version and nbProc
: ${Version:="3.1"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de libedit-$Version"

module load ncurses/5.9_gcc_${gccVersion}
#module load vim/7.3_gcc-4.6.2
 
tar zxvf ~/Downloads/libedit_3.1-20140213.orig.tar.gz 

cd libedit-20140213-3.1/
mkdir -p build
cd build/
../configure --prefix=/home/trophime/modules/libs/libedit-${Version}_gcc_${gccVersion} --enable-widec
make -j ${nbProc} 
make install
 
