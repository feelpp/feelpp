#!/bin/bash -x

installDir=${1}
patchDir=${2}
nbProc=${3}


wget -c http://geuz.org/gmsh/src/gmsh-svn-source.tgz
tar xzf gmsh-svn-source.tgz

cd `tar -tvf gmsh-svn-source.tgz | head -n 1 | sed "s/\// /g" | awk '{print $7}'`

echo "Apply patches from ${patchDir}"
for file in $(ls -1 ${patchDir}); do
  echo "Apply patch $file";
  patch -p1 <  $file;
done
mkdir build
cd build

cmake .. \
  -DCMAKE_INSTALL_PREFIX=${installDir} \
  -DCMAKE_BUILD_TYPE=release \
  -DENABLE_MPI=ON \
  -DENABLE_BUILD_LIB=1 \
  -DENABLE_BUILD_SHARED=1 \
  -DENABLE_SLEPC=0

make -j${nbProc} lib
make -j${nbProc} shared
make -j${nbProc} install

