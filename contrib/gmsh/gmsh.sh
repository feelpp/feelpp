#!/bin/bash -x

installDir=${1}
patchDir=${2}
nbProc=${3}
CxxCompiler=${4}

## should be able to choice a version or svn
VERSION="svn"
## VERSION="2.12.0"
wget -c http://geuz.org/gmsh/src/gmsh-$VERSION-source.tgz
tar xzf gmsh-$VERSION-source.tgz

cd `tar -tvf gmsh-$VERSION-source.tgz | head -n 1 | sed "s/\// /g" | awk '{print $7}'`
echo "$PWD"

echo "Apply patches from ${patchDir}"
for file in $(ls -1 ${patchDir}); do
  echo "Apply patch $file";
  patch -p1 <   ${patchDir}/$file;
done
mkdir -p build
cd build
echo "$PWD"

echo "Install to ${installDir}"
cmake .. \
  -DCMAKE_INSTALL_PREFIX=${installDir} \
  -DCMAKE_CXX_COMPILER=${CxxCompiler} \
  -DCMAKE_BUILD_TYPE=release \
  -DENABLE_MPI=ON \
  -DENABLE_BUILD_LIB=1 \
  -DENABLE_BUILD_SHARED=1 \
  -DENABLE_SLEPC=0

make -j${nbProc} lib
make -j${nbProc} shared
make -j${nbProc} install

