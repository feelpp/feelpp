#!/bin/bash -x

installDir=${1}
patchDir=${2}
nbProc=${3}
CxxCompiler=${4}

## should be able to choice a version or svn
#VERSION="svn"
VERSION="2.16.0"
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

if [ -z ${PETSC_DIR} ]; then
  export PETSC_DIR=$HOME/modules/libs/petsc-3.6.3.dfsg1_mpicc_4.9.3
  export PETSC_ARCH=arch-linux2-mpicc-opt
fi
if [ -z ${SLEPC_DIR} ]; then
  export SLEPC_DIR=$HOME/modules/libs/slepc-3.6.2.dfsg_mpicc_4.9.3
fi

echo "PETSC_DIR=$PETSC_DIR"
echo "SLEPC_DIR=$SLEPC_DIR"

echo "Install to ${installDir}"
cmake .. \
  -DCMAKE_INSTALL_PREFIX=${installDir} \
  -DCMAKE_CXX_COMPILER=${CxxCompiler} \
  -DCMAKE_BUILD_TYPE=release \
  -DENABLE_MPI=ON \
  -DENABLE_BUILD_LIB=1 \
  -DENABLE_BUILD_SHARED=1 \
  -DENABLE_KBIPACK:BOOL=ON \
  -DENABLE_METIS:BOOL=ON \
  -DENABLE_TAUCS:BOOL=OFF \
  -DENABLE_PETSC:BOOL=ON \
  -DENABLE_SLEPC:BOOL=ON

make -j${nbProc} lib
make -j${nbProc} shared
make -j${nbProc} install

