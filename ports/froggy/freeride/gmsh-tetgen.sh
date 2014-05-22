#! /bin/bash

#profile should include feelpprc from V.C
. ~/.bash_profile

VERSION=2.7.1~svn15108

# to get gmsh
module load openmpi/1.4.4_gcc-4.6.2
module load swig/2.0.9_gcc-4.6.2
module load python/2.7.3_gcc-4.6.2
module load Mesa/8.0.3_gcc-4.6.2

module unload bison/2.6.5_gcc-4.6.2 
export PATH=$PATH:${HOME}/packages/bison/bin

# optional packages numpy, scipy, matplotlib
if [ ! -d "${HOME}/packages/build" ]; then
  mkdir -p ${HOME}/packages/build
fi
pushd ${HOME}/packages/build

#apply patches if any
if [ -d "${HOME}/packages/patches/gmsh-tetgen" ]; then
 pushd ${HOME}/packages/gmsh-tetgen-$VERSION
 for patch in $( ls ${HOME}/packages/patches/gmsh-tetgen/*.patch ); do
      patch -p1 < $patch || exit
 done
 popd
fi

#
export cflags="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export cxxflags="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export fflags="-g -O2"
export ldflags="-Wl,-z,relro -Wl,--as-needed"

cmake \
	${HOME}/packages/gmsh-tetgen-$VERSION \
-DENABLE_KBIPACK:BOOL=ON \
-DENABLE_METIS:BOOL=ON \
-DENABLE_TAUCS:BOOL=OFF \
-DENABLE_GRAPHICS:BOOL=OFF \
-DENABLE_MPI:BOOL=ON \
-DENABLE_WRAP_PYTHON:BOOL=ON \
-DENABLE_BUILD_SHARED:BOOL=ON \
-DCMAKE_INCLUDE_PATH:STRING="/applis/ciment/v2/stow/x86_64/gcc_4.6.2/openmpi_1.4.4/include" \
	    -DCMAKE_INSTALL_PREFIX:PATH=${HOME}/packages/gmsh-tetgen \
	    -DINSTALL_LIB_DIR:PATH=${HOME}/packages/gmsh-tetgen/lib \
	    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
            -DCMAKE_CXX_COMPILER:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/gcc_4.6.2/bin/g++ \
	-DCMAKE_CXX_FLAGS="-DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -fopenmp -lmpi -fPIC -Wall $cxxflags" \
            -DCMAKE_Fortran_COMPILER:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/gcc_4.6.2/bin/gfortran \
            -DCMAKE_C_COMPILER:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/gcc_4.6.2/bin/gcc \
	    -DCMAKE_INSTALL_PREFIX:PATH=${HOME}/packages/gmsh-tetgen \
	    -DCMAKE_INSTALL_LIBDIR:PATH=${HOME}/packages/gmsh-tetgen/lib \
	    -DBLAS_LIBRARIES:FILEPATH=/home/chabannes/packages/blas/BLAS/blas_LINUX.a \
            -DPYTHON_EXECUTABLE:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/python_2.7.3/bin/python \
            -DPYTHON_INCLUDE_DIR:PATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/python_2.7.3/include/python2.7 \
            -DPYTHON_LIBRARY:PATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/python_2.7.3/lib/libpython2.7.a \
 || exit


make -j 4 || exit
make install  || exit

popd
rm -rf ${HOME}/packages/build

#make all lib shared
