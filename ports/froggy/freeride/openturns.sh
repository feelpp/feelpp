#! /bin/bash

#profile should include feelpprc from V.C
. ~/.bash_profile

# to get openturns
module load R/2.14.1_gcc-4.6.2
module load swig/2.0.9_gcc-4.6.2
module load libxml2/2.9.0_gcc-4.6.2
module load python/2.7.3_gcc-4.6.2

module unload bison/2.6.5_gcc-4.6.2 
export PATH=$PATH:${HOME}/packages/bison/bin

# optional packages numpy, scipy, matplotlib
if [ ! -d "${HOME}/packages/build" ]; then
  mkdir -p ${HOME}/packages/build
fi
pushd ${HOME}/packages/build

#apply patches if any
if [ -d "${HOME}/packages/patches/openturns-1.1" ]; then
 pushd ${HOME}/packages/openturns-1.1
 for patch in $( ls ${HOME}/packages/patches/openturns-1.1/*.patch ); do
      patch -p1 < $patch || exit
 done
 popd
fi

#R_LIBS_USER is set to directory "R/R.version$platform-library/x.y" for R x.y.z.
export R_LIBS_USER=${HOME}/packages/R/lib
export PYTHONPATH=${HOME}/packages/lib/python2.7/site-packages/

#
export cflags="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export cxxflags="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export fflags="-g -O2"
export ldflags="-Wl,-z,relro -Wl,--as-needed"

# Pre-install To allow OpenTURNS to communicate with R
R CMD INSTALL R --library=${R_LIBS_USER} ${HOME}/packages/openturns-1.1/utils/rot_1.4.5.tar.gz

cmake \
	${HOME}/packages/openturns-1.1 \
            -DUSE_TBB:BOOL=OFF \
	    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
            -DCMAKE_CXX_COMPILER:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/gcc_4.6.2/bin/g++ \
            -DCMAKE_Fortran_COMPILER:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/gcc_4.6.2/bin/gfortran \
            -DCMAKE_C_COMPILER:FILEPATH=/applis/ciment/v2/stow/x86_64/gcc_4.6.2/gcc_4.6.2/bin/gcc \
            -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
	    -DCMAKE_INSTALL_PREFIX:PATH=${HOME}/packages/openturns \
	    -DCMAKE_INSTALL_LIBDIR:PATH=${HOME}/packages/openturns/lib \
	    -DOPENTURNS_LIBRARY_PATH:PATH=lib \
            -DOPENTURNS_CONFIG_CMAKE_PATH:PATH=lib/cmake/openturns \
            -DOPENTURNS_WRAPPER_PATH:PATH=lib/openturns-1.1/wrappers \
            -DOPENTURNS_SYSTEM_INSTALL:BOOL=ON \
	    -DOPENTURNS_SYSCONFIG_PATH:PATH=${HOME}/packages/openturns/etc \
	    -DOPENTURNS_DATA_PATH:PATH=${HOME}/packages/openturns/share \
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
