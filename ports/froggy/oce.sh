#!/bin/bash

#profile should include feelpprc from V.C
. ~/.bash_profile

VERSION=0.12

module load gcc/4.6.2_gcc-4.6.2
module load libX11/1.4.0_gcc-4.6.2
module load libXmu/1.1.0_gcc-4.6.2
module load libXext/1.2.0_gcc-4.6.2
module load Mesa/8.0.3_gcc-4.6.2
module load tcl/8.5.11_gcc-4.6.2
module load tk/8.5.11_gcc-4.6.2

export LD_LIBRARY_PATH=$gccDir/libX11_1.4.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$gccDir/libXmu_1.1.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$gccDir/libXext_1.2.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$gccDir/Mesa_8.0.3/lib:$LD_LIBRARY_PATH

X11_INC_SEARCH_PATH
X11_LIB_SEARCH_PATH

export CFLAGS="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export CXXFLAGS="-g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security"
export FFLAGS="-g -O2"
export LDFLAGS="-Wl,-z,relro -Wl,--as-needed"

# 
if [ ! -d "${HOME}/packages/build" ]; then
  mkdir -p ${HOME}/packages/build
fi
pushd ${HOME}/packages/build

#apply patches if any
#if [ -d "${HOME}/packages/patches/oce" ]; then
# pushd ${HOME}/packages/oce-OCE-$VERSION
# for patch in $( ls ${HOME}/packages/patches/oce/*.patch ); do
#      patch -p1 < $patch || exit
# done
# popd
#fi

cmake \
	${HOME}/packages/oce-OCE-$VERSION \
        -DOCE_BUILD_SHARED_LIB:BOOL=ON \
        -DOCE_TESTING:BOOL=ON \
        -DOCE_BUILD_TYPE:STRING=RelWithDebInfo \
        -DOCE_INSTALL_PREFIX:PATH=$HOME/packages/oce \
        -DOCE_INSTALL_LIB_DIR:PATH=lib \
        -DOCE_INSTALL_CMAKE_DATA_DIR:PATH=share \
        -DOCE_DRAW:BOOL=OFF \
        -DOCE_RPATH_FILTER_SYSTEM_PATHS:BOOL=ON \
        -DCMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES:STRING='/lib;/usr/lib' \
        -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING='$(CFLAGS)' \
        -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING='$(CXXFLAGS)' \
        -DCMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO:STRING='$(LDFLAGS)' \
        -DCMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO:STRING='$(LDFLAGS)' \
        -DOCE_MULTITHREAD_LIBRARY:STRING=NONE \
        -DOCE_WITH_FREEIMAGE:BOOL=ON \
        -DOCE_WITH_GL2PS:BOOL=OFF

make -j 4
make install
