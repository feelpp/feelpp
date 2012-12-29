#!/bin/sh

# LAPACK
# tested on Fri, Dec 14, 2012

# file lapack-3.4.1.tgz
# size 6147915
# md5  44c3869c38c8335c2b9c2a8bb276eb55
# sha1 c115223ac1bac9ab971aae865d3e95442bc979bc

export PREFIX=$CCCWORKDIR/local-titane

module unload intel/11.1.056
module unload bullxmpi/1.0.2
module load cmake/2.8.0
module load gcc/4.6.3
module load openmpi/1.4.2_gnu

mkdir build
cd build

sed 's/^\(. *WORKING_DIRECTORY\)/#\1/' ../TESTING/CMakeLists.txt > tmp.txt
mv tmp.txt ../TESTING/CMakeLists.txt

cmake \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_Fortran_COMPILER=/applications/gcc-4.6.3/bin/gfortran \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  ..

make -j8 all install
