#! /bin/bash

# BOOST
# tested on Fri, Dec 14, 2012

# file boost_1_49_0.tar.gz
# size 59136116
# md5  e0defc8c818e4f1c5bbb29d0292b76ca
# sha1 fce1a7d8d9866e39e4797975b5132d33ad919511

export PREFIX=$CCCWORKDIR/local-titane
export MPI_ROOT=/applications/openmpi-1.4.2_gnu

module unload intel/11.1.056
module unload bullxmpi/1.0.2
module load cmake/2.8.0
module load gcc/4.6.3
module load openmpi/1.4.2_gnu

echo "using mpi : $MPI_ROOT/bin/mpic++ : : $MPI_ROOT/bin/mpirun ;" > user-config.jam
./bootstrap.sh
./b2 install \
  -j8 -d2 --prefix=$PREFIX \
  --layout=tagged \
  --debug-configuration \
  --user-config=user-config.jam \
  -sBZIP2_INCLUDE=$PREFIX/include \
  -sBZIP2_LIBPATH=$PREFIX/lib \
  -sEXPAT_INCLUDE=$PREFIX/include \
  -sEXPAT_LIBPATH=$PREFIX/lib \
  -sZLIB_INCLUDE=$PREFIX/include \
  -sZLIB_LIBPATH=$PREFIX/lib \
  -sICU_PATH=$PREFIX \
  variant=release \
  threading=single,multi \
  link=static,shared
