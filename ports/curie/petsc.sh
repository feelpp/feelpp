#! /bin/bash

. ~/.bash_profile
module load c++/gnu/4.5.1   
prefix=$WORKDIR/local-gcc45

./configure --prefix=$prefix --with-cc=/usr/local/gcc-4.5.1/bin/gcc --with-cxx=/usr/local/gcc-4.5.1/bin/g++ --with-mpi-dir=$MPI_ROOT
