#!/bin/sh

module load cmake/2.8.7   
module load gcc/4.6.3  
module load openmpi/gcc/1.4.5  

export prefix=/usr
export MPI_ROOT=/softs/openmpi/gcc/1.4.5/
export BOOST_DIR=/softs/cemracs/boost/1.49.0

cat > user-config.jam << EOF
using mpi : $MPI_ROOT/bin/mpic++ : : $MPI_ROOT/bin/mpirun ;
EOF
./bootstrap.sh
./bjam install -j20 -d2 --prefix=$BOOST_DIR \
                --layout=tagged \
                --debug-configuration \
                --user-config=user-config.jam \
                -sBZIP2_INCLUDE=${prefix}/include \
                -sBZIP2_LIBPATH=${prefix}/lib \
                -sEXPAT_INCLUDE=${prefix}/include \
                -sEXPAT_LIBPATH=${prefix}/lib \
                -sZLIB_INCLUDE=${prefix}/include \
                -sZLIB_LIBPATH=${prefix}/lib \
                -sICU_PATH=${prefix} \
                variant=release \
                threading=single,multi \
                link=static,shared
