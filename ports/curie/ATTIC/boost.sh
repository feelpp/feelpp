#! /bin/bash

. ~/.bash_profile
module load c++/gnu/4.6.3   
prefix=$WORKDIR/local/gcc46
cat > user-config.jam << EOF
using mpi : $MPI_ROOT/bin/mpic++ : : $MPI_ROOT/bin/mpirun ;
EOF
./bjam install -j20 -d2 --prefix=$prefix \
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
