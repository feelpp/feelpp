#! /bin/bash

. ~/.bash_profile
export PATH=$WORK/bgq/gnu-linux-4.7.2/bin:$PATH
#module load compilers/bgq-gcc47.feelpp
prefix=$WORK/local/boost/1.54-bgq
cat > user-config.jam << EOF
using mpi : mpic++ : :  ;
EOF
export CXX=powerpc64-bgq-linux-g++
export CC=powerpc64-bgq-linux-gcc
./b2 clean
#./b2 toolset=clang-bg cxxflags="-stdlib=libc++ -std=c++11" linkflags="-stdlib=libc++"
#./b2 install -j20 -d2 --prefix=$prefix toolset=clang-bg cxxflags="-stdlib=libc++ -std=c++11" linkflags="-stdlib=libc++" \
./b2 install -j20 -d2 --prefix=$prefix toolset=gcc-bgq cxxflags="-std=c++11" \
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
                link=static
