#! /bin/bash

. ~/.bash_profile

prefix=$COMMONDIR/local/boost/1.56-bgq-static
cat > user-config.jam << EOF
using mpi : mpicxx : :  ;
using gcc : 4.7.2 :  powerpc64-bgq-linux-g++ ;  
EOF

export CXX=powerpc64-bgq-linux-g++
export CC=powerpc64-bgq-linux-gcc-4.7.2
#./b2 clean
#./b2 toolset=clang-bg cxxflags="-stdlib=libc++ -std=c++11" linkflags="-stdlib=libc++"
#./b2 install -j20 -d2 --prefix=$prefix toolset=clang-bg cxxflags="-stdlib=libc++ -std=c++11" linkflags="-stdlib=libc++" \
./b2 install -j20 -d2 --prefix=$prefix  cxxflags=" -fPIC -Bdynamic -std=c++0x"  linkflags="-fPIC -dynamic -Bdynamic -std=c++0x "\
                --layout=tagged \
                --debug-configuration \
                --without-python \
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
