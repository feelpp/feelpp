export LD_LIBRARY_PATH=/data/software/install/gcc-4.9.0:/data/software/install/gcc-4.9.0/lib:/data/software/install/gcc-4.9.0/lib32:/data/software/install/gcc-4.9.0/lib64:$LD_LIBRARY_PATH
export CC=/data/software/install/gcc-4.9.0/bin/gcc
export CXX=/data/software/install/gcc-4.9.0/bin/g++
../configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=/data/software/install/openmpi-1.8.5/gcc-4.9.0/
make -j all
make -j install
