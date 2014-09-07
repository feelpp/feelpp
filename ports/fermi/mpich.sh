#!/ bin/sh

. ~/.bash_profile
module load fermi/bgqgcc472.feelpp
module load compilers/bgclang.feelpp

#CC=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-gcc \
#CXX=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-g++ \
#F77=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-gfortran \
#FC=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-gfortran \
unset FC F77

AR=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-ar \
LD=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-ld \
RANLIB=$WORK/bgq/gnu-linux-4.7.2/bin/powerpc64-bgq-linux-ranlib \
CPPFLAGS="-I/bgsys/drivers/ppcfloor -I/bgsys/drivers/ppcfloor/spi/include/kernel/cnk" \
LDFLAGS="-dynamic -L/bgsys/drivers/ppcfloor/cnk/lib/ -Wl,--allow-multiple-definition"  \
../configure --prefix=$WORK/local/mpich/3.1.1-bgclang  --host=powerpc64-bgq-linux --target=powerpc64-bgq-linux --build=powerpc64-linux-gnu \
           --with-pm=none --with-device=pamid --with-file-system=gpfs:BGQ --with-file-system=bg+bglockless \
           --enable-timing=no --disable-collchk --disable-graphics --disable-rlog --disable-sample --disable-rpath \
           --with-aint-size=8 --with-assert-level=2 --enable-fast=O3 --enable-error-messages --disable-debuginfo \
           --enable-thread-cs=per-object --enable-handle-allocation=tls --enable-refcount=lock-free --disable-predefined-refcount \
           --with-bgq-install-dir=$WORK/bgq/gnu-linux-4.7.2  --enable-shared --disable-static --disable-fortran LDFLAGS="-dynamic -L/bgsys/drivers/ppcfloor/cnk/lib/ -Wl,--allow-multiple-definition"  MPI_SIZEOF_OFFSET="8" MPI_OFFSET_TYPE="long long"




make -j12 -k;make -k install

#--disable-wrapper-rpath
