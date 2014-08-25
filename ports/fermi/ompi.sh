#!/ bin/sh

. ~/.bash_profile
module load fermi/bgqgcc472.feelpp
#module load compilers/bgclang.feelpp

../configure --prefix=$WORK/local/openmpi/1.8.1-bggcc472  --host=powerpc64-bgq-linux  --enable-static --disable-shared CPPFLAGS="-I/bgsys/drivers/ppcfloor -I/bgsys/drivers/ppcfloor/spi/include/kernel/cnk/" LDFLAGS="-dynamic -L/bgsys/drivers/ppcfloor/cnk/lib/"




make -j12 -k;make -k install

#--disable-wrapper-rpath
