#! /bin/bash
# Optionnaly set version and nbProc
: ${Version:="1.1"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de ppl-$Version"

#if [ -f $HOME/github/feelpp.config/etc/feelpprc.sh ]; then
#   . $HOME/github/feelpp.config/etc/feelpprc.sh
#fi

module load libs/gmp-6.0.0_gcc_${gccVersion}.feelpp
export LD_PRELOAD=/applis/ciment/v2/x86_64/lib/libc.so.6

tar jxvf ~/Downloads/ppl-${VERSION}.tar.bz2
cd ppl-${VERSION}
mkdir -p build
cd build/
../configure --prefix=/home/trophime/modules/libs/ppl-${Version}_gcc_${gccVersion} \
	--with-gmp=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion}
make -j${nbProc}
make -j${nbProc} install
#make check
