#! /bin/bash
# Optionnaly set version and nbProc
: ${Version:="0.16.1"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de cloog-parma-$Version"

#if [ -f $HOME/github/feelpp.config/etc/feelpprc.sh ]; then
#   . $HOME/github/feelpp.config/etc/feelpprc.sh
#fi

module load libs/gmp-6.0.0_gcc_${gccVersion}.feelpp
module load libs/ppl-1.1_gcc_${gccVersion}.feelpp
export LD_PRELOAD=/applis/ciment/v2/x86_64/lib/libc.so.6

tar zxf ~/Downloads/cloog-parma-${Version}.tar.gz
cd cloog-parma-${Version}
mkdir -p build
cd build/
../configure --prefix=/home/trophime/modules/libs/cloog-parma-${Version}_gcc_${gccVersion} \
	--with-gmp-prefix=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion} \
	--with-ppl-prefix=/home/trophime/modules/libs/ppl-1.1_gcc_${gccVersion}
make -j${nbProc}
make -j${nbProc} install
#make check
