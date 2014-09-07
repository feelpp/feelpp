 #! /bin/bash
 # Optionnaly set version and nbProc
: ${Version:="6.0.0"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de gmp-$Version"

# Should add this in your $HOME/.bashrc
#if [ -f $HOME/github/feelpp.config/etc/feelpprc.sh ]; then
#   . $HOME/github/feelpp.config/etc/feelpprc.sh
#fi

export LD_PRELOAD=/applis/ciment/v2/x86_64/lib/libc.so.6

tar jxvf ~/Downloads/gmp-${Version}a.tar.bz2 

 cd gmp-${Version}
 mkdir -p build
 cd build/
 ../configure --prefix=/home/trophime/modules/libs/gmp-${Version}_gcc_${gccVersion} --enable-cxx
 make -j$nbProc
 make -j$nbProc install
 make check
