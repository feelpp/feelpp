#! /bin/bash -x

# Optionnaly set version and nbProc
: ${Version:="1.0.2"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de mpc-$Version"

#if [ -f $HOME/github/feelpp.config/etc/feelpprc.sh ]; then
#   . $HOME/github/feelpp.config/etc/feelpprc.sh
#fi

module load libs/gmp-6.0.0_gcc_${gccVersion}.feelpp
module load libs/mpfr-3.1.2_gcc_${gccVersion}.feelpp
export LD_PRELOAD=/applis/ciment/v2/x86_64/lib/libc.so.6

#Création d'un repertoire temporaire
  tar xzf ~/Downloads/mpc-$Version.tar.gz
  cd mpc-$Version
  mkdir -p build
  cd  build
  #installation en mode utilisateur
  # see http://www.linuxfromscratch.org/lfs/view/development/chapter06/gcc.html

  ../configure --prefix=/home/trophime/modules/libs/mpc-${Version}_gcc_${gccVersion} \
     --with-mpfr-include=/home/trophime/modules/libs/mpfr-3.1.2_gcc_${gccVersion}/include \
     --with-mpfr-lib=/home/trophime/modules/libs/mpfr-3.1.2_gcc_${gccVersion}/lib \
     --with-gmp-include=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion}/include \
     --with-gmp-lib=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion}/lib
  make -j$nbProc
  make check
  make -j$nbProc install
