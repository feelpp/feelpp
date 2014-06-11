#! /bin/bash -x

# Optionnaly set version and nbProc
: ${Version:="3.1.2"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de mpfr-$Version"

# Should add this in your $HOME/.bashrc
#if [ -f $HOME/github/feelpp.config/etc/feelpprc.sh ]; then
#   . $HOME/github/feelpp.config/etc/feelpprc.sh
#fi

module load libs/gmp-6.0.0_gcc_${gccVersion}.feelpp
export LD_PRELOAD=/applis/ciment/v2/x86_64/lib/libc.so.6

#Création d'un repertoire temporaire
  tar xjf ~/Downloads/mpfr-$Version.tar.bz2
  cd mpfr-$Version
  mkdir -p build
  cd  build
  #installation en mode utilisateur
  # see http://www.linuxfromscratch.org/lfs/view/development/chapter06/gcc.html

  ../configure --prefix=/home/trophime/modules/libs/mpfr-${Version}_gcc_${gccVersion} \
     --with-gmp-include=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion}/include \
     --with-gmp-lib=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion}/lib
  make -j$nbProc
  make check
  make -j$nbProc install
