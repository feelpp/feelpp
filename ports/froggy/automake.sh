#! /bin/bash -x

# Optionnaly set version and nbProc
: ${Version:="1.14"}
: ${gccVersion:="4.4.6"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de automake-$Version"

#if [ -f $HOME/github/feelpp.config/etc/feelpprc.sh ]; then
#   . $HOME/github/feelpp.config/etc/feelpprc.sh
#fi

#export LD_PRELOAD=/applis/ciment/v2/x86_64/lib/libc.so.6

#Création d'un repertoire temporaire
  tar xzf ~/Downloads/automake-$Version.tar.gz
  cd automake-$Version
  mkdir -p build
  cd  build
  #installation en mode utilisateur
  # see http://www.linuxfromscratch.org/lfs/view/development/chapter06/gcc.html

  ../configure --prefix=/home/trophime/modules/tools/automake-${Version}_gcc_${gccVersion}
  make -j$nbProc
#  make check
  make -j$nbProc install
