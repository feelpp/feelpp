#! /bin/bash -x

# Optionnaly set version and nbProc
: ${Version:="4.7.3"}
: ${gccVersion:="4.6.2"}
: ${nbProc:=$(echo `/usr/bin/getconf _NPROCESSORS_ONLN`)}
echo "Installation de gcc_$gccVersion"

#source ~/github/feelpp.config/etc/feelpprc.sh

module load libs/gmp-6.0.0_gcc_${gccVersion}.feelpp
module load libs/mpfr-3.1.2_gcc_${gccVersion}.feelpp
module load libs/mpc-1.0.2_gcc_${gccVersion}.feelpp
module load libs/ppl-1.1_gcc_${gccVersion}.feelpp
module load libs/cloog-parma-0.16.1_gcc_${gccVersion}.feelpp

module load zlib/1.2.7_gcc-${gccVersion}

module load tools/autoconf-2.64_gcc_${gccVersion}
module load automake/1.12_gcc-${gccVersion}
module load libtool/2.4_gcc-${gccVersion}

#Création d'un repertoire temporaire
tar xjf ~/Downloads/gcc-$Version.tar.bz2
cd gcc-$Version

# Eventually apply patches
for file in $(ls ~/Patches/gcc-${Version}*); do
   echo "Applying " $file
   patch -p1 < $file
done

autoreconf

#installation en mode utilisateur
# see http://www.linuxfromscratch.org/lfs/view/development/chapter06/gcc.html
mkdir -p build
cd  build

../configure --prefix=/home/trophime/modules/compilers/gcc-$Version \
  --enable-bootstrap \
  --enable-shared \
  --enable-threads=posix \
  --enable-checking=release \
  --with-system-zlib \
  --enable-__cxa_atexit \
  --disable-libunwind-exceptions \
  --enable-gnu-unique-object \
  --enable-languages=c,c++,fortran \
  --disable-dssi \
  --with-mpc=/home/trophime/modules/libs/mpc-1.0.2_gcc_${gccVersion} \
  --with-gmp=/home/trophime/modules/libs/gmp-6.0.0_gcc_${gccVersion} \
  --with-mpfr=/home/trophime/modules/libs/mpfr-3.1.2_gcc_${gccVersion} \
  --with-ppl=/home/trophime/modules/libs/ppl-1.1_gcc_${gccVersion} \
  --with-cloog=/home/trophime/modules/libs/cloog-parma-0.16.1_gcc_${gccVersion} \
  --with-tune=generic \
  --with-arch_32=i686 
  make -j$nbProc
  make -j$nbProc install
