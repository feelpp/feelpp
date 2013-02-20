# Installation de gcc 4.7
#Création d'un repertoire temporaire
if [ ! -d "gcc-$gccdir" ]; then
  mkdir $workdir/_gcc
  cd $workdir/_gcc
  if [ ! -d "gcc-$gccVersion" ]; then
    wget -c ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-$gccVersion/gcc-$gccVersion.tar.bz2
    tar xjf gcc-$gccVersion.tar.bz2
  fi
  mkdir BUILD_DIR
  cd BUILD_DIR
  #installation en mode utilisateur
  ../gcc-$gccVersion/configure --prefix=$gccDir
  make -j$nbProc
  make -j$nbProc install
fi
