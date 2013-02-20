# Installation de gcc_$gccVersion
#Création d'un repertoire temporaire
  mkdir $workdir/_gcc/
  mkdir $workdir/_gcc/gcc_$gccVersion
  cd $workdir/_gcc/gcc_$gccVersion
  if [ ! -d "gcc-$gccVersion" ]; then
    wget -c ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-$gccVersion/gcc-$gccVersion.tar.bz2
    tar xjf gcc-$gccVersion.tar.bz2
  fi
  mkdir BUILD_DIR_$gccVersion
  cd BUILD_DIR_$gccVersion
  #installation en mode utilisateur
  ../gcc-$gccVersion/configure --prefix=$gccDir
  make -j$nbProc
  make -j$nbProc install
