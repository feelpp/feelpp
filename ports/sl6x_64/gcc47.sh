# Installation de gcc 4.7
#Création d'un repertoire temporaire
mkdir $workdir/_gcc
cd $workdir/_gcc
wget ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-$gccVersion/gcc-$gccVersion.tar.bz2
tar xjf gcc-$gccVersion.tar.bz2
mkdir BUILD_DIR
cd BUILD_DIR
#installation en mode utilisateur
../gcc-$gccVersion/configure --prefix=$gccDir
make -j4
make -j4 install
export PATH=$gccDir/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$gccDir/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$gccDir/lib64
export CC=$gccDir/bin/gcc
export CXX=$gccDir/bin/g++
export FC=$gccDir/bin/gfortran
export F90=$gccDir/bin/gfortran
