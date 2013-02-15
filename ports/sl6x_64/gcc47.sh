# Installation de gcc 4.7
mkdir $workdir/_gcc
cd $workdir/_gcc
wget ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-4.7.2/gcc-4.7.2.tar.bz2
tar xjf gcc-4.7.2.tar.bz2
mkdir BUILD_DIR
cd BUILD_DIR
#installation en mode utilisateur
../gcc-4.7.2/configure --prefix=$workdir/gcc
make -j4
make -j4 install
export LD_LIBRARY_PATH=$workdir/gcc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vhuber/work/gcc/lib/../lib64
export CC=$workdir/gcc/bin/gcc
export CXX=$workdir/gcc/bin/g++
export FC=$workdir/gcc/bin/gfortran
export F90=$workdir/gcc/bin/gfortran
