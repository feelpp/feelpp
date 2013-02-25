#installation d'openmPI
mkdir $workdir/_openmpi
mkdir $workdir/_openmpi/gcc_$gccVersion
cd $workdir/_openmpi/gcc_$gccVersion
if [ ! -d "openmpi-1.6.3" ]; then
wget -c http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.bz2
tar xjf openmpi-1.6.3.tar.bz2
fi
cd openmpi-1.6.3
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$openmpiDir
make -j$nbProcs all
make install
