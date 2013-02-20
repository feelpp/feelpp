#installation d'openmPI
mkdir $workdir/_openmpi
cd $workdir/_openmpi
wget -c http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.bz2
tar xjf openmpi-1.6.3.tar.bz2
cd openmpi-1.6.3
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$openmpiDir
make -j$nbProcs all install
