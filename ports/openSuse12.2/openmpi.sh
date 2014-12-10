#installation d'openmPI
mkdir $workdir/_openmpi
cd $workdir/_openmpi
if [ ! -d "openmpi-1.8.3" ]; then
wget -c http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.3.tar.bz2
tar xjf openmpi-1.8.3.tar.bz2
fi
cd openmpi-1.8.3
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$openmpiDir
make -j$nbProcs all
make install
