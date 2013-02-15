#installation d'openmPI
mkdir $workdir/_openmpi
cd $workdir/_openmpi
wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.bz2
tar xjf openmpi-1.6.3.tar.bz2
cd openmpi-1.6.3
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64 --prefix=$workdir/openmpi
make -j all install
export PATH=$workdir/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$workdir/openmpi/libmake:$LD_LIBRARY_PATH
