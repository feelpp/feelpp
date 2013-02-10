# On clone GIT
git clone https://github.com/feelpp/feelpp.git
mkdir feelpp_build_dir

# On "allume" les dépots
sudo rm /etc/yum.repos.d/rpm*
sudo rm /etc/yum.repos.d/slc6-cernonly*
sudo sed -i "s/enabled=0/enabled=1/g" -e /etc/yum.repos.d
sudo yum upgrade

# Install cd cmake28
sudo yum install cmake28
sudo ln -s /usr/bin/cmake28 /usr/bin/cmake

# Installation de gcc 4.7
mkdir /tmp/gcc
cd /tmp/gcc
#wget http://www.mpfr.org/mpfr-current/mpfr-3.1.1.tar.xz
#wget ftp://ftp.gmplib.org/pub/gmp-5.1.0/gmp-5.1.0a.tar.bz2
#wget http://www.multiprecision.org/mpc/download/mpc-1.0.1.tar.gz
wget ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-4.7.2/gcc-4.7.2.tar.bz2
tar xjf gcc-4.7.2.tar.bz2
#tar xjf gmp-5.1.0a.tar.bz2
#mv gmp-5.1.0 gcc-4.7.2/gmp
#tar xzf mpc-1.0.1.tar.gz
#mv mpc-1.0.1 gcc-4.7.2/mpc 
#unxz mpfr-3.1.1.tar.xz
#tar xf mpfr-3.1.1.tar
#mv mpfr-3.1.1 gcc-4.7.2/mpfr
mkdir BUILD_DIR
cd BUILD_DIR
../gcc-4.7.2/configure
make -j 
sudo make install
export LD_LIBRARY_PATH=/usr/local/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib/openmpi:$LD_LIBRARY_PATH
#installation d'openmPI
wget http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.bz2
tar xjf openmpi-1.6.3.tar.bz2
cd openmpi-1.6.3
./configure
./configure CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64
su -c "make -j all install"

export  PETSC_ARCH="arch-linux2-c-opt"
export  PETSC_DIR="$HOME/petsc-3.3-p5"
export GMSH_DIR="$HOME/gmsh-2.6.2-svn-Linux/"

#PETSC - Il y a des gags de droits, d'ou les sudo
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p5.tar.gz
tar xzf petsc-3.3-p5.tar.gz
cd petsc-3.3-p5
./configure \
  --with-shared-libraries=1 \
  --with-fortran=0 \
  --with-mpi-include=$MPI_INCLUDE \
  --with-mpi-lib=[$MPI_LIB/libmpi_cxx.so,$MPI_LIB/libmpi.so] \
  --with-blas-lib=$LIB_DIR/libblas.so \
  --with-lapack-lib=$LIB_DIR/lib/liblapack.so \
	--download-umfpack=1 \
        --download-ml \
        --download-metis \
        --download-parmetis \
        --download-blacs \
        --download-scalapack \
        --download-mumps \
        --download-pastix \
        --download-ptscotch \
  #--with-debugging=0 COPTFLAGS='-O3 -march=p4 -mtune=p4' FOPTFLAGS='-O3 -qarch=p4 -qtune=p4'
make all
make test

wget http://geuz.org/gmsh/bin/Linux/gmsh-svn-Linux64.tgz
tar xzf gmsh-svn-Linux64.tgz


